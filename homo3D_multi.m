function CH = homo3D_multi(lx, ly, lz, lambda_vec, mu_vec, voxel,reg_coeff)
% homo3D_multi  三维周期复合材料均匀化（六面体单元，周期边界条件）
%
% 输入：
%   lx,ly,lz     - 单元胞物理尺寸
%   lambda_vec   - 各相拉梅常数 lambda（第k个值对应 voxel==k 的相）
%   mu_vec       - 各相拉梅常数 mu
%   voxel        - 相标记矩阵，nelx×nely×nelz，整数（0=空洞,1=相1,2=相2,...）
% 输出：
%   CH           - 6×6 均匀化弹性刚度矩阵

    [nelx, nely, nelz] = size(voxel);
    nel  = nelx * nely * nelz;
    ndof = 3 * nel;   % 周期DOF总数 = 3 × 周期节点数

    % ---- 单元刚度矩阵与等效载荷 ----
    [keLambda, keMu, feLambda, feMu] = hexahedron(lx/(2*nelx), ly/(2*nely), lz/(2*nelz));

    % ================================================================
    % edof 构造（核心修复）
    % ----------------------------------------------------------------
    % 原代码用 nodenrs（非周期，31x31x31）+ addrow 偏移量构造 edof，
    % 导致角点单元的 edof 最大值（89376）超过 dofVector 长度（89373），
    % 差值恰好3，即最后一个周期节点的3个DOF越界。
    %
    % 正确做法：直接用周期节点编号构造 edof。
    % 周期节点编号：node(ix,iy,iz) = mod(ix,nelx) + mod(iy,nely)*nelx
    %                                + mod(iz,nelz)*nelx*nely + 1
    % 共 nel 个周期节点，DOF 范围 [1, 3*nel]，无越界可能。
    % ================================================================

    % 各单元左下角节点的0-indexed坐标（列主序，与MATLAB reshape一致）
    [IX, IY, IZ] = ndgrid(0:nelx-1, 0:nely-1, 0:nelz-1);
    IX = IX(:);  IY = IY(:);  IZ = IZ(:);   % nel×1

    % 周期节点编号函数（1-indexed）
    pnode = @(dx,dy,dz) mod(IX+dx, nelx) + mod(IY+dy, nely)*nelx ...
                       + mod(IZ+dz, nelz)*nelx*nely + 1;

    % 8节点，顺序与 hexahedron 中节点坐标矩阵一致：
    % n1(-x,-y,-z)  n2(+x,-y,-z)  n3(+x,+y,-z)  n4(-x,+y,-z)
    % n5(-x,-y,+z)  n6(+x,-y,+z)  n7(+x,+y,+z)  n8(-x,+y,+z)
    n1=pnode(0,0,0); n2=pnode(1,0,0); n3=pnode(1,1,0); n4=pnode(0,1,0);
    n5=pnode(0,0,1); n6=pnode(1,0,1); n7=pnode(1,1,1); n8=pnode(0,1,1);

    % 展开为 nel×24 的 DOF 矩阵
    dof  = @(n) [3*n-2, 3*n-1, 3*n];
    edof = [dof(n1), dof(n2), dof(n3), dof(n4), ...
            dof(n5), dof(n6), dof(n7), dof(n8)];   % nel×24

    % ---- 材料属性场 ----
    lambda = zeros(nel, 1);
    mu     = zeros(nel, 1);
    vf = voxel(:);
    for k = 1:length(lambda_vec)
        mask      = (vf == k);
        lambda(mask) = lambda_vec(k);
        mu(mask)     = mu_vec(k);
    end

    % ---- 全局刚度矩阵 K (ndof×ndof) ----
    iK = kron(edof, ones(24,1))';
    jK = kron(edof, ones(1,24))';
    sK = keLambda(:) * lambda.' + keMu(:) * mu.';
    K  = sparse(iK(:), jK(:), sK(:), ndof, ndof);
    K  = 0.5 * (K + K');

    % ---- 全局载荷矩阵 F (ndof×6) ----
    % sF_raw: (24*6) × nel，重组为 (24*nel) × 6
    sF_raw = feLambda(:) * lambda.' + feMu(:) * mu.';  % (144) × nel
    sF_mat = reshape(sF_raw, 24, 6, nel);   % 24 × 6 × nel
    sF_mat = permute(sF_mat, [1, 3, 2]);    % 24 × nel × 6
    sF_mat = reshape(sF_mat, 24*nel, 6);    % (24*nel) × 6

    iF = repmat(edof(:), 1, 6);             % (24*nel) × 6
    jF = repmat(1:6, 24*nel, 1);            % (24*nel) × 6
    F  = sparse(iF(:), jF(:), sF_mat(:), ndof, 6);

    
    % ---- 求解特征位移 X (ndof×6) ----
    activedofs = sort(unique(edof(vf ~= 0, :)));
    X = zeros(ndof, 6);
    if length(activedofs) > 3
        fixeddofs = activedofs(1:3);
        freedofs  = setdiff(1:ndof, fixeddofs);
        
        K_free = K(freedofs, freedofs);
        F_free = F(freedofs, :);
        
        % ================== 正则化系数（新增可调功能） ==================
        if nargin < 7 || isempty(reg_coeff) || ~isnumeric(reg_coeff)
            reg_coeff = 1e-10;                 % 默认安全值
        end
        K_free = K_free + reg_coeff * speye(length(freedofs));
        % =================================================================
        
        % 求解（直接法，适合 30^3 网格）
        X(freedofs, :) = K_free \ F_free;
    end

    % ---- 基础应变场 X0 (nel×24×6) ----
    X0   = zeros(nel, 24, 6);
    X0_e = zeros(24, 6);
    ke   = keLambda + keMu;
    fe   = feLambda + feMu;
    idx  = [4, 7:11, 13:24];
    X0_e(idx, :) = ke(idx, idx) \ fe(idx, :);
    for c = 1:6
        X0(:,:,c) = repmat(X0_e(:,c)', nel, 1);
    end

    % ---- 均匀化刚度矩阵 CH (6×6) ----
    CH     = zeros(6);
    volume = lx * ly * lz;
    for i = 1:6
        for j = 1:6
            Xi = reshape(X(edof, i), nel, 24);
            Xj = reshape(X(edof, j), nel, 24);
            dXi = X0(:,:,i) - Xi;
            dXj = X0(:,:,j) - Xj;
            sum_L = sum(dXi * keLambda .* dXj, 2);
            sum_M = sum(dXi * keMu     .* dXj, 2);
            CH(i,j) = (1/volume) * sum(lambda .* sum_L + mu .* sum_M);
        end
    end
end


function [keLambda, keMu, feLambda, feMu] = hexahedron(a, b, c)
% 六面体单元矩阵（3×3×3 Gauss积分）
% 节点顺序：n1(-,-,-) n2(+,-,-) n3(+,+,-) n4(-,+,-)
%           n5(-,-,+) n6(+,-,+) n7(+,+,+) n8(-,+,+)

    CMu     = diag([2 2 2 1 1 1]);
    CLambda = zeros(6);  CLambda(1:3,1:3) = 1;

    xx = [-sqrt(3/5), 0, sqrt(3/5)];
    ww = [5/9, 8/9, 5/9];

    keLambda = zeros(24,24);  keMu = zeros(24,24);
    feLambda = zeros(24,6);   feMu = zeros(24,6);

    nodeCoords = [-a  a  a -a -a  a  a -a;
                  -b -b  b  b -b -b  b  b;
                  -c -c -c -c  c  c  c  c]';

    for ii = 1:3
        for jj = 1:3
            for kk = 1:3
                xi = xx(ii);  eta = xx(jj);  zeta = xx(kk);

                qx = [ -((eta-1)*(zeta-1))/8,  ((eta-1)*(zeta-1))/8, ...
                       -((eta+1)*(zeta-1))/8,  ((eta+1)*(zeta-1))/8, ...
                        ((eta-1)*(zeta+1))/8, -((eta-1)*(zeta+1))/8, ...
                        ((eta+1)*(zeta+1))/8, -((eta+1)*(zeta+1))/8 ];

                qy = [ -((xi-1)*(zeta-1))/8,  ((xi+1)*(zeta-1))/8, ...
                       -((xi+1)*(zeta-1))/8,  ((xi-1)*(zeta-1))/8, ...
                        ((xi-1)*(zeta+1))/8, -((xi+1)*(zeta+1))/8, ...
                        ((xi+1)*(zeta+1))/8, -((xi-1)*(zeta+1))/8 ];

                qz = [ -((xi-1)*(eta-1))/8,  ((xi+1)*(eta-1))/8, ...
                       -((xi+1)*(eta+1))/8,  ((xi-1)*(eta+1))/8, ...
                        ((xi-1)*(eta-1))/8, -((xi+1)*(eta-1))/8, ...
                        ((xi+1)*(eta+1))/8, -((xi-1)*(eta+1))/8 ];

                J    = [qx; qy; qz] * nodeCoords;
                qxyz = J \ [qx; qy; qz];

                B = zeros(6, 24);
                for ni = 1:8
                    col = (ni-1)*3 + (1:3);
                    B(:,col) = [ qxyz(1,ni)  0           0;
                                 0           qxyz(2,ni)  0;
                                 0           0           qxyz(3,ni);
                                 0           qxyz(3,ni)  qxyz(2,ni);
                                 qxyz(3,ni)  0           qxyz(1,ni);
                                 qxyz(2,ni)  qxyz(1,ni)  0 ];
                end

                w = det(J) * ww(ii) * ww(jj) * ww(kk);
                keLambda = keLambda + w * (B' * CLambda * B);
                keMu     = keMu     + w * (B' * CMu     * B);
                feLambda = feLambda + w * (B' * CLambda);
                feMu     = feMu     + w * (B' * CMu);
            end
        end
    end
end