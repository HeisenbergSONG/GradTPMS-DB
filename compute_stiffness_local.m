function C = compute_stiffness_local(w, t, L, mat, topology_types)
% compute_stiffness_local  计算局部等效刚度矩阵 (Pa)
%
% 输入：
%   w              - 拓扑混合权重向量 [1×n_topo]，和为1
%   t              - 各拓扑水平集常数 [1×n_topo]
%   L              - 单元胞尺寸 (mm)
%   mat            - [E1(GPa), nu1, E2(GPa), nu2]
%   topology_types - 拓扑库 cell array，{名称, @(X,Y,Z)函数句柄}
%
% 输出：
%   C              - 6×6 等效刚度矩阵 (Pa)

    res = [30, 30, 30];

    % ---- 生成网格并计算混合TPMS场 ----
    xv = linspace(0, 2*pi, res(1));
    yv = linspace(0, 2*pi, res(2));
    zv = linspace(0, 2*pi, res(3));
    [XX, YY, ZZ] = meshgrid(xv, yv, zv);

    phi_mixed = compute_mixed_tpms(w, t, XX, YY, ZZ, topology_types);

    % ---- 体素化：phi<=0 为固相（材料1），其余为空洞 ----
    voxel = zeros(res);
    voxel(phi_mixed <= 0) = 1;

    % ---- 拉梅常数（Pa）----
    E1 = mat(1); nu1 = mat(2);
    lambda1 = E1*1e9 * nu1  / ((1+nu1) *(1-2*nu1));
    mu1     = E1*1e9        / (2*(1+nu1));

    E2 = mat(3); nu2 = mat(4);
    lambda2 = E2*1e9 * nu2  / ((1+nu2) *(1-2*nu2));
    mu2     = E2*1e9        / (2*(1+nu2));

    % ---- 均匀化 ----
    C_unit = homo3D_multi(1, 1, 1, [lambda1, lambda2], [mu1, mu2], voxel);

    % ---- 尺寸缩放（幂律 α=1.5）----
    C = C_unit * (L ^ (-1.5));
end


function phi_mixed = compute_mixed_tpms(w_list, t_list, X, Y, Z, topology_types)
% 计算归一化加权混合TPMS场
% 各拓扑先归一化到 [-1,1] 再加权，消除不同TPMS值域差异

    phi_mixed = zeros(size(X));
    for i = 1:length(w_list)
        if w_list(i) < 1e-10
            continue   % 跳过权重为0的拓扑，节省计算
        end
        % 取拓扑函数句柄
        func  = topology_types{i, 2};
        phi_i = func(X, Y, Z);

        % 归一化到 [-1, 1]
        pmin = min(phi_i(:));
        pmax = max(phi_i(:));
        if (pmax - pmin) > 1e-10
            phi_i = 2*(phi_i - pmin)/(pmax - pmin) - 1;
        end

        phi_mixed = phi_mixed + w_list(i) * (phi_i - t_list(i));
    end
end