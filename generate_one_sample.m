function generate_one_sample(s, N_total, materials, topology_types, ...
    grad_types, grad_directions, grid_size, L_range, t_range)
% generate_one_sample  生成单个梯度TPMS样本并保存为 .mat 文件
%
% param_field 第4维布局（共 2*n_topo+1 个通道）：
%   槽 1        ~ n_topo     : 混合权重 w（每种拓扑）
%   槽 n_topo+1 ~ 2*n_topo   : 水平集常数 t（每种拓扑）
%   槽 2*n_topo+1            : 单元胞尺寸 L

    % ---- 随机选材料对 ----
    mat_idx = randi(size(materials,1));
    mat = materials(mat_idx,:);
    E1 = mat{1}; nu1 = mat{2}; E2 = mat{3}; nu2 = mat{4};

    % ---- 随机选梯度类型和方向 ----
    grad_type = grad_types{randi(length(grad_types))};
    grad_dir  = grad_directions{randi(length(grad_directions))};

    % ---- 随机梯度参数 ----
    L_min = L_range(1) + rand * (L_range(2)-L_range(1)) * 0.3;
    L_max = L_min + rand * (L_range(2)-L_min);
    beta  = 2 + rand * 8;

    % ---- 拓扑权重采样（90%连续混合 + 10%单一拓扑）----
    n_topo = length(topology_types);
    if rand < 0.1
        topo_idx = randi(n_topo);
        w = zeros(1, n_topo);
        w(topo_idx) = 1;
    else
        w_temp = rand(1, n_topo);
        w = w_temp / sum(w_temp);
    end

    % 水平集常数（每种拓扑独立）
    t_vec = t_range(1) + (t_range(2)-t_range(1)) * rand(1, n_topo);

    % ---- 初始化场数组 ----
    % 修复：第4维正确大小为 2*n_topo+1（w×n + t×n + L×1）
    n_channels = 2*n_topo + 1;
    param_field = zeros([grid_size, n_channels]);
    perf_field  = zeros([grid_size, 9]);

    % 槽位索引（常量，不随体素变化）
    idx_w = 1:n_topo;
    idx_t = n_topo+1 : 2*n_topo;
    idx_L = 2*n_topo+1;

    % ---- 遍历空间网格点 ----
    for ix = 1:grid_size(1)
        for iy = 1:grid_size(2)
            for iz = 1:grid_size(3)
                x = (ix-1) / max(grid_size(1)-1, 1);
                y = (iy-1) / max(grid_size(2)-1, 1);
                z = (iz-1) / max(grid_size(3)-1, 1);

                % 局部单元胞尺寸
                L = compute_gradient_field(grad_type, grad_dir, x, y, z, L_min, L_max, beta);

                % 写入参数场（w和t在本样本内为全局常数，L随空间变化）
                param_field(ix, iy, iz, idx_w) = w;
                param_field(ix, iy, iz, idx_t) = t_vec;
                param_field(ix, iy, iz, idx_L) = L;

                % 计算局部等效刚度
                C = compute_stiffness_local(w, t_vec, L, [E1, nu1, E2, nu2], topology_types);
                C_gpa = C / 1e9;
                perf_field(ix, iy, iz, :) = [C_gpa(1,1), C_gpa(2,2), C_gpa(3,3), ...
                                              C_gpa(1,2), C_gpa(1,3), C_gpa(2,3), ...
                                              C_gpa(4,4), C_gpa(5,5), C_gpa(6,6)];
            end
        end
    end

    % ---- 保存文件 ----
    filename = sprintf('sample_%04d.mat', s);
    save(filename, 'perf_field', 'param_field', '-v7.3');

    if mod(s, 10) == 0 || s == N_total || s == 1
        fprintf('进度: %d/%d\n', s, N_total);
    end
end