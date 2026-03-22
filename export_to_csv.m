% export_to_csv.m
% 将所有 dataset_*/sample_*.mat 文件汇总为一个结构化 CSV
% 适用于你的梯度TPMS多材料项目（2026年3月版本）

clear; clc;

%% 1. 找到最新/所有 dataset_ 目录
dirs = dir('dataset_*');
if isempty(dirs)
    error('没有找到任何以 "dataset_" 开头的目录');
end

% 按时间降序排序（最新在前）
[~, idx] = sort([dirs.datenum], 'descend');
dirs = dirs(idx);

% 输出文件名（可按需修改）
output_csv = 'tpms_gradient_dataset_2026.csv';

%% 2. 定义拓扑名称顺序（必须与 generate_one_sample.m 中一致！）
topology_names = {'Neovius', 'SchoenIWP', 'SchwartzP', 'Gyroid', ...
                  'Diamond', 'Lidinoid', 'FischerKoch'};
n_topo = length(topology_names);

%% 3. 材料映射表（示例，根据你的 materials cell 数组完善）
% 键：E1,E2 的组合字符串；值：可读名称
material_map = containers.Map(...
    { '120_2',   '75_2',   '70_2',   '193_2',   '120_70',  '190_2',   '400_70',  '230_2' }, ...
    { 'Ti6Al4V-PLA', 'NiTi-PLA', 'AlSi10Mg-PLA', '316L-PLA', ...
      'Ti6Al4V-AlSi10Mg', 'ZrO2-PLA', 'Cu-AlSi10Mg', 'CoCr-PLA' });

%% 4. 准备表结构
headers = [ ...
    {'sample_id', 'material_pair', 'E1_GPa', 'nu1', 'E2_GPa', 'nu2'}, ...
    strcat('w_', topology_names), ...
    strcat('t_', topology_names), ...
    {'L_mean', 'L_min', 'L_max', 'L_gradient_type', 'L_gradient_dir', 'beta'}, ...
    {'C11_mean','C22_mean','C33_mean','C12_mean','C13_mean','C23_mean', ...
     'C44_mean','C55_mean','C66_mean'}, ...
    {'volume_fraction_est', 'timestamp'} ];

data_table = table('Size', [0, length(headers)], ...
                   'VariableTypes', repmat({'double'},1,length(headers)), ...
                   'VariableNames', headers);

% 把前几列改为 cell 以存字符串
data_table.sample_id = cell(0,1);
data_table.material_pair = cell(0,1);
data_table.L_gradient_type = cell(0,1);
data_table.L_gradient_dir  = cell(0,1);
data_table.timestamp = cell(0,1);

fprintf('开始扫描并转换数据集...\n');

%% 5. 遍历所有 dataset 目录和 mat 文件
total_files = 0;
for d = 1:length(dirs)
    data_dir = fullfile(dirs(d).folder, dirs(d).name);
    mat_files = dir(fullfile(data_dir, 'sample_*.mat'));
    
    fprintf('处理目录: %s  (找到 %d 个样本)\n', dirs(d).name, length(mat_files));
    
    for f = 1:length(mat_files)
        total_files = total_files + 1;
        filename = fullfile(data_dir, mat_files(f).name);
        % fprintf('  %3d | %s\n', total_files, mat_files(f).name);
        
        try
            loaded = load(filename);
            param = loaded.param_field;   % [nx ny nz  2*n_topo+1]
            perf  = loaded.perf_field;    % [nx ny nz  9]
            
            % 提取全局参数（所有体素相同）
            w_vec = squeeze(param(1,1,1,1:n_topo));
            t_vec = squeeze(param(1,1,1,n_topo+1:2*n_topo));
            L_field = squeeze(param(:,:,:,end));
            
            % 材料参数（这里假设你保存了或能推断）
            % 实际项目中建议在 generate_one_sample 中 save 时加入：
            % save(..., 'material_idx', 'grad_type', 'grad_dir', 'beta', ...)
            E1 = 120; nu1 = 0.34; E2 = 2.0; nu2 = 0.40;   % 默认示例
            material_key = sprintf('%g_%g', E1, E2);
            if isKey(material_map, material_key)
                mat_pair = material_map(material_key);
            else
                mat_pair = sprintf('E1=%.1f_E2=%.1f', E1, E2);
            end
            
            % L 统计
            L_mean = mean(L_field(:));
            L_min  = min(L_field(:));
            L_max  = max(L_field(:));
            
            % 刚度平均（GPa）
            C_mat = mean(reshape(perf, [], 9), 1);
            
            % 体积分数粗估（可替换为真实 phi<=0 比例）
            vf_est = 0.25 + 0.5 * exp(-mean(t_vec)) * sum(w_vec.^2);  % 示例
            
            % 梯度信息（当前 .mat 缺少，需补充保存）
            grad_type = 'unknown';   % ← 待完善
            grad_dir  = 'unknown';
            beta      = NaN;
            
            % 构建一行数据
            row_data = { ...
                mat_files(f).name(1:end-4), ...   % sample_id
                mat_pair, E1, nu1, E2, nu2, ...
                w_vec', t_vec', ...
                L_mean, L_min, L_max, grad_type, grad_dir, beta, ...
                C_mat(1), C_mat(2), C_mat(3), C_mat(4), C_mat(5), C_mat(6), ...
                C_mat(7), C_mat(8), C_mat(9), ...
                vf_est, datestr(now, 'yyyy-mm-dd HH:MM:SS') };
            
            % 追加到 table
            new_row = cell2table(row_data, 'VariableNames', headers);
            data_table = [data_table; new_row];
            
        catch ME
            fprintf('  跳过出错文件 %s → %s\n', mat_files(f).name, ME.message);
        end
    end
end

%% 6. 保存 CSV
writetable(data_table, output_csv);
fprintf('\n转换完成！\n');
fprintf('总样本数：%d\n', height(data_table));
fprintf('输出文件：%s\n', output_csv);
fprintf('列数：%d\n', width(data_table));