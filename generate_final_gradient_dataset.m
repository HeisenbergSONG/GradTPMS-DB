% generate_final_gradient_dataset.m
% 完整版：三维梯度多材料混合TPMS数据库生成（带监控、邮件、错误捕获）
% 设计维度：8种材料对 × 5种梯度函数 × 5种梯度方向 × 连续拓扑混合 × 连续尺寸/壁厚

clear; clc;

% ========== 配置 ==========
email_to   = '你的监控邮箱@qq.com';           % ← 改成你自己的
email_from = '你的gmail@gmail.com';            % ← 改成你自己的

run_timestamp = datestr(now, 'yyyymmdd_HHMMSS');
data_dir = fullfile(pwd, ['dataset_', run_timestamp]);
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

manifest_file = fullfile(pwd, 'manifest.txt');

% 确保当前目录在路径中（worker 继承 addpath，但不继承 cd）
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
addpath(script_dir);

% ========== 1.并行池（28核） ==========

pool = gcp('nocreate');   % 获取当前池（如果有）

if isempty(pool)
    fprintf('没有找到现有并行池，正在启动 %d 个 worker...\n', 28);
    try
        pool = parpool('Processes', 28);   % 显式使用 'Processes' profile（更稳定）
        % 或者 pool = parpool('local', 28);  % 如果你有自定义 local profile
        fprintf('并行池启动成功，已连接 %d 个 worker\n', pool.NumWorkers);
    catch ME
        fprintf('并行池启动失败：%s\n', ME.message);
        fprintf('将使用单线程模式继续（性能会大幅下降）\n');
        pool = [];   % 明确置空，避免后续误用
    end
else
    fprintf('复用已有并行池，已有 %d 个 worker\n', pool.NumWorkers);
end

% ---- 只有当 pool 非空时才附加文件 ----
if ~isempty(pool) && isa(pool, 'parallel.Pool')
    % 你的文件列表收集逻辑（已有代码很好）
    needed_files = { ...
        fullfile(script_dir, 'homo3D_multi.m'), ...
        fullfile(script_dir, 'compute_stiffness_local.m'), ...
        fullfile(script_dir, 'compute_gradient_field.m'), ...
        fullfile(script_dir, 'compute_L.m'), ...
        fullfile(script_dir, 'compute_tpms.m'), ...
        fullfile(script_dir, 'generate_one_sample.m') };
    
    existing = {};
    for k = 1:length(needed_files)
        if exist(needed_files{k}, 'file')
            existing{end+1} = needed_files{k};
        else
            fprintf('警告：文件缺失，无法附加 → %s\n', needed_files{k});
        end
    end
    
    if ~isempty(existing)
        addAttachedFiles(pool, existing);
        fprintf('成功附加 %d 个依赖文件到并行池\n', length(existing));
    else
        fprintf('没有有效文件需要附加\n');
    end
else
    fprintf('无有效并行池，跳过文件附加（parfor 可能因找不到函数而失败）\n');
end

% ========== 2. 材料库+(按实际的来) ==========
materials = {
    120, 0.34,   2, 0.40, "Ti6Al4V",    "PLA",      "biomedical_aerospace";
    75,  0.33,   2, 0.40, "NiTi",        "PLA",      "smart_structure";
    70,  0.33,   2, 0.40, "AlSi10Mg",    "PLA",      "lightweight_thermal";
    193, 0.30,   2, 0.40, "316L",        "PLA",      "high_strength";
    120, 0.34,  70, 0.33, "Ti6Al4V",    "AlSi10Mg", "aerospace_bimetal";
    190, 0.31,   2, 0.40, "ZrO2",        "PLA",      "dental_implant";
    400, 0.35,  70, 0.33, "Cu",          "AlSi10Mg", "thermal_management";
    230, 0.31,   2, 0.40, "CoCr",        "PLA",      "orthopedic_implant";
};

% ========== 3. 参数设置(按需调整) ==========
N_samples      = 200;
grid_size      = [48, 48, 48];
L_range        = [0.5, 10.0];
t_range        = [-1.5, 1.5];
grad_types     = {'linear', 'quadratic', 'sigmoid', 'radial', 'multi_axial'};
grad_directions = {'z', 'x', 'y', 'radial', 'coupled'};

% ========== 4. 拓扑类型定义 (按实际的来)==========
topology_types = {
    'Neovius',     @(X,Y,Z) 3*(cos(X)+cos(Y)+cos(Z)) + 4*cos(X).*cos(Y).*cos(Z);
    'SchoenIWP',   @(X,Y,Z) 2*(cos(X).*cos(Y)+cos(Y).*cos(Z)+cos(Z).*cos(X)) - (cos(2*X)+cos(2*Y)+cos(2*Z));
    'SchwartzP',   @(X,Y,Z) cos(X) + cos(Y) + cos(Z);
    'Gyroid',      @(X,Y,Z) sin(X).*cos(Y) + sin(Y).*cos(Z) + sin(Z).*cos(X);
    'Diamond',     @(X,Y,Z) sin(X).*sin(Y).*sin(Z) + sin(X).*cos(Y).*cos(Z) + cos(X).*sin(Y).*cos(Z) + cos(X).*cos(Y).*sin(Z);
    'Lidinoid',    @(X,Y,Z) 0.5*(sin(2*X).*cos(Y).*sin(Z) + sin(2*Y).*cos(Z).*sin(X) + sin(2*Z).*cos(X).*sin(Y)) - 0.5*(cos(2*X).*cos(2*Y) + cos(2*Y).*cos(2*Z) + cos(2*Z).*cos(2*X)) + 0.15;
    'FischerKoch', @(X,Y,Z) cos(2*X).*sin(Y).*cos(Z) + cos(2*Y).*sin(Z).*cos(X) + cos(2*Z).*sin(X).*cos(Y);
};

% ========== 输出配置信息 ==========
fprintf('开始生成 - 时间戳: %s\n', run_timestamp);
fprintf('输出目录: %s\n', data_dir);
fprintf('拓扑库包含 %d 种基础类型\n', length(topology_types));
fprintf('\n========================================\n');
fprintf('三维梯度多材料混合TPMS数据库生成\n');
fprintf('========================================\n');
fprintf('材料对数: %d\n', size(materials,1));
fprintf('拓扑类型: %d\n', length(topology_types));
fprintf('梯度函数: %d种 (%s)\n', length(grad_types), strjoin(grad_types, ', '));
fprintf('梯度方向: %d种 (%s)\n', length(grad_directions), strjoin(grad_directions, ', '));
fprintf('网格分辨率: %d x %d x %d\n', grid_size(1), grid_size(2), grid_size(3));
fprintf('单元胞尺寸范围: %.1f - %.1f mm\n', L_range(1), L_range(2));
fprintf('目标样本数: %d\n', N_samples);
fprintf('========================================\n');
fprintf('开始生成数据...\n\n');

tic;

% ========== 主循环（加错误捕获） ==========
try
    parfor s = 1:N_samples
        generate_one_sample(s, N_samples, materials, topology_types, ...
            grad_types, grad_directions, grid_size, L_range, t_range, ...
            data_dir);   % ← 新增 data_dir 参数
    end

    elapsed_time = toc;

    % 正常完成记录 + 邮件
    fid = fopen(manifest_file, 'a');
    fprintf(fid, 'Run started: %s\n  Directory: %s\n  Samples: %d\n  Elapsed: %.1f min\n----------------------------------------\n', ...
            run_timestamp, data_dir, N_samples, elapsed_time/60);
    fclose(fid);

    sendmail(email_to, 'TPMS数据集生成完成', ...
        sprintf('Run %s 完成\n样本数: %d\n耗时: %.1f 分钟\n目录: %s', ...
                run_timestamp, N_samples, elapsed_time/60, data_dir));

catch ME
    err_msg = getReport(ME, 'extended', 'hyperlinks', 'off');
    error_save_file = fullfile(pwd, sprintf('ERROR_state_%s.mat', datestr(now,'yyyymmdd_HHMMSS')));
    save(error_save_file, '-v7.3');

    body = sprintf('【严重】TPMS生成出错！\n时间: %s\n错误:\n%s\n\n已保存状态: %s', ...
                   datestr(now), err_msg, error_save_file);
    sendmail(email_to, 'TPMS生成异常告警', body);

    rethrow(ME);
end

fprintf('生成完成！目录: %s\n', data_dir);