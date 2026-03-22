function [latest_dir, manifest_info] = load_latest_m()
    % 加载最新生成的TPMS数据集目录

    dirs = dir('dataset_*');
    if isempty(dirs)
        error('没有找到任何 dataset_ 开头的目录');
    end

    [~, idx] = max([dirs.datenum]);
    latest_dir = fullfile(pwd, dirs(idx).name);

    % 读取manifest最后一段
    if exist('manifest.txt', 'file')
        manifest = fileread('manifest.txt');
        parts = strsplit(manifest, '----------------------------------------');
        latest_block = strtrim(parts{end-1});
        manifest_info = latest_block;
    else
        manifest_info = '无 manifest 记录';
    end

    fprintf('最新数据目录: %s\n', latest_dir);
    disp(manifest_info);
end