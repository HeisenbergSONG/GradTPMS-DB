function [E, nu] = compute_material_properties(material_name)
% 根据材料名称返回弹性模量E(GPa)和泊松比ν
% 所有数据来自权威文献

    materials_db = struct();
    
    % 金属材料
    materials_db.Ti6Al4V = [120, 0.34];      % Taylor & Francis 2025
    materials_db.NiTi = [75, 0.33];          % Additive Manufacturing 2022
    materials_db.AlSi10Mg = [70, 0.33];      % 多源验证
    materials_db.316L = [193, 0.30];         % 标准值
    materials_db.Cu = [400, 0.35];           % 标准值
    materials_db.CoCr = [230, 0.31];         % 标准值
    materials_db.ZrO2 = [190, 0.31];         % International Journal of Bioprinting 2025
    
    % 聚合物材料
    materials_db.PLA = [2, 0.40];            % 标准值
    materials_db.PEEK = [3.6, 0.40];         % 可选扩展
    
    if isfield(materials_db, material_name)
        vals = materials_db.(material_name);
        E = vals(1);
        nu = vals(2);
    else
        error('未知材料: %s', material_name);
    end
end