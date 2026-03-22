function phi = compute_tpms(topo_type, X, Y, Z, t, w)
% 计算指定拓扑类型的TPMS隐函数值
% topo_type: 拓扑名称字符串
% X,Y,Z: 网格坐标矩阵
% t: 水平集常数
% w: 该拓扑的混合权重（可选，用于后续混合）

    switch topo_type
        case 'Neovius'
            phi = 3*(cos(X)+cos(Y)+cos(Z)) + 4*cos(X).*cos(Y).*cos(Z);
            
        case 'SchoenIWP'
            phi = 2*(cos(X).*cos(Y) + cos(Y).*cos(Z) + cos(Z).*cos(X)) ...
                  - (cos(2*X) + cos(2*Y) + cos(2*Z));
            
        case 'SchwartzP'
            phi = cos(X) + cos(Y) + cos(Z);
            
        case 'Gyroid'
            phi = sin(X).*cos(Y) + sin(Y).*cos(Z) + sin(Z).*cos(X);
            
        case 'Diamond'
            phi = sin(X).*sin(Y).*sin(Z) + sin(X).*cos(Y).*cos(Z) ...
                + cos(X).*sin(Y).*cos(Z) + cos(X).*cos(Y).*sin(Z);
            
        case 'Lidinoid'
            phi = 0.5*(sin(2*X).*cos(Y).*sin(Z) + sin(2*Y).*cos(Z).*sin(X) + sin(2*Z).*cos(X).*sin(Y)) ...
                - 0.5*(cos(2*X).*cos(2*Y) + cos(2*Y).*cos(2*Z) + cos(2*Z).*cos(2*X)) + 0.15;
            
        case 'FischerKoch'
            phi = cos(2*X).*sin(Y).*cos(Z) + cos(2*Y).*sin(Z).*cos(X) + cos(2*Z).*sin(X).*cos(Y);
            
        otherwise
            error('未知拓扑类型: %s', topo_type);
    end
    
    % 应用水平集偏移
    phi = phi - t;
end

function phi_mixed = compute_mixed_tpms(w_list, t_list, X, Y, Z, topology_types)
% 计算混合TPMS的隐函数值
% w_list: 各拓扑权重向量（和为1）
% t_list: 各拓扑水平集常数向量
% topology_types: 拓扑名称元胞数组

    phi_mixed = zeros(size(X));
    for i = 1:length(w_list)
        phi_i = compute_tpms(topology_types{i}, X, Y, Z, t_list(i), w_list(i));
        phi_mixed = phi_mixed + w_list(i) * phi_i;
    end
end