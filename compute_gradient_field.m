function val = compute_gradient_field(grad_type, grad_dir, x, y, z, val_min, val_max, beta)
% 计算指定位置(x,y,z)的梯度场值
% 支持多种梯度类型和方向，值域[val_min, val_max]

    range = val_max - val_min;
    
    % 根据梯度方向计算基础位置参数t (0~1)
    switch grad_dir
        case 'z'
            t = z;
        case 'x'
            t = x;
        case 'y'
            t = y;
        case 'radial'
            r = sqrt((x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2) * 2;
            t = min(1, r);
        case 'coupled'
            t = (x + y + z) / 3;
        otherwise
            error('未知梯度方向: %s', grad_dir);
    end
    
    % 根据梯度类型计算实际值
    switch grad_type
        case 'linear'
            val = val_min + range * t;
            
        case 'quadratic'
            val = val_min + range * (t^2);
            
        case 'sigmoid'
            % t在0~1之间，sigmoid中心在0.5
            val = val_min + range * (tanh(beta*(t-0.5)) + 1)/2;
            
        case 'radial'
            % 径向梯度从中心向外增加
            r = sqrt((x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2) * 2;
            val = val_min + range * min(1, r^2);
            
        case 'multi_axial'
            % 多轴耦合：三个方向非线性组合
            t_coupled = x * y * z * 8;
            val = val_min + range * min(1, max(0, t_coupled));
            
        otherwise
            error('未知梯度类型: %s', grad_type);
    end
    
    % 确保在范围内
    val = max(val_min, min(val_max, val));
end