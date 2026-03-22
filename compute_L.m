function L = compute_L(grad_type, grad_dir, x, y, z, L_min, L_max, beta)
% 根据梯度类型和方向计算局部尺寸
    range = L_max - L_min;
    
    switch grad_type
        case 'linear'
            if contains(grad_dir, 'z')
                t = z;
            elseif contains(grad_dir, 'x')
                t = x;
            elseif contains(grad_dir, 'y')
                t = y;
            else
                t = (x + y + z)/3;
            end
            L = L_min + range * t;
            
        case 'quadratic'
            if contains(grad_dir, 'z')
                t = z^2;
            elseif contains(grad_dir, 'x')
                t = x^2;
            else
                t = (x^2 + y^2 + z^2)/3;
            end
            L = L_min + range * t;
            
        case 'sigmoid'
            if contains(grad_dir, 'z')
                t = (tanh(beta*(z-0.5)) + 1)/2;
            else
                t = (tanh(beta*((x+y+z)/3-0.5)) + 1)/2;
            end
            L = L_min + range * t;
            
        case 'radial'
            r = sqrt((x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2) * 2;
            L = L_min + range * min(1, r);
            
        case 'multi_axial'
            L = L_min + range * (x * y * z * 8);
    end
    L = max(L_min, min(L_max, L));
end