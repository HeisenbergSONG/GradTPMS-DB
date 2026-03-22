function [dL_dx, dL_dy, dL_dz] = compute_gradient(grad_type, grad_dir, x, y, z, L_min, L_max, beta, grid_size)
% 计算局部梯度（用于高梯度区域标记）
    dx = 1/(grid_size(1)-1);
    dy = 1/(grid_size(2)-1);
    dz = 1/(grid_size(3)-1);
    
    dL_dx = (compute_L(grad_type, grad_dir, x+dx, y, z, L_min, L_max, beta) - ...
             compute_L(grad_type, grad_dir, x, y, z, L_min, L_max, beta)) / dx;
    dL_dy = (compute_L(grad_type, grad_dir, x, y+dy, z, L_min, L_max, beta) - ...
             compute_L(grad_type, grad_dir, x, y, z, L_min, L_max, beta)) / dy;
    dL_dz = (compute_L(grad_type, grad_dir, x, y, z+dz, L_min, L_max, beta) - ...
             compute_L(grad_type, grad_dir, x, y, z, L_min, L_max, beta)) / dz;
end