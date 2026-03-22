function C = compute_equivalent_stiffness(voxel, E1, nu1, E2, nu2)
% 直接调用homo3D计算等效刚度
% 多材料通过体积分数平均处理

    V1 = sum(voxel(:) == 1) / numel(voxel);
    V2 = 1 - V1;
    
    % Voigt-Reuss-Hill平均
    E_avg = (V1*E1 + V2*E2 + 1/(V1/E1 + V2/E2)) / 2;
    nu_avg = V1*nu1 + V2*nu2;
    
    lambda = E_avg*1e9 * nu_avg / ((1+nu_avg)*(1-2*nu_avg));
    mu = E_avg*1e9 / (2*(1+nu_avg));
    
    C = homo3D(1, 1, 1, lambda, mu, voxel);
end