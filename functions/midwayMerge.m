function [w_ij, Tc_ij, To_ij, cov_ij] = midwayMerge(uw_i, Tc_i, To_i, cov_i, uw_j, Tc_j, To_j, cov_j)

w_i = uw_i / (uw_i + uw_j);
w_j = uw_j / (uw_i + uw_j);
w_ij = uw_i + uw_j;

diff_Tc = Log_se3(Tc_j/Tc_i);
dTc = w_j*diff_Tc;
Tc_ij = Exp_se3(dTc) * Tc_i;

% initialize memory
No = size(To_i,3);
diff_x = zeros(6*No,1);
To_ij = zeros(4,4,No);
Jl_inv_i = zeros(size(cov_i));
Jl_inv_j = zeros(size(cov_i));

diff_Tc_i = -w_j*diff_Tc;
diff_Tc_j = w_i*diff_Tc;
Jl_inv_i(1:6,1:6) = JacobianRightInv_se3(-diff_Tc_i);
Jl_inv_j(1:6,1:6) = JacobianRightInv_se3(-diff_Tc_j);
diff_x(1:6,1) = diff_Tc;

for n = 1:No
    diff_To = Log_se3(To_j(:,:,n)/To_i(:,:,n));
    dTo = w_j*diff_To;
    To_ij(:,:,n) = Exp_se3(dTo) * To_i(:,:,n);
    
    ridx = 6*n+1:6*n+6;
    
    diff_To_i = -w_j*diff_To;
    diff_To_j = w_i*diff_To;
    
    Jl_inv_i(ridx,ridx) = JacobianRightInv_se3(-diff_To_i);
    Jl_inv_j(ridx,ridx) = JacobianRightInv_se3(-diff_To_j);
    diff_x(ridx,1) = diff_To;
end

cov_ij = w_i*(Jl_inv_i*cov_i*Jl_inv_i') + w_j*(Jl_inv_j*cov_j*Jl_inv_j') + w_i*w_j*(diff_x*diff_x');

end