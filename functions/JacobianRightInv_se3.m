function Jr_inv = JacobianRightInv_se3(xi)
    % This function computes the inverse of right Jacobian of SE(3)
    % upto 5th order
    ad_xi = [skew(xi(1:3)), zeros(3); 
             skew(xi(4:6)), skew(xi(1:3))];
    B = [1, -1/2, 1/6, 0, -1/30]; % Bernoulli numbers
    Jr_inv = zeros(6,6);
    for n = 0:length(B)-1
        if B(n+1) == 0
            continue;
        end
        Jr_inv = Jr_inv + B(n+1)/factorial(n)*(-ad_xi)^n;
    end
end

