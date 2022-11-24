function L = Exp_se3( xi )
% input : xi = [phi, rho] \in R^6
% output : L \in SE(3)
% Closed form reference at T.D. Barfoot, State Estimation for Robotics.
if sum(xi) ~= 0
    phi = sqrt(sum(xi(1:3).^2, 1));
    xi_hat = zeros(4,4);
    xi_hat(1:3,1:3) = skew(xi(1:3));
    xi_hat(1:3,4) = xi(4:6);
    if phi == 0
        L = eye(4) + xi_hat;
    else
        L = eye(4) + xi_hat + ((1-cos(phi))/phi^2)*xi_hat*xi_hat ...
            + ((phi-sin(phi))/phi^3)*xi_hat*xi_hat*xi_hat;
    end
else
    L = eye(4);
end
end

