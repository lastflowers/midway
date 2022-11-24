function xi = Log_se3(X)
% input : L \in SE(3)
% output : xi = [phi, rho] \in R^6
% Closed form reference at T.D. Barfoot, State Estimation for Robotics.
phi = dcm2rvec(X(1:3,1:3));
m = sqrt(sum(phi.^2, 1));
a = phi/m;
a_hat = [0, -a(3), a(2);
         a(3), 0, -a(1);
         -a(2), a(1), 0];
J = (sin(m)/m)*eye(3) + (1-sin(m)/m)*(a*a') + (1-cos(m))/m*a_hat;
if m ~= 0
    rho = J\X(1:3,4);
else 
    rho = X(1:3,4);
end
xi = [phi; rho];

end

