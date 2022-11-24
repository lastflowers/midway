function phi_x = skew(phi)
% validateattributes(phi,{'double'}, {'size',[3,1]});
phi_x = [0 , -phi(3), phi(2);
    phi(3), 0, -phi(1);
    -phi(2), phi(1), 0];
end