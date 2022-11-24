function Ad_T = Adj_se3(T)
% This function calculates adjoint of T in SE(3)
% T X^ T^{-1} = (Ad_T X)^
% X = [rot; transl] in R^6 (Lie algebra)
R = T(1:3,1:3);
t = T(1:3,4);
Ad_T = zeros(6,6);
Ad_T(1:3,1:3) = R;
Ad_T(4:6,4:6) = R;
Ad_T(4:6,1:3) = skew(t)*R;
end

