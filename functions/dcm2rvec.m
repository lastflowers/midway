function rvec = dcm2rvec(R)
% from direction cosine matrix to rotation vector
% if you have dcm(C_{2}^{1}), this function calculates rotation vector({1}->{2}).
%%Input
% R: C_{2}^{1}	[3x3] matrix
%%Output
% rvec: rho_{2}^{1} [3x1] vector, from {1} to {2}

tr  = trace(R);
if abs(tr - 3) == 0
    rvec = zeros(3,1);
else
    tmp = (tr-1)/2;
    theta = acos(tmp);
    k = 1/2/sin(theta)*[R(3,2) - R(2,3); R(1,3) - R(3,1); R(2,1) - R(1,2)];
    rvec = theta *k ;
end
end