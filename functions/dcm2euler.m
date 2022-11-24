function euler = dcm2euler(dcm)
% from direction cosine matrix to euler angle vector
% if you have dcm(C_{2}^{1}), this function calculates euler angle
% that rotates frame {1} by (z,psi), (y,theta), (x,phi) to frame {2}.
%%Input
% dcm  : C_{2}^{1}	[3x3] matrix
%%Output
% euler: [roll; pitch; heading] (rad) [3x1] vector, from {1} to {2}
% -pi/2 <= pitch <= pi/2
% if dcm(3,1) <= -0.999 then
%    (heading - roll) will be stored in heading and NaN in roll
% if dcm(3,1) >= 0.999
%    (heading + roll) will be stored in heading and NaN in roll
% Check Starpdown Analytics by Savage.

pitch = atan(-dcm(3,1)/sqrt(dcm(3,2)^2 + dcm(3,3)^2));

if dcm(3,1) <= -0.999
    roll = NaN;
    heading = atan2((dcm(2,3)-dcm(1,2)),(dcm(1,3)+dcm(2,2)));
elseif dcm(3,1) >= 0.999
    roll = NaN;
    heading = pi + atan2((dcm(2,3)+dcm(1,2)),(dcm(1,3)-dcm(2,2)));
else
    roll = atan2(dcm(3,2), dcm(3,3));
    heading = atan2(dcm(2,1), dcm(1,1));
end

euler = [roll; pitch; heading];
end

