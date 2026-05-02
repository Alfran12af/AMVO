function [r1,r2,theta1,theta2] = rtheta(plj,xcpan,zcpan)
% Outputs
%   r1: distance between nodal and control point
%   r2:distance between control point and next nodal point
%   theta1: angle between r1 and panel
%   theta2: angle between r2 and panel

r1=sqrt(xcpan^2+zcpan^2);

r2=sqrt((xcpan-plj)^2+zcpan^2);

theta1 = atan2(zcpan,xcpan);

theta2 = atan2(zcpan,xcpan-plj);
end
