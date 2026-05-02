function [r1, r2, theta1, theta2] = rtheta(plj, xcpan, zcpan)

% -------------------------------------------------------------------------
% COMPUTE GEOMETRIC TERMS FOR PANEL INFLUENCE
%
% Computes distances and angular terms between a control point and a panel
% expressed in the local coordinate system of the panel.
%
% Inputs:
%   plj    -> length of panel j
%   xcpan  -> local x-coordinate of control point relative to panel j
%   zcpan  -> local z-coordinate of control point relative to panel j
%
% Outputs:
%   r1     -> distance from control point to panel start node
%   r2     -> distance from control point to panel end node
%   theta1 -> angle between r1 and panel axis
%   theta2 -> angle between r2 and panel axis
% -------------------------------------------------------------------------

% --- Distances
r1 = sqrt(xcpan^2 + zcpan^2);
r2 = sqrt((xcpan - plj)^2 + zcpan^2);

% --- Angles (using atan2 for correct quadrant)
theta1 = atan2(zcpan, xcpan);
theta2 = atan2(zcpan, xcpan - plj);

end