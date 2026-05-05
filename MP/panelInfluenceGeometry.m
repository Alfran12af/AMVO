function [r1, r2, theta1, theta2] = panelInfluenceGeometry(plj, xcpan, zcpan)

% -------------------------------------------------------------------------
% PANEL INFLUENCE GEOMETRY
%
% Computes geometric terms required for panel influence calculations:
%   - Distances from control point to panel endpoints
%   - Angular terms for vortex/source integrals
%
% Inputs:
%   plj    -> length of panel j
%   xcpan  -> local x-coordinate of control point (panel reference frame)
%   zcpan  -> local z-coordinate of control point (panel reference frame)
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

% --- Angles (robust quadrant handling)
theta1 = atan2(zcpan, xcpan);
theta2 = atan2(zcpan, xcpan - plj);

end