function [nv, tv, pl] = nortanlen(N, x, z)

% -------------------------------------------------------------------------
% COMPUTE NORMAL AND TANGENT VECTORS + PANEL LENGTHS
%
% Computes the geometric properties of each panel:
%   - Panel length
%   - Unit normal vector
%   - Unit tangent vector
%
% Inputs:
%   N   -> number of panels
%   x   -> x-coordinates of nodes (N+1)
%   z   -> z-coordinates of nodes (N+1)
%
% Outputs:
%   nv  -> normal vectors at each panel (Nx2)
%   tv  -> tangent vectors at each panel (Nx2)
%   pl  -> panel lengths (Nx1)
% -------------------------------------------------------------------------

nv = zeros(N,2);
tv = zeros(N,2);
pl = zeros(N,1);

for i = 1:N
    
    % --- Panel length
    dx = x(i+1) - x(i);
    dz = z(i+1) - z(i);
    pl(i) = sqrt(dx^2 + dz^2);
    
    % --- Unit normal vector (pointing outward)
    nx =  (z(i) - z(i+1)) / pl(i);
    nz =  (x(i+1) - x(i)) / pl(i);
    nv(i,:) = [nx, nz];
    
    % --- Unit tangent vector (aligned with panel)
    tx =  nz;
    tz = -nx;
    tv(i,:) = [tx, tz];

end

end