function F = halo_update(F)
% HALO_UPDATE Updates halo nodes using periodic boundary conditions
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This function updates the halo (ghost cells) of a 2D field assuming
%   periodic boundary conditions in both spatial directions.
%
%   The values at the boundaries are copied from the opposite side of the
%   domain, ensuring continuity across periodic interfaces.
%
%   This operation is required to correctly evaluate derivatives and fluxes
%   near the domain boundaries in a finite control-volume (FCV) framework.
%
% Input:
%   F : scalar or vector field including halo nodes (size: (N+2) x (N+2))
%
% Output:
%   F : field with updated halo values
%
% Notes:
%   - Assumes periodic domain in both x and y directions
%   - Interior nodes are indexed from 2 to N+1
%   - Halo nodes are located at indices 1 and N+2
%

    %% 1. Mesh Size

    N = size(F,1) - 2;


    %% 2. Horizontal Halo Update (x-direction)

    % Left boundary ← right interior
    F(1,:)   = F(N+1,:);

    % Right boundary ← left interior
    F(N+2,:) = F(2,:);


    %% 3. Vertical Halo Update (y-direction)

    % Bottom boundary ← top interior
    F(:,1)   = F(:,N+1);

    % Top boundary ← bottom interior
    F(:,N+2) = F(:,2);

end