function [dif_u, dif_v] = diffusive_2D(N, u, v, h)
% DIFFUSIVE_2D Computes the diffusive (viscous) terms of a 2D velocity field
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This function computes the numerical diffusive terms of the velocity
%   field using a finite control-volume (FCV) discretization on a staggered
%   grid.
%
%   The Laplacian operator is approximated using second-order central
%   differences, consistent with the formulation described in Part A.
%
% Input:
%   N : number of control volumes per direction
%   u : horizontal velocity field (staggered in x)
%   v : vertical velocity field (staggered in y)
%   h : grid spacing (uniform mesh)
%
% Output:
%   dif_u : diffusive term for u-velocity component
%   dif_v : diffusive term for v-velocity component
%
% Notes:
%   - Assumes periodic boundary conditions (halo nodes must be updated)
%   - Output corresponds to Laplacian(u) and Laplacian(v)
%

    %% 1. Memory Allocation

    dif_u = zeros(N+2, N+2);
    dif_v = zeros(N+2, N+2);


    %% 2. Diffusive Term Computation (Interior Nodes)

    for j = 2:N+1
        for i = 2:N+1

            % --- 2.1 East derivatives ---
            dudx_e = u(i+1,j) - u(i,j);
            dvdx_e = v(i+1,j) - v(i,j);

            % --- 2.2 West derivatives ---
            dudx_w = u(i,j) - u(i-1,j);
            dvdx_w = v(i,j) - v(i-1,j);

            % --- 2.3 North derivatives ---
            dudy_n = u(i,j+1) - u(i,j);
            dvdy_n = v(i,j+1) - v(i,j);

            % --- 2.4 South derivatives ---
            dudy_s = u(i,j) - u(i,j-1);
            dvdy_s = v(i,j) - v(i,j-1);

            % --- 2.5 Laplacian discretization ---
            dif_u(i,j) = dudx_e + dudy_n - dudx_w - dudy_s;
            dif_v(i,j) = dvdx_e + dvdy_n - dvdx_w - dvdy_s;

        end
    end


    %% 3. Normalization by Grid Spacing

    dif_u = dif_u ./ h^2;
    dif_v = dif_v ./ h^2;

end