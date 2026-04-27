function [dif_u, dif_v] = diffusive_2D(N, u, v, h)
%
% Computes the numerical diffusive terms of a velocity field
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Inputs:
%   N : number of cells
%   u : horizontal velocity field
%   v : vertical velocity field
%   h : grid spacing
%
% Outputs:
%   dif_u : numerical diffusive term for u
%   dif_v : numerical diffusive term for v

    % Allocate memory
    dif_u = zeros(N+2, N+2);
    dif_v = zeros(N+2, N+2);

    % Evaluate diffusive terms on the mesh
    for j = 2:1:N+1
        for i = 2:1:N+1

            % --- East derivatives ---
            dudx_e = u(i+1,j) - u(i,j);
            dvdx_e = v(i+1,j) - v(i,j);

            % --- West derivatives ---
            dudx_w = u(i,j) - u(i-1,j);
            dvdx_w = v(i,j) - v(i-1,j);

            % --- North derivatives ---
            dudy_n = u(i,j+1) - u(i,j);
            dvdy_n = v(i,j+1) - v(i,j);

            % --- South derivatives ---
            dudy_s = u(i,j) - u(i,j-1);
            dvdy_s = v(i,j) - v(i,j-1);

            % --- Diffusive terms (Laplacian discretization) ---
            dif_u(i,j) = dudx_e + dudy_n - dudx_w - dudy_s;
            dif_v(i,j) = dvdx_e + dvdy_n - dvdx_w - dvdy_s;

        end
    end

    % Normalize by grid spacing
    dif_u = dif_u./h^2;
    dif_v = dif_v./h^2;

end