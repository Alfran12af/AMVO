function [cu, cv] = convective_2D(u, v, h, N)
% CONVECTIVE_2D Computes the convective terms of a 2D velocity field (FCV)
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This function computes the numerical convective terms of a velocity
%   field using a finite control-volume (FCV) discretization on a staggered
%   grid.
%
%   The convective term is evaluated in conservative (divergence) form:
%
%       ∇·(u_i u)
%
%   Fluxes are computed at control-volume faces using midpoint
%   interpolation (second-order accurate), ensuring consistency with the
%   analytical formulation used in the Method of Manufactured Solutions (MMS).
%
% Input:
%   u : horizontal velocity field (staggered in x)
%   v : vertical velocity field (staggered in y)
%   h : grid spacing (uniform mesh)
%   N : number of control volumes per direction
%
% Output:
%   cu : convective term for u-velocity component
%   cv : convective term for v-velocity component
%
% Notes:
%   - Conservative formulation (flux-based) is used
%   - Face values obtained via central interpolation
%   - Periodic boundary conditions assumed (halo nodes required)
%

    %% 1. Memory Allocation

    cu = zeros(N+2, N+2);
    cv = zeros(N+2, N+2);


    %% 2. Convective Term Computation (Interior Nodes)

    for i = 2:N+1
        for j = 2:N+1

            %% 2.1 Face Velocities (Midpoint Interpolation)

            ue = (u(i+1,j) + u(i,j)) / 2;
            uw = (u(i-1,j) + u(i,j)) / 2;
            un = (u(i,j+1) + u(i,j)) / 2;
            us = (u(i,j-1) + u(i,j)) / 2;

            ve = (v(i+1,j) + v(i,j)) / 2;
            vw = (v(i-1,j) + v(i,j)) / 2;
            vn = (v(i,j+1) + v(i,j)) / 2;
            vs = (v(i,j-1) + v(i,j)) / 2;


            %% 2.2 Mass Fluxes at Faces

            % Horizontal momentum control volume
            Feh = h * (u(i+1,j) + u(i,j)) / 2;
            Fwh = h * (u(i-1,j) + u(i,j)) / 2;
            Fnh = h * (v(i,j) + v(i+1,j)) / 2;
            Fsh = h * (v(i,j-1) + v(i+1,j-1)) / 2;

            % Vertical momentum control volume
            Fev = h * (u(i,j+1) + u(i,j)) / 2;
            Fwv = h * (u(i-1,j+1) + u(i-1,j)) / 2;
            Fnv = h * (v(i,j+1) + v(i,j)) / 2;
            Fsv = h * (v(i,j) + v(i,j-1)) / 2;


            %% 2.3 Convective Flux Balance

            cu(i,j) = ue * Feh - uw * Fwh + un * Fnh - us * Fsh;
            cv(i,j) = ve * Fev - vw * Fwv + vn * Fnv - vs * Fsv;

        end
    end


    %% 3. Normalization by Control Volume Area

    cu = cu ./ h^2;
    cv = cv ./ h^2;

end