function div = diverg(u, v, h, N)
% DIVERG Computes the divergence of a 2D velocity field (FCV formulation)
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This function computes the numerical divergence of a velocity field
%   using a finite control-volume (FCV) formulation on a staggered grid.
%
%   The divergence is evaluated in integral form:
%
%       ∫_V ∇·u dV ≈ (u_e - u_w + v_n - v_s) · Δ
%
%   where fluxes are evaluated at the faces of each control volume.
%
%   The result corresponds to the integral of the divergence over each
%   control volume (not pointwise divergence).
%
% Input:
%   u : horizontal velocity field (staggered in x)
%   v : vertical velocity field (staggered in y)
%   h : grid spacing (uniform mesh)
%   N : number of control volumes per direction
%
% Output:
%   div : divergence field (integral form, centered grid)
%
% Notes:
%   - Assumes periodic boundary conditions (halo nodes required)
%   - Output is used as RHS of Poisson equation in projection method
%

    %% 1. Memory Allocation

    div = zeros(N+2, N+2);


    %% 2. Divergence Computation (Interior Nodes)

    for j = 2:N+1
        for i = 2:N+1

            %% 2.1 Face Velocities

            u_p = u(i,j);       % East face
            u_w = u(i-1,j);     % West face

            v_p = v(i,j);       % North face
            v_s = v(i,j-1);     % South face


            %% 2.2 Finite Volume Divergence

            div(i,j) = h * (u_p - u_w + v_p - v_s);

        end
    end


    %% 3. Halo Update (Periodic Boundary Conditions)

    div = halo_update(div);

end