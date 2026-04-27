function [out_x, out_y] = grad(q, h, N)
% GRAD Computes the gradient of a scalar field on a staggered grid
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This function computes the numerical gradient of a scalar field using
%   a finite difference approximation on a structured grid.
%
%   The gradient is evaluated at staggered locations:
%     - ∂q/∂x is evaluated at u-velocity nodes
%     - ∂q/∂y is evaluated at v-velocity nodes
%
%   This is consistent with the projection method, where the pressure
%   gradient is used to correct the velocity field.
%
% Input:
%   q : scalar field (with halo, defined at cell centers)
%   h : grid spacing (uniform mesh)
%   N : number of control volumes per direction
%
% Output:
%   out_x : gradient in x-direction (located at u-nodes)
%   out_y : gradient in y-direction (located at v-nodes)
%
% Notes:
%   - Forward difference approximation is used
%   - Assumes periodic boundary conditions (halo nodes required)
%   - Output is used in velocity correction step (projection method)
%

    %% 1. Memory Allocation

    out_x = zeros(N+2, N+2);
    out_y = zeros(N+2, N+2);


    %% 2. Gradient Computation (Interior Nodes)

    for j = 2:N+1
        for i = 2:N+1

            % --- 2.1 Gradient in x-direction (u-nodes) ---
            out_x(i,j) = (q(i+1,j) - q(i,j)) / h;

            % --- 2.2 Gradient in y-direction (v-nodes) ---
            out_y(i,j) = (q(i,j+1) - q(i,j)) / h;

        end
    end


    %% 3. Halo Update (Periodic Boundary Conditions)

    out_x = halo_update(out_x);
    out_y = halo_update(out_y);

end