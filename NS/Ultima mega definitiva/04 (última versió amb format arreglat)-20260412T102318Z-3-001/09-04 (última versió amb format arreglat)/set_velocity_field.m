function [u_field, v_field] = set_velocity_field(N, u, v, ux, uy, vx, vy)
% SET_VELOCITY_FIELD Evaluates velocity field from symbolic expressions
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This function evaluates the analytical (manufactured) velocity field
%   on a staggered grid by converting symbolic expressions into numerical
%   functions.
%
%   The velocity components are evaluated at their respective locations:
%     - u is evaluated at u-nodes (vertical faces)
%     - v is evaluated at v-nodes (horizontal faces)
%
%   This is used in the Method of Manufactured Solutions (MMS) to generate
%   consistent velocity fields for verification.
%
% Input:
%   N  : number of control volumes per direction
%   u  : symbolic expression of horizontal velocity
%   v  : symbolic expression of vertical velocity
%   ux : x-coordinates of u-velocity nodes
%   uy : y-coordinates of u-velocity nodes
%   vx : x-coordinates of v-velocity nodes
%   vy : y-coordinates of v-velocity nodes
%
% Output:
%   u_field : evaluated horizontal velocity field (staggered in x)
%   v_field : evaluated vertical velocity field (staggered in y)
%
% Notes:
%   - Uses MATLAB symbolic toolbox (matlabFunction)
%   - Assumes periodic boundary conditions (halo update required)
%   - Fields include halo nodes (size: (N+2) x (N+2))
%

    %% 1. Symbolic Variables

    syms x y


    %% 2. Symbolic → Numerical Functions

    u_fun = matlabFunction(u, 'Vars', [x y]);
    v_fun = matlabFunction(v, 'Vars', [x y]);


    %% 3. Memory Allocation

    u_field = zeros(N+2, N+2);
    v_field = zeros(N+2, N+2);


    %% 4. Field Evaluation (Interior Nodes)

    for i = 2:N+1
        for j = 2:N+1

            % --- 4.1 u-component (u-nodes) ---
            u_field(i,j) = u_fun(ux(i,j), uy(i,j));

            % --- 4.2 v-component (v-nodes) ---
            v_field(i,j) = v_fun(vx(i,j), vy(i,j));

        end
    end


    %% 5. Halo Update (Periodic Boundary Conditions)

    u_field = halo_update(u_field);
    v_field = halo_update(v_field);

end