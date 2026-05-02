function [u_field, v_field] = set_velocity_field(N, u, v, ux, uy, vx, vy)
%
% Evaluates the velocity field from symbolic expressions on the mesh
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Inputs:
%   N  : number of cells
%   u  : symbolic expression of horizontal velocity
%   v  : symbolic expression of vertical velocity
%   ux : x-coordinates of u velocity nodes
%   uy : y-coordinates of u velocity nodes
%   vx : x-coordinates of v velocity nodes
%   vy : y-coordinates of v velocity nodes
%
% Outputs:
%   u_field : evaluated horizontal velocity field
%   v_field : evaluated vertical velocity field

    % Symbolic variables
    syms x y

    % --- Convert symbolic expressions to MATLAB functions ---
    u_fun = matlabFunction(u, 'Vars', [x y]);
    v_fun = matlabFunction(v, 'Vars', [x y]);

    % --- Evaluate velocity fields on the mesh ---
    u_field = u_fun(ux, uy);
    v_field = v_fun(vx, vy);

    % --- Update halo (ghost cells) ---
    u_field = halo_update(u_field);
    v_field = halo_update(v_field);

end