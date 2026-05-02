function dt = timestep(N, L, u_field, v_field, visc)
%
% Computes the stable time step based on CFL condition
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Inputs:
%   N        : number of cells
%   L        : domain length
%   u_field  : horizontal velocity field
%   v_field  : vertical velocity field
%   visc     : kinematic viscosity
%
% Outputs:
%   dt       : time step satisfying stability conditions
%

    % Grid spacing
    h = L/N;

    % --- Maximum velocities ---
    u_max = max(abs(u_field(:)));
    v_max = max(abs(v_field(:)));

    % --- Convective time limit ---
    dtconv_x = h / u_max;
    dtconv_y = h / v_max;
    dtconv = min(dtconv_x, dtconv_y);

    % --- Diffusive time limit ---
    dtdiff = 0.5 * h^2 / visc;

    % --- Final time step (CFL condition) ---
    f = 0.1; % safety factor
    dt = f * min(dtconv, dtdiff);

end