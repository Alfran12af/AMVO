function [u_an, v_an, px_an, py_an] = analytic_temporal( ...
    ux, uy, vx, vy, px, py, time, visc_nu, fu, fv, rho, N, h )
%
% Computes the analytical solution of velocity and pressure fields in time
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Inputs:
%   ux       : x-coordinates for horizontal velocity (u)
%   uy       : y-coordinates for horizontal velocity (u)
%   vx       : x-coordinates for vertical velocity (v)
%   vy       : y-coordinates for vertical velocity (v)
%   px       : x-coordinates for pressure field
%   py       : y-coordinates for pressure field
%   time     : current time instant
%   visc_nu  : kinematic viscosity
%   u_sym    : symbolic expression for horizontal velocity
%   v_sym    : symbolic expression for vertical velocity
%   rho      : fluid density
%   N        : number of cells
%   h        : grid spacing
%
% Outputs:
%   u_an  : analytical horizontal velocity field
%   v_an  : analytical vertical velocity field
%   px_an : analytical pressure gradient (x-direction)
%   py_an : analytical pressure gradient (y-direction)
%

    % Allocate memory
    u_an = zeros(N+2, N+2);
    v_an = zeros(N+2, N+2);
    p_field_an = zeros(N+2, N+2);

    % --- Temporal decay factor ---
    F = exp(-8*pi^2*visc_nu*time);

    % --- Compute analytical velocity and pressure fields ---
    for j = 1:N+2
        for i = 1:N+2

            % Velocity fields
            u_an(i,j) = F * fu(ux(i,j), uy(i,j));
            v_an(i,j) = F * fv(vx(i,j), vy(i,j));

            % Pressure field (analytical expression)
            p_field_an(i,j) = -rho * F^2 * ...
                (cos(2*pi*px(i,j))^2 + cos(2*pi*py(i,j))^2) / 2;

        end
    end

    % --- Compute analytical pressure gradient ---
    [px_an, py_an] = grad(p_field_an, h, N);

end