function [u_an, v_an, px_an, py_an] = analytic_temporal( ...
    ux, uy, vx, vy, px, py, time, visc_nu, fu, fv, rho, N, h )
% ANALYTIC_TEMPORAL Computes the analytical solution of velocity and pressure fields in time
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This function evaluates the analytical (manufactured) solution of the
%   2D incompressible Navier–Stokes equations at a given time.
%
%   The velocity field follows an exponential temporal decay, while the
%   pressure gradient is derived analytically from the exact solution.
%
%   The solution is evaluated on a staggered grid, consistent with the
%   numerical discretization used in Part C.
%
% Input:
%   ux       : x-coordinates of u-velocity nodes (staggered grid)
%   uy       : y-coordinates of u-velocity nodes
%   vx       : x-coordinates of v-velocity nodes
%   vy       : y-coordinates of v-velocity nodes
%   px       : x-coordinates of pressure nodes
%   py       : y-coordinates of pressure nodes
%   time     : current time instant
%   visc_nu  : kinematic viscosity
%   fu       : function handle for horizontal velocity (symbolic → numeric)
%   fv       : function handle for vertical velocity (symbolic → numeric)
%   rho      : fluid density
%   N        : number of control volumes per direction
%   h        : grid spacing (uniform mesh)
%
% Output:
%   u_an  : analytical horizontal velocity field
%   v_an  : analytical vertical velocity field
%   px_an : analytical pressure gradient in x-direction
%   py_an : analytical pressure gradient in y-direction
%
% Notes:
%   - Assumes periodic boundary conditions (halo update required)
%   - Analytical solution corresponds to decaying Taylor–Green vortex
%   - Pressure gradient is evaluated directly from analytical expression
%

    %% 1. Memory Allocation

    u_an     = zeros(N+2, N+2);
    v_an     = zeros(N+2, N+2);
    px_an    = zeros(N+2, N+2);
    py_an    = zeros(N+2, N+2);


    %% 2. Temporal Decay Factor

    F = exp(-8*pi^2 * visc_nu * time);


    %% 3. Analytical Field Evaluation (Interior Nodes)

    for j = 2:N+1
        for i = 2:N+1

            % --- 3.1 Velocity field ---
            u_an(i,j) = F * fu(ux(i,j), uy(i,j));
            v_an(i,j) = F * fv(vx(i,j), vy(i,j));

            % --- 3.2 Pressure gradient (analytical) ---
            % p = -rho * F^2 * (cos(2*pi*x)^2 + cos(2*pi*y)^2) / 2
            % dp/dx = rho*pi*F^2*sin(4*pi*x)
            % dp/dy = rho*pi*F^2*sin(4*pi*y)

            px_an(i,j) = rho * pi * F^2 * sin(4*pi * px(i,j));
            py_an(i,j) = rho * pi * F^2 * sin(4*pi * py(i,j));

        end
    end


    %% 4. Halo Update (Periodic Boundary Conditions)

    u_an  = halo_update(u_an);
    v_an  = halo_update(v_an);
    px_an = halo_update(px_an);
    py_an = halo_update(py_an);

end