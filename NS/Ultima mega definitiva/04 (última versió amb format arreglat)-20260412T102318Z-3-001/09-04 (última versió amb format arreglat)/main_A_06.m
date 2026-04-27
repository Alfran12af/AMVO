%% PART A - Verification of Convective and Diffusive Terms (2D)
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This script verifies the implementation of the convective and diffusive
%   terms of the 2D incompressible Navier–Stokes equations using the
%   Method of Manufactured Solutions (MMS).
%
%   A manufactured analytical velocity field is defined symbolically and
%   evaluated on a staggered mesh. The numerical discretization of the
%   convective and diffusive operators is computed using a finite
%   control-volume (FCV) formulation.
%
%   Analytical expressions of the operators are obtained via symbolic
%   differentiation and compared with the numerical results. The error is
%   evaluated for several mesh resolutions to verify the expected
%   second-order convergence of the scheme.
%

clear
close all
clc


%% 1. Mesh Definition

cells = [8 16 32 64];      % Number of control volumes per direction
L     = 1;                 % Domain length
h     = L ./ cells;        % Grid spacing


%% 2. Error Storage Initialization

error_cu = zeros(1, length(cells));   % Convective error (u component)
error_cv = zeros(1, length(cells));   % Convective error (v component)
error_du = zeros(1, length(cells));   % Diffusive error (u component)
error_dv = zeros(1, length(cells));   % Diffusive error (v component)


%% 3. Symbolic Velocity Field

syms x y

% Manufactured analytical velocity field (divergence-free)
u_sym =  cos(2*pi*x) * sin(2*pi*y);
v_sym = -sin(2*pi*x) * cos(2*pi*y);


%% 4. Mesh Loop

for i = 1:length(cells)
    
    % Current mesh size
    N = cells(i);

    
    %% 4.1 Staggered Grid Coordinates
    
    % Coordinates for u, v and pressure nodes (Arakawa-C grid)
    [coords.ux, coords.uy, ...
     coords.vx, coords.vy, ...
     coords.px, coords.py] = get_coords_2D(L, N);


    %% 4.2 Velocity Field Initialization
    
    % Evaluate analytical solution on staggered grid
    [u_field, v_field] = set_velocity_field( ...
        N, u_sym, v_sym, ...
        coords.ux, coords.uy, ...
        coords.vx, coords.vy);


    %% 4.3 Numerical Convective Terms
    
    % Discrete convective operator (FCV formulation)
    [conv_u, conv_v] = convective_2D(u_field, v_field, h(i), N);


    %% 4.4 Numerical Diffusive Terms
    
    % Discrete diffusive operator (Laplacian)
    [dif_u, dif_v] = diffusive_2D(N, u_field, v_field, h(i));


    %% 4.5 Analytical Operators
    
    % Analytical convective and diffusive terms
    [conv_u_an, conv_v_an, diff_u_an, diff_v_an] = ...
        analytic_terms(u_sym, v_sym, ...
        N, ...
        coords.ux, coords.uy, ...
        coords.vx, coords.vy);


    %% 4.6 Error Computation
    
    % Error between numerical and analytical operators
    [error_cu(i), error_cv(i), ...
     error_du(i), error_dv(i)] = ...
        errors_2D(N, ...
        conv_u, conv_v, ...
        dif_u, dif_v, ...
        conv_u_an, conv_v_an, ...
        diff_u_an, diff_v_an);

end


%% 5. Logarithmic Values for Convergence Analysis

logh = log(h);

log_error_cu = log(error_cu);
log_error_cv = log(error_cv);
log_error_du = log(error_du);
log_error_dv = log(error_dv);


%% 6. Convergence Plot - Horizontal Velocity (u)

figure

loglog(h, error_cu, 'o-', 'LineWidth', 1.5)
hold on
loglog(h, error_du, 'o-', 'LineWidth', 1.5)

% Reference second-order slope
loglog(h, h.^2, 'LineWidth', 1.5)

xlim([h(end) h(1)])

xlabel('h [m]')
ylabel('error')

title('Convective and Diffusive Error (u velocity)')

legend('error convection u', ...
       'error diffusion u', ...
       'h^2', ...
       'Location', 'NorthWest')


%% 7. Convergence Plot - Vertical Velocity (v)

figure

loglog(h, error_cv, 'o-', 'LineWidth', 1.5)
hold on
loglog(h, error_dv, 'o-', 'LineWidth', 1.5)

% Reference second-order slope
loglog(h, h.^2, 'LineWidth', 1.5)

xlim([h(end) h(1)])

xlabel('h [m]')
ylabel('error')

title('Convective and Diffusive Error (v velocity)')

legend('error convection v', ...
       'error diffusion v', ...
       'h^2', ...
       'Location', 'NorthWest')