%% PART A - Verification of convective and diffusive terms (2D)
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
% This script verifies the implementation of the convective and diffusive
% terms of the 2D Navier-Stokes equations using the Method of Manufactured
% Solutions (MMS).
%
% A manufactured analytical velocity field is defined symbolically and
% evaluated on the numerical mesh. The numerical discretization of the
% convective and diffusive operators is then computed using a finite
% control-volume formulation.
%
% Analytical expressions of the operators are obtained through symbolic
% differentiation and compared with the numerical results. The error is
% evaluated for several mesh resolutions in order to verify the expected
% second-order convergence of the scheme.
%

clear
close all
clc

%% Mesh definition

cells = [8 16 32 64];      % Number of control volumes per direction
L = 1;                     % Domain length
h = L./cells;              % Grid spacing for each mesh

%% Error storage initialization

error_cu = zeros(1, length(cells));   % Convective error (u component)
error_cv = zeros(1, length(cells));   % Convective error (v component)
error_du = zeros(1, length(cells));   % Diffusive error (u component)
error_dv = zeros(1, length(cells));   % Diffusive error (v component)

%% Symbolic velocity field

syms x y

% Analytical velocity field
u_sym = cos(2*pi*x)*sin(2*pi*y);
v_sym = -sin(2*pi*x)*cos(2*pi*y);

%% Mesh loop

for i = 1:length(cells)
    
    % Current mesh size
    N = cells(i);

    %% Generate staggered grid coordinates

    % Coordinates of the staggered mesh for u, v and pressure nodes
    [coords.ux, coords.uy, coords.vx, coords.vy, coords.px, coords.py] = get_coords_2D(L, N);

    %% Generate velocity fields from analytical solution

    % Evaluate the symbolic velocity field on the staggered grid
    [u_field, v_field] = set_velocity_field(N, u_sym, v_sym, coords.ux, coords.uy, coords.vx, coords.vy);

    %% Numerical convective terms

    % Compute discrete convective operator
    [conv_u, conv_v] = convective_2D(u_field, v_field, h(i), N);

    %% Numerical diffusive terms

    % Compute discrete Laplacian operator
    [dif_u, dif_v] = diffusive_2D(N, u_field, v_field, h(i));

    %% Analytical operators

    % Compute analytical convective and diffusive terms
    [conv_u_an, conv_v_an, diff_u_an, diff_v_an] = analytic_terms(u_sym, v_sym, N, coords.ux, coords.uy, coords.vx, coords.vy);

    %% Error computation

    % Compute error between numerical and analytical operators
    [error_cu(i), error_cv(i), error_du(i), error_dv(i)] = errors_2D(N, conv_u, conv_v, dif_u, dif_v, conv_u_an, conv_v_an, diff_u_an, diff_v_an);

end

%% Logarithmic values for convergence analysis

logh = log(h);

log_error_cu = log(error_cu);
log_error_cv = log(error_cv);

log_error_du = log(error_du);
log_error_dv = log(error_dv);

%% Convergence plot for the horizontal velocity

figure

loglog(h, error_cu,'o-','LineWidth',1.5)
hold on
loglog(h, error_du,'o-','LineWidth',1.5)

% Reference second-order slope
loglog(h, h.^2,'LineWidth',1.5)

xlim([h(end) h(1)])

xlabel('h [m]')
ylabel('error')

title('Convective and diffusive error (u velocity)')

legend('error convection u', ...
       'error diffusion u', ...
       'h^2', ...
       'Location','NorthWest')

%% Convergence plot for the vertical velocity

figure

loglog(h, error_cv,'o-','LineWidth',1.5)
hold on
loglog(h, error_dv,'o-','LineWidth',1.5)

% Reference second-order slope
loglog(h, h.^2,'LineWidth',1.5)

xlim([h(end) h(1)])

xlabel('h [m]')
ylabel('error')

title('Convective and diffusive error (v velocity)')

legend('error convection v', ...
       'error diffusion v', ...
       'h^2', ...
       'Location','NorthWest')