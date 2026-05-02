%% PART B - Pressure-velocity coupling (2D)
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
% This script implements the pressure-velocity coupling for the 2D
% incompressible Navier-Stokes equations using a projection method.
%
% An intermediate velocity field is first defined and its divergence is
% computed. A Poisson equation for the pseudo-pressure is then assembled
% and solved in order to enforce mass conservation.
%
% The pressure gradient is used to correct the intermediate velocity field,
% obtaining a divergence-free velocity at the new time step.
%
% This script serves as a basic validation of the pressure correction step
% and the consistency of the discrete operators involved.
%

clear
close all
clc

%% Constants

N = 8;                      % Number of control volumes per direction
L = 1;                      % Domain length
h = L/N;                    % Grid spacing
rho = 1.225;                % Fluid density
dt = 0.01;                  % Time step

%% Initial velocity field

% Initialize provisional velocity field
up = zeros(N+2);
vp = zeros(N+2);

% Assign arbitrary non-zero values
up(3,3) = 1;
vp(3,3) = 1;

% Apply halo update
up = halo_update(up);
vp = halo_update(vp);

%% Divergence of provisional velocity

% Compute divergence of the intermediate field
d = diverg(up, vp, h, N);

%% Right-hand side of Poisson equation

% Convert divergence field into vector form
b = field2vector(d, N);

%% Laplacian matrix assembly

% Build discrete Laplacian operator
A = laplacian_matrix(N);

%% Poisson equation solution

% Solve for pseudo-pressure
p = A \ b;

%% Pressure field reconstruction

% Convert pressure vector back to field form
pfield = vector2field(p, N);

%% Pressure gradient computation

% Compute gradient of pseudo-pressure
[gx, gy] = grad(pfield, h, N);

%% Velocity correction

% Enforce incompressibility
ufut = up - gx;
vfut = vp - gy;

%% Validation

% Divergence of provisional velocity (should NOT be zero)
d_up = diverg(up, vp, h, N);
max_div_up = max(max(abs(d_up)));

fprintf('Max divergence of provisional velocity (up,vp): %e\n', max_div_up);

if max_div_up == 0
    disp('ERROR: Divergence of provisional velocity is zero');
end

% Global mass conservation (sum of RHS must be zero)
sum_b = sum(b);

fprintf('Sum of RHS vector (b): %e\n', sum_b);

if abs(sum_b) > 1e-12
    disp('ERROR: Sum of b is not zero (global mass conservation violated)');
end

% Divergence of corrected velocity (should be ~0)
d_ufut = diverg(ufut, vfut, h, N);
max_div_ufut = max(max(abs(d_ufut)));

fprintf('Max divergence of corrected velocity (ufut,vfut): %e\n', max_div_ufut);

if max_div_ufut > 1e-10
    disp('ERROR: Corrected velocity is not divergence-free');
end

%% Plots

% 