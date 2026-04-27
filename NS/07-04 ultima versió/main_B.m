%% PART B - Pressure-velocity coupling (2D)
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
% This script implements and verifies the pressure-velocity coupling for
% the 2D incompressible Navier-Stokes equations using a projection method.
%
% An arbitrary intermediate velocity field is first defined on a staggered
% grid. Its divergence is computed in order to quantify the lack of mass
% conservation at the discrete level.
%
% A Poisson equation for the pseudo-pressure is then assembled and solved.
% The solution of this equation allows correcting the velocity field by
% subtracting the gradient of the pseudo-pressure, enforcing a divergence-
% free velocity field.
%
% The correctness of the implementation is verified by checking that:
%   - The divergence of the provisional velocity is non-zero
%   - The right-hand side of the Poisson equation satisfies global
%     conservation (sum equals zero)
%   - The corrected velocity field is divergence-free
%

clear
close all
clc

%% Mesh definition

N = 8;                      % Number of control volumes per direction
L = 1;                      % Domain length
h = L/N;                    % Grid spacing

%% Physical parameters

rho = 1.225;                % Fluid density
dt = 0.01;                  % Time step (used in Part C)

%% Initial velocity field (provisional)

% Initialize velocity fields (including halo)
up = zeros(N+2, N+2);
vp = zeros(N+2, N+2);

% Assign arbitrary perturbation
up(3,3) = 1;
%vp(3,3) = 1;

% Apply halo update (periodic boundary conditions)
up = halo_update(up);
vp = halo_update(vp);

    %% Generate staggered grid coordinates

    % Coordinates of the staggered mesh for u, v and pressure nodes
    [coords.ux, coords.uy, coords.vx, coords.vy, coords.px, coords.py] = get_coords_2D(L, N);

%% Divergence of provisional velocity

% Compute divergence of the intermediate velocity field
d = diverg(up, vp, h, N);

%% Right-hand side of Poisson equation

% Convert divergence field into vector form
b = field2vector(d, N);

%% Laplacian matrix assembly

% Build discrete Laplacian operator
A = laplacian_matrix(N);

%% Poisson equation solution

% Solve for pseudo-pressure
p_pseudo = A \ b;

%% Pressure field reconstruction

% Convert pressure vector back to field form
pfield = vector2field(p_pseudo, N);

%% Pressure gradient computation

% Compute gradient of pseudo-pressure
[gx, gy] = grad(pfield, h, N);

%% Velocity correction

% Correct velocity field to enforce incompressibility
ufut = up - gx;
vfut = vp - gy;

%% Validation

% --- Divergence of provisional velocity (should NOT be zero) ---
d_up = diverg(up, vp, h, N);
max_div_up = max(abs(d_up(:)));

fprintf('Max divergence of provisional velocity (up,vp): %e\n', max_div_up);

if max_div_up == 0
    disp('ERROR: Divergence of provisional velocity is zero');
end

% --- Global mass conservation (sum of right hand side must be zero) ---
sum_b = sum(b);

fprintf('Sum of right hand side vector (b): %e\n', sum_b);

if abs(sum_b) > 1e-12
    disp('ERROR: Sum of b is not zero (global mass conservation violated)');
end

% --- Divergence of corrected velocity (should be approx 0) ---
d_ufut = diverg(ufut, vfut, h, N);
max_div_ufut = max(abs(d_ufut(:)));

fprintf('Max divergence of corrected velocity (ufut,vfut): %e\n', max_div_ufut);

if max_div_ufut > 1e-10
    disp('ERROR: Corrected velocity is not divergence-free');
end

%% Plots

% --- Figure 1: Provisional velocity field ------------------------------
% This figure displays the arbitrary intermediate velocity field used
% as input for the pressure-correction step. The horizontal component
% (u) is represented in red on the u-staggered nodes, while the vertical
% component (v) is represented in blue on the v-staggered nodes.
% The plot provides a visual interpretation of the imposed perturbation
% before applying the incompressibility correction.
figure(1);
hold on
% Plot of the horizontal velocity component u on its staggered locations
quiver(coords.ux, coords.uy, up, zeros(size(up)), 'r');
% Plot of the vertical velocity component v on its staggered locations
quiver(coords.vx, coords.vy, zeros(size(vp)), vp, 'b');
title('Distribution of the provisional velocity field');
xlabel('X');
ylabel('Y');
% Draw the grid lines to highlight the staggered arrangement
for xg = 0:h:L
    xline(xg, '-k');
end
for yg = 0:h:L
    yline(yg, '-k');
end
axis equal;
xlim([0 L]);
ylim([0 L]);
hold off


% --- Figure 2: Corrected velocity field --------------------------------
% This figure shows the velocity field after applying the pressure
% correction. The corrected field should satisfy the discrete mass
% conservation condition, meaning that its divergence is expected to
% be approximately zero. Comparing this figure with Figure 1 allows
% visualizing the effect of the projection step on the velocity field.
figure(2);
hold on
% Plot of the corrected horizontal velocity component
quiver(coords.ux, coords.uy, ufut, zeros(size(ufut)), 'r');
% Plot of the corrected vertical velocity component
quiver(coords.vx, coords.vy, zeros(size(vfut)), vfut, 'b');
title('Distribution of the corrected velocity field at n+1');
xlabel('X');
ylabel('Y');
% Draw the grid lines for reference
for xg = 0:h:L
    xline(xg, '-k');
end
for yg = 0:h:L
    yline(yg, '-k');
end
axis equal;
xlim([0 L]);
ylim([0 L]);
hold off


% --- Figure 3: Pseudo-pressure field -----------------------------------
% This figure represents the pseudo-pressure distribution obtained from
% the solution of the Poisson equation. The pressure values are plotted
% at the pressure-control-volume centers, which correspond to the cell
% centers of the staggered grid. Marker color indicates the magnitude
% of the pseudo-pressure, allowing identification of the regions where
% the pressure correction is strongest.
figure(3);
clf
hold on
% Extract pressure values only at the interior pressure nodes
p_plot = pfield(2:N+1, 2:N+1);
% Use the pressure-node coordinates. If the coordinate arrays include
% halo nodes, keep only the interior values; otherwise, use them directly.
if isequal(size(coords.px), size(pfield))
    px_plot = coords.px(2:N+1, 2:N+1);
    py_plot = coords.py(2:N+1, 2:N+1);
else
    px_plot = coords.px;
    py_plot = coords.py;
end
% Use a fixed marker size so that all pressure points remain clearly visible
marker_size = 140;
% Scatter plot of pseudo-pressure at cell centers
scatter(px_plot(:), py_plot(:), marker_size, p_plot(:), ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
colorbar;
colormap(jet);
title('Distribution of the pseudo-pressure field');
xlabel('X');
ylabel('Y');
% Draw the control-volume grid
for xg = 0:h:L
    xline(xg, '-k');
end
for yg = 0:h:L
    yline(yg, '-k');
end

axis equal;
xlim([0 L]);
ylim([0 L]);
hold off

toc