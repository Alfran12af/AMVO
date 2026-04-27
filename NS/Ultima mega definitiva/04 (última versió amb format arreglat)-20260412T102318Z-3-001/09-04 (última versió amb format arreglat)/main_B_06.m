%% PART B - Pressure-Velocity Coupling (2D)
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This script implements and verifies the pressure–velocity coupling
%   for the 2D incompressible Navier–Stokes equations using a projection
%   method.
%
%   An arbitrary intermediate (provisional) velocity field is defined on
%   a staggered grid. Its divergence is computed to quantify the lack of
%   mass conservation at the discrete level.
%
%   A Poisson equation for the pseudo-pressure is assembled and solved.
%   The resulting pressure field is used to correct the velocity field by
%   subtracting its gradient, enforcing a divergence-free condition.
%
%   Verification checks:
%     - Divergence of provisional velocity ≠ 0
%     - Global conservation of RHS (sum(b) = 0)
%     - Corrected velocity is divergence-free
%

clear
close all
clc


%% 1. Mesh Definition

N = 8;                      % Number of control volumes per direction
L = 1;                      % Domain length
h = L / N;                  % Grid spacing


%% 2. Physical Parameters

rho = 1.225;                % Fluid density
dt  = 0.01;                 % Time step (used in Part C)


%% 3. Initial Velocity Field (Provisional)

% Initialize velocity fields (including halo nodes)
up = zeros(N+2, N+2);
vp = zeros(N+2, N+2);

% Arbitrary perturbation
up(3,3) = 1;
% vp(3,3) = 1;

% Apply halo update (periodic boundary conditions)
up = halo_update(up);
vp = halo_update(vp);


%% 4. Staggered Grid Coordinates

% Coordinates for u, v and pressure nodes (Arakawa-C grid)
[coords.ux, coords.uy, ...
 coords.vx, coords.vy, ...
 coords.px, coords.py] = get_coords_2D(L, N);


%% 5. Divergence of Provisional Velocity

% Discrete divergence (integral form)
d = diverg(up, vp, h, N);


%% 6. Right-Hand Side of Poisson Equation

% Convert divergence field into vector form
b = field2vector(d, N);


%% 7. Laplacian Matrix Assembly

% Discrete Laplacian operator (Poisson system)
A = laplacian_matrix(N);


%% 8. Poisson Equation Solution

% Solve for pseudo-pressure
p_pseudo = A \ b;


%% 9. Pressure Field Reconstruction

% Convert vector to scalar field (centered grid)
pfield = vector2field(p_pseudo, N);


%% 10. Pressure Gradient Computation

% Gradient on staggered grid
[gx, gy] = grad(pfield, h, N);


%% 11. Velocity Correction

% Enforce incompressibility
ufut = up - gx;
vfut = vp - gy;


%% 12. Validation

% --- 12.1 Divergence of provisional velocity (should NOT be zero) ---
d_up = diverg(up, vp, h, N);
max_div_up = max(abs(d_up(:)));

fprintf('Max divergence of provisional velocity (up,vp): %e\n', max_div_up);

if max_div_up == 0
    disp('ERROR: Divergence of provisional velocity is zero');
end


% --- 12.2 Global mass conservation (sum of RHS must be zero) ---
sum_b = sum(b);

fprintf('Sum of right hand side vector (b): %e\n', sum_b);

if abs(sum_b) > 1e-12
    disp('ERROR: Sum of b is not zero (global mass conservation violated)');
end


% --- 12.3 Divergence of corrected velocity (should be ~ 0) ---
d_ufut = diverg(ufut, vfut, h, N);
max_div_ufut = max(abs(d_ufut(:)));

fprintf('Max divergence of corrected velocity (ufut,vfut): %e\n', max_div_ufut);

if max_div_ufut > 1e-10
    disp('ERROR: Corrected velocity is not divergence-free');
end


%% 13. Plots

%% 13.1 Provisional Velocity Field

% Visualization of the initial (uncorrected) velocity field
figure(1);
hold on

% Horizontal component (u) - red
quiver(coords.ux, coords.uy, up, zeros(size(up)), 'r');

% Vertical component (v) - blue
quiver(coords.vx, coords.vy, zeros(size(vp)), vp, 'b');

title('Distribution of the Provisional Velocity Field');
xlabel('X');
ylabel('Y');

% Grid lines (staggered mesh visualization)
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


%% 13.2 Corrected Velocity Field

% Visualization after projection step
figure(2);
hold on

quiver(coords.ux, coords.uy, ufut, zeros(size(ufut)), 'r');
quiver(coords.vx, coords.vy, zeros(size(vfut)), vfut, 'b');

title('Distribution of the Corrected Velocity Field (n+1)');
xlabel('X');
ylabel('Y');

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


%% 13.3 Pseudo-Pressure Field

figure(3);
clf
hold on

% Extract interior pressure values
p_plot = pfield(2:N+1, 2:N+1);

% Handle coordinates (with/without halo)
if isequal(size(coords.px), size(pfield))
    px_plot = coords.px(2:N+1, 2:N+1);
    py_plot = coords.py(2:N+1, 2:N+1);
else
    px_plot = coords.px;
    py_plot = coords.py;
end

% Scatter plot
marker_size = 140;

scatter(px_plot(:), py_plot(:), marker_size, p_plot(:), ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);

colorbar;
colormap(jet);

title('Distribution of the Pseudo-Pressure Field');
xlabel('X');
ylabel('Y');

% Grid
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


%% 14. Execution Time

toc