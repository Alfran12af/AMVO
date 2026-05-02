%% PART C - Time integration and full Navier-Stokes solution (2D)
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
% This script solves the incompressible Navier-Stokes equations in a 2D
% periodic domain by combining the numerical methods developed in Parts A
% and B.
%
% A time integration scheme is implemented using a projection method.
% At each time step:
%   - A predictor velocity is computed from convective and diffusive terms
%   - A Poisson equation is solved to obtain the pseudo-pressure
%   - The velocity field is corrected to enforce incompressibility
%
% The numerical solution is compared with an analytical solution over time.
% The maximum errors of velocity and pressure gradient are evaluated for
% different mesh resolutions to assess the accuracy of the method.
%

clear
close all
clc

%% Mesh definition

cells = [8 16 32];         % Number of control volumes per direction
L = 1;                     % Domain length
t_end = 2;                 % Final simulation time
Re = 100;                  % Reynolds number
rho = 1.225;               % Fluid density

%% Error storage initialization

% Mirar esto, esto no guardara bien el error
mesh_err_u = zeros(1, length(cells));   % Velocity error
mesh_err_p = zeros(1, length(cells));   % Pressure gradient error

%% Symbolic velocity field (manufactured solution)

syms x y

u_sym = cos(2*pi*x)*sin(2*pi*y);
v_sym = -sin(2*pi*x)*cos(2*pi*y);

%% Study point (only for first mesh)

posx = 4;
posy = 4;

first_mesh = true;

%% Mesh loop

for i = 1:length(cells)

    % Current mesh size
    N = cells(i);
    h = L/N;

    
    time_vec = [];
    err_u_time = [];
    err_v_time = [];
    

    %% Temporal plots initialization (first mesh only)

    if i == 1

        % --- Velocity vs time ---
        figure(1); clf;
        hold on; grid on;
        title('Velocity vs time');
        xlabel('t'); ylabel('u, v');

        % --- Pressure gradient vs time ---
        figure(2); clf;
        hold on; grid on;
        title('Pressure gradient vs time');
        xlabel('t'); ylabel('\partial p / \partial x');

    end

    %% Time initialization

    t = 0;

    %% Residual initialization (Adams-Bashforth scheme)

    Ru_prev = zeros(N+2, N+2);
    Rv_prev = zeros(N+2, N+2);

    %% Generate staggered grid coordinates

    [coords.ux, coords.uy, coords.vx, coords.vy, ...
     coords.px, coords.py] = get_coords_2D(L, N);

    %% Initial velocity field

    [u, v] = set_velocity_field(N, u_sym, v_sym, ...
        coords.ux, coords.uy, coords.vx, coords.vy);

    %% Viscosity from Reynolds number

    u_max = max(abs(u(:)));
    visc = u_max * L / Re;

    %% Convert symbolic expressions to numerical functions

    fu = matlabFunction(u_sym, 'Vars', [sym('x'), sym('y')]);
    fv = matlabFunction(v_sym, 'Vars', [sym('x'), sym('y')]);

    %% Time iteration

    while t <= t_end

        %% Time step (CFL condition)

        dt = timestep(N, L, u, v, visc);

        %% Convective and diffusive operators

        [conv_u, conv_v] = convective_2D(u, v, h, N);
        [dif_u, dif_v] = diffusive_2D(N, u, v, h);

        %% Predictor velocity (Adams-Bashforth)

        [u_p, v_p, Ru, Rv] = predictor_velocity( ...
            u, v, conv_u, conv_v, dif_u, dif_v, ...
            h, N, Ru_prev, Rv_prev, dt, visc);

        %% Divergence of predictor velocity

        d = diverg(u_p, v_p, h, N);

        %% Poisson equation assembly

        b = field2vector(d, N);
        
        A = laplacian_matrix(N);

        %% Pressure computation

        % --- Pseudo-pressure ---
        p_pseudo = A \ b;

        % --- Physical pressure ---
        p_real = (rho/dt) * p_pseudo;

        % --- Pressure field ---
        p_field = vector2field(p_real, N);

        % --- Gradient of physical pressure ---
        [px, py] = grad(p_field, h, N);

        % --- Gradient of pseudo-pressure ---
        p_pseudo_field = vector2field(p_pseudo, N);
        [gx, gy] = grad(p_pseudo_field, h, N);

        %% Velocity correction (projection step)

        u_new = u_p - gx;
        v_new = v_p - gy;

        %% Validation: divergence-free condition

        d_new = diverg(u_new, v_new, h, N);
        max_div = max(abs(d_new(:)));

        if max_div > 1e-10
            disp('ERROR: Corrected velocity is not divergence-free');
        end

        %% Analytical solution

        [u_an, v_an, px_an, py_an] = analytic_temporal( ...
            coords.ux, coords.uy, coords.vx, coords.vy, ...
            coords.px, coords.py, t, visc, ...
            fu, fv, rho, N, h);

        %% Error computation

        err_u = max(max(abs(u_new - u_an)));
        err_v = max(max(abs(v_new - v_an)));

        err_px = max(max(abs(px - px_an)));
        err_py = max(max(abs(py - py_an)));

        % Store maximum errors
        mesh_err_u(i) = max(mesh_err_u(i), max(err_u, err_v));

        if t > 0
            mesh_err_p(i) = max(mesh_err_p(i), max(err_px, err_py));
        end

        time_vec(end+1) = t;
        err_u_time(end+1) = err_u;
        err_v_time(end+1) = err_v;

        %% Temporal plots (first mesh only)

        if i == 1

            % --- Velocity comparison ---
            figure(1)
            plot(t, u_an(posx,posy),'bo','MarkerSize',3)
            plot(t, v_an(posx,posy),'ko','MarkerSize',3)
            plot(t, u_new(posx,posy),'rs','MarkerSize',3)
            plot(t, v_new(posx,posy),'gs','MarkerSize',3)

            % --- Pressure gradient comparison ---
            if t > 0
                figure(2)
                plot(t, px_an(posx,posy),'bo','MarkerSize',3)
                plot(t, px(posx,posy),'rs','MarkerSize',3)
            end

        end

        %% Update variables

        u = u_new;
        v = v_new;

        Ru_prev = Ru;
        Rv_prev = Rv;

        t = t + dt;

    end
    
    if i == 1
        err_u_N8 = err_u_time;
        err_v_N8 = err_v_time;
        time_N8 = time_vec;
    elseif i == 2
        err_u_N16 = err_u_time;
        err_v_N16 = err_v_time;
        time_N16 = time_vec;
    elseif i == 3
        err_u_N32 = err_u_time;
        err_v_N32 = err_v_time;
        time_N32 = time_vec;
    end

    %% Legends (only first mesh)

    if i == 1

        figure(1)
        legend('u_{an}','v_{an}','u_{num}','v_{num}')
        hold off

        figure(2)
        legend('px_{an}','px_{num}')
        hold off

    end

end

%% Mesh convergence plot

h_vals = L ./ cells;

figure(3)
loglog(h_vals, mesh_err_u, 'ob-', 'LineWidth',1.5)
hold on
loglog(h_vals, mesh_err_p, 'or-', 'LineWidth',1.5)
loglog(h_vals, h_vals.^2, 'k--', 'LineWidth',1.5)

grid on
xlabel('h')
ylabel('Error')
title('Mesh convergence')
legend('Velocity error','Pressure error','h^2')


figure(4)
hold on
grid on

set(gca,'YScale','log')

title('Error vs time for different N')
xlabel('Time (s)')
ylabel('Abs error')

% N = 8
plot(time_N8, err_u_N8, '-', 'LineWidth',1.5)
plot(time_N8, err_v_N8, '--', 'LineWidth',1.5)

% N = 16
plot(time_N16, err_u_N16, '-', 'LineWidth',1.5)
plot(time_N16, err_v_N16, '--', 'LineWidth',1.5)

% N = 32
plot(time_N32, err_u_N32, '-', 'LineWidth',1.5)
plot(time_N32, err_v_N32, '--', 'LineWidth',1.5)

legend( ...
    'u error, N=8','v error, N=8', ...
    'u error, N=16','v error, N=16', ...
    'u error, N=32','v error, N=32')
