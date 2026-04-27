%% PART C - Time Integration and Full Navier–Stokes Solution (2D)
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This script solves the incompressible Navier–Stokes equations in a
%   2D periodic domain by combining the numerical methods developed in
%   Parts A and B.
%
%   A projection-based time integration scheme is implemented. At each
%   time step:
%     - A predictor velocity is computed from convective and diffusive terms
%     - A Poisson equation is solved for the pseudo-pressure
%     - The velocity field is corrected to enforce incompressibility
%
%   The numerical solution is compared with an analytical solution over
%   time. The maximum errors of velocity and pressure gradient are evaluated
%   for different mesh resolutions to assess the accuracy of the method.
%

clear
close all
clc


% 1. Mesh Definition

cells = [8 16 32];         % Number of control volumes per direction
L     = 1;                 % Domain length
t_end = 2;                 % Final simulation time
Re    = 100;               % Reynolds number
rho   = 1.225;             % Fluid density


% 2. Error Storage Initialization

mesh_err_u = zeros(1, length(cells));   % Velocity error


% 3. Symbolic Velocity Field (Manufactured Solution)

syms x y

u_sym =  cos(2*pi*x) * sin(2*pi*y);
v_sym = -sin(2*pi*x) * cos(2*pi*y);


% 4. Study Point (First Mesh Only)

posx = 4;
posy = 4;

first_mesh = true;
iterations = zeros(1, length(cells));



% 5. Mesh Loop

for i = 1:length(cells)

    % 5.1 Mesh Parameters
    
    N = cells(i);
    h = L / N;

    time_vec   = [];
    err_u_time = [];
    err_v_time = [];


    % 5.2 Temporal Plots Initialization (First Mesh Only)

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


    % 5.3 Time Initialization

    t = 0;
    iter = 0;


    % 5.4 Residual Initialization (Adams–Bashforth)

    Ru_prev = zeros(N+2, N+2);
    Rv_prev = zeros(N+2, N+2);


    % 5.5 Staggered Grid Coordinates

    [coords.ux, coords.uy, ...
     coords.vx, coords.vy, ...
     coords.px, coords.py] = get_coords_2D(L, N);


    % 5.6 Initial Velocity Field

    [u, v] = set_velocity_field( ...
        N, u_sym, v_sym, ...
        coords.ux, coords.uy, ...
        coords.vx, coords.vy);


    % 5.7 Viscosity from Reynolds Number

    u_max = max(abs(u(:)));
    visc  = u_max * L / Re;


    % 5.8 Symbolic → Numerical Functions

    fu = matlabFunction(u_sym, 'Vars', [sym('x'), sym('y')]);
    fv = matlabFunction(v_sym, 'Vars', [sym('x'), sym('y')]);



    % Animation set up
    max_steps = 5000;
    U_hist = zeros(N, N, max_steps); 
    V_hist = zeros(N, N, max_steps);
    P_hist = zeros(N, N, max_steps);
    T_hist = zeros(1, max_steps);
    k_save = 1;


    % 5.9 Time Iteration

    while t <= t_end

        iter = iter + 1;

        % 5.9.1 Time Step (CFL Condition)

        dt    = timestep(N, L, u, v, visc);
        t_new = t + dt;


        % 5.9.2 Convective and Diffusive Operators

        [conv_u, conv_v] = convective_2D(u, v, h, N);
        [dif_u,  dif_v]  = diffusive_2D(N, u, v, h);


        % 5.9.3 Predictor Velocity (Adams–Bashforth)

        [u_p, v_p, Ru, Rv] = predictor_velocity( ...
            u, v, ...
            conv_u, conv_v, ...
            dif_u,  dif_v, ...
            h, N, ...
            Ru_prev, Rv_prev, ...
            dt, visc);


        % 5.9.4 Divergence of Predictor Velocity

        d = diverg(u_p, v_p, h, N);


        % 5.9.5 Poisson Equation Assembly

        b = field2vector(d, N);
        A = laplacian_matrix(N);


        % 5.9.6 Pressure Computation

        % Pseudo-pressure
        p_pseudo = A \ b;

        % Physical pressure
        p_real = (rho / dt) * p_pseudo;

        % Pressure field
        p_field = vector2field(p_real, N);

        % Gradient of physical pressure
        [px, py] = grad(p_field, h, N);

        % Gradient of pseudo-pressure
        p_pseudo_field = vector2field(p_pseudo, N);
        [gx, gy] = grad(p_pseudo_field, h, N);


        % 5.9.7 Velocity Correction (Projection Step)

        u_new = u_p - gx;
        v_new = v_p - gy;


        % 5.9.8 Validation: Divergence-Free Condition

        d_new   = diverg(u_new, v_new, h, N);
        max_div = max(abs(d_new(:)));

        if max_div > 1e-10
            disp('ERROR: Corrected velocity is not divergence-free');
        end


        % 5.9.9 Analytical Solution

        [u_an, v_an, px_an, py_an] = analytic_temporal( ...
            coords.ux, coords.uy, ...
            coords.vx, coords.vy, ...
            coords.px, coords.py, ...
            t_new, visc, ...
            fu, fv, rho, ...
            N, h);


        % 5.9.10 Error Computation

        err_u = max(max(abs(u_new(2:N+1,2:N+1) - u_an(2:N+1,2:N+1))));
        err_v = max(max(abs(v_new(2:N+1,2:N+1) - v_an(2:N+1,2:N+1))));

        err_px = max(max(abs(px(2:N+1,2:N+1) - px_an(2:N+1,2:N+1))));
        err_py = max(max(abs(py(2:N+1,2:N+1) - py_an(2:N+1,2:N+1))));

        % Store maximum error
        mesh_err_u(i) = max(mesh_err_u(i), max(err_u, err_v));

        time_vec(end+1)   = t;
        err_u_time(end+1) = err_u;
        err_v_time(end+1) = err_v;


        % 5.9.11 Temporal Plots (First Mesh Only)

        if i == 1

            % --- Velocity ---
            figure(1)
            plot(t, u_an(posx,posy),'bo','MarkerSize',3)
            plot(t, v_an(posx,posy),'ko','MarkerSize',3)
            plot(t, u_new(posx,posy),'rs','MarkerSize',3)
            plot(t, v_new(posx,posy),'gs','MarkerSize',3)

            % --- Pressure gradient ---
            if t > 0
                figure(2)
                plot(t, px_an(posx,posy),'bo','MarkerSize',3)
                plot(t, px(posx,posy),'rs','MarkerSize',3)
            end

        end


        % Animation set up
        if i==length(cells)
            if mod(iter, 5) == 0 && k_save <= max_steps
                U_hist(:,:,k_save) = u_new(2:N+1, 2:N+1);
                V_hist(:,:,k_save) = v_new(2:N+1, 2:N+1);
                P_hist(:,:,k_save) = p_field(2:N+1, 2:N+1);
                T_hist(k_save)     = t;
                k_save = k_save + 1;
            end
        end


        % 5.9.12 Update Variables

        u = u_new;
        v = v_new;

        Ru_prev = Ru;
        Rv_prev = Rv;

        t = t + dt;

    end

    iterations(i) = iter;


    % 5.10 Store Time Histories

    if i == 1
        err_u_N8  = err_u_time;
        err_v_N8  = err_v_time;
        time_N8   = time_vec;
    elseif i == 2
        err_u_N16 = err_u_time;
        err_v_N16 = err_v_time;
        time_N16  = time_vec;
    elseif i == 3
        err_u_N32 = err_u_time;
        err_v_N32 = err_v_time;
        time_N32  = time_vec;
    end


    % 5.11 Legends (First Mesh Only)

    if i == 1

        figure(1)
        legend('u_{an}','v_{an}','u_{num}','v_{num}')
        hold off

        figure(2)
        legend('px_{an}','px_{num}')
        hold off

    end

end


% 6. Mesh Convergence Plot

h_vals = L ./ cells;

figure(3)
loglog(h_vals, mesh_err_u, 'ob-', 'LineWidth',1.5)
hold on
loglog(h_vals, h_vals.^2, 'k--', 'LineWidth',1.5)

grid on
xlabel('h')
ylabel('Error')
title('Mesh convergence')
legend('Velocity error','h^2')


% 7. Error vs Time

figure(4)
hold on
grid on

set(gca,'YScale','log')

title('Error vs time for different N')
xlabel('Time (s)')
ylabel('Abs error')

% N = 8
plot(time_N8,  err_u_N8,  '-',  'LineWidth',1.5)
plot(time_N8,  err_v_N8,  '--', 'LineWidth',1.5)

% N = 16
plot(time_N16, err_u_N16, '-',  'LineWidth',1.5)
plot(time_N16, err_v_N16, '--', 'LineWidth',1.5)

% N = 32
plot(time_N32, err_u_N32, '-',  'LineWidth',1.5)
plot(time_N32, err_v_N32, '--', 'LineWidth',1.5)

legend( ...
    'u error, N=8','v error, N=8', ...
    'u error, N=16','v error, N=16', ...
    'u error, N=32','v error, N=32')



%% Animations
% Centered mesh values

x_int = coords.px(2:N+1, 2:N+1);
y_int = coords.py(2:N+1, 2:N+1);


U_hist = U_hist(:,:, 1:k_save-1);
V_hist = V_hist(:,:, 1:k_save-1);
P_hist = P_hist(:,:, 1:k_save-1);
T_hist = T_hist(1:k_save-1);

% Y ahora lanzas la magia:
animate_NS(x_int, y_int, U_hist, V_hist, P_hist, T_hist, rho);
