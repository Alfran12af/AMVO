%% Part A
% Authors: Joan Duro | Axl Francia | Pol Jimenez

% 
%
%
%

clear
close all
clc

cells = [8 16 32 64]; % Mesh sizes
L = 1;
h = L./cells;

% Error definition
error_cu = zeros(1, length(cells));
error_cv = zeros(1, length(cells));
error_du = zeros(1, length(cells));
error_dv = zeros(1, length(cells));

% Definition of velocity field
syms x y
u_sym = cos(2*pi*x)*sin(2*pi*y);
v_sym = -sin(2*pi*x)*cos(2*pi*y);
%velocity = [u_sym v_sym];

for i=1:length(cells)
    
    N = cells(i);

    % Mapping
    [coords.ux, coords.uy, coords.vx, coords.vy, coords.px, coords.py] = get_COORDS_2D(L, N);

    % Field generation
    [u_field, v_field] = set_velocity_field(N, u_sym, v_sym, coords.ux, coords.uy, coords.vx, coords.vy);
    
    % Numeric convective form
    [conv_u, conv_v] = convective_2D(u_field, v_field, h(i), L, N);

    % Numeric diffusive form
    [dif_u, dif_v] = diffusive_2D(N, u_field, v_field, h(i));

    % Analytical validation
    [conv_u_an, conv_v_an, diff_u_an, diff_v_an] = analytic_terms(u_sym, v_sym, N, coords.ux, coords.uy, coords.vx, coords.vy);

    % Calculate the error between the numerical and analytical part
    [error_cu(i), error_cv(i), error_du(i), error_dv(i)] = errors_2D(conv_u, conv_v, dif_u, dif_v, conv_u_an, conv_v_an, diff_u_an, diff_v_an);

end

logh = log(h);
log_error_cu = log(error_cu);
log_error_cv = log(error_cv);
log_error_du = log(error_du);
log_error_dv = log(error_dv);



% Convective and diffusive horizontal velocity error
figure;
plot(logh, log_error_cu,'o-', 'LineWidth', 1.5);
hold on;
plot(logh, log_error_du,'o-', 'LineWidth', 1.5);
loglog(logh,log(h.^2), 'LineWidth', 1.5);
xlabel('h [m]'); 
ylabel('error'); 
title('Convective and diffusive error');
legend('error convection u', 'error diffusion u', 'h^2', 'Location', 'NorthWest');

% Convective and diffusive vertical velocity error
figure;
plot(logh,log_error_cv,'o-', 'LineWidth', 1.5);
hold on;
plot(logh,log_error_dv,'o-', 'LineWidth', 1.5);
loglog(logh,log(h.^2), 'LineWidth', 1.5);
xlabel('h [m]'); 
%xlim([1e-4, 1]); 
ylabel('error'); 
%ylim([1e-7, 1e1]); 
title('Convective and diffusive error');
legend('error convection v', 'error diffusion v', 'h^2', 'Location', 'NorthWest');