%% Part A
% Authors: Joan Duro | Axl Francia | Pol Jimenez

% 
%
%
%

clear
close all
clc

cells = [10]; % Mesh sizes
L = 1;
h = L/cells;

% Definition of velocity field
syms x y
u_sym = x;%cos(2*pi*x)*sin(2*pi*y);
v_sym = x;%-sin(2*pi*x)*cos(2*pi*y);
%velocity = [u_sym v_sym];


for i=1:length(cells)
    N = cells(i);

    [coords.ux, coords.uy, coords.vx, coords.vy, coords.px, coords.py] = get_COORDS_2D(L, N);

    % Field generation
    [u_field, v_field] = set_velocity_field(N, u_sym, v_sym, coords.ux, coords.uy, coords.vx, coords.vy);


    % Numeric convective form
    % cu = convective_2D(u_field, v_field, L)

    % Numeric diffusive form
    [dif_u, dif_v] = diffusive_2D(N, u_field, v_field);



    % Analytical validation





end
