%% PART C - (2D)
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%

clear
close all
clc

syms x y

%% Constants

cells = [8 16 32 64];       % Number of control volumes per direction
L = 1;                      % Domain length
Re = 100;                   % Viscosity will be calculated as Re = 100 
time_end = 3;               % Total time to test
rho = 1.225;                % Fluid density

%% Initializations

h = zeros (1,length(cells));    % Espai entre nodes
vp = zeros(N+2);                

%% Loops
% Loop d'elements
for i = 1:length(cells)
    time = 0;                      % Inicialització del temps
    N = cells(i);
    h(i) = L/N;
    Ruanterior = zeros(N+2);
    Rvanterior = zeros(N+2);

    % Coordinates of the staggered mesh for u, v and pressure nodes
    [coords.ux, coords.uy, coords.vx, coords.vy, coords.px, coords.py] = get_coords_2D(L, N);

    % Inicialitzem el camp velocitats per t = 0
    
    u_sym = cos(2*pi*x)*sin(2*pi*y);
    v_sym = -sin(2*pi*x)*cos(2*pi*y);

    [u_field, v_field] = set_velocity_field(N, u_sym, v_sym, coords.ux, coords.uy, coords.vx, coords.vy);
    
    % Càlcul de la viscositat per el primer valor màxim de velocitat de cada malla
    visc_nu = max(max(abs(u_field(:))))*L/Re;

    % Loop de temps
    while time <= time_end

        % Càlcul del timestep
        [dt] = timestep(N, L, u_field, v_field, visc_nu);

        % Compute discrete convective operator
        [conv_u, conv_v] = convective_2D(u_field, v_field, h(i), N);

        % Compute discrete Laplacian operator
        [dif_u, dif_v] = diffusive_2D(N, u_field, v_field, h(i));
    

        % Calculation of the predictor velocity
        [u_p, v_p, Ru, Rv] = predictor_velocity(u_field, v_field, conv_u, conv_v, dif_u, dif_v, h(i), Ruanterior, Rvanterior, dt, visc_nu); 

        % Compute divergence of the intermediate field and Validation
        d = diverg(u_p, v_p, h(i), N); % The divergence of u_p is not zero
        if d(:) < 10^(-5)
            disp('ERROR div(u_p) = 0'); % It may be zero at some nodes
        end
        
        % Build discrete Laplacian operator
        A = laplacian_matrix(N);
        
        % Convert divergence field into vector form
        b = field2vector(d, N);
        
        % Solve for pseudo-pressure
        p = A \ b;
        
        % Convert pressure vector back to field form
        p_pseudo = vector2field(p, N);

        % Obtain real pressure
        p_real = p*rho/dt;
        p_field = vector2field(p_real, N);
        
        % Compute gradient of pseudo-pressure and pressure
        [px, py] = grad(p_field, h(i), N);
        [gx, gy] = grad(p_pseudo, h(i), N);

        
        % Enforce incompressibility
        u_fut = u_p - gx;
        v_fut = v_p - gy;





        % Validation
        [d] = diverg(u_fut, v_fut, h(i), N);
        if d(:) > 10^(-5)
            disp('ERROR div(u_f) should be 0');
        end
        
        % Calculation of the analytic pressures and velocities  
        [u_an, v_an, px_an, py_an] = analytic_temporal(coords.ux, coords.uy, coords.vx, coords.vy, coords.px, coords.py, time, visc_nu, u_sym, v_sym, rho, N, h(i)); 
        
        
        % Update of variables 
        u = u_fut;
        v = v_fut;
        Ruanterior = Ru;
        Rvanterior = Rv;
        
        % Store the maximum error from the current time step
        err_u = abs(u - u_an);
        err_v = abs(v - v_an);
        err_px = abs(px - px_an);
        err_py = abs(py - py_an);
        max_erru = max(max(err_u(:)), max(err_v(:)));
        max_errp = max(max(err_px(:)), max(err_py(:)));
        
        % Store the maximum error from the current mesh size
        if max_erru > mesh_erru(i)
            mesh_erru(i) = max_erru;
        end
        
        if max_errp > mesh_errp(i) && t~=0
            mesh_errp(i) = max_errp;
        end
        
        
        % Temporal plot (first mesh only)
        if i == 1
            figure (1)
            plot(t, u_an(posx_study, posy_study),'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerSize', 3);  
            plot(t, v_an(posx_study, posy_study), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 3);
            plot(t, u(posx_study,posy_study), 's', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 3);
            plot(t, v(posx_study,posy_study), 's', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'MarkerSize', 3);
        if t ~= 0
            figure (5)
            plot(t, px_an(posx_study, posy_study),'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerSize', 3);  
            plot(t, px(posx_study,posy_study), 's', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 3);
        end
        drawnow;
        end
        
        t=t+dt; % Time update
    end % End temporal iteration

    if i == 1 % Only plot the graph for the first mesh size
        figure (1)
        legend('u_{analytic}', 'v_{analytic}','u_{numeric}' , 'v_{numeric}');
        hold off
        figure (5)
        legend('px_{analytic}','px_{numeric}' )
        hold off;  
    end


end




%% Validation


%% Plots
