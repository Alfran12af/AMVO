function [u_an, v_an, px_an, py_an] = analytic_temporal(ux, uy, vx, vy, px, py, time, visc_nu, u_sym, v_sym, rho, N, h)
    % Obtains the analytic solutions of the problem
    % Input:
    %   ux: x coordinates for horizontal velocity
    %   uy: y coordinates for horizontal velocity
    %   vx: x coordinates for vertical velocity
    %   vy: y coordinates for vertical velocity
    %   px: x coordinates for pressure
    %   py: y coordinates for pressure
    %   time: value of the time of the iteration
    %   visc_nu: viscosity of the problem
    %   u_sym: equation for the horizontal velocity
    %   v_sym: equation for the vertical velocity
    %   N: mesh size
    %   h: space between nodes
    % Output: 
    %   u_an: analytic result of the velocity for direction-x
    %   v_an: analytic result of the velocity for direction-y
    %   px_an: analytic result of the pressure for direction-x
    %   py_an: analytic result of the pressure for direction-y


    % Preallocating values
    u_an = zeros(N+2);
    v_an = zeros(N+2);
    p_field_an = zeros(N+2);
    
    % syms x y
    
    fu = matlabFunction(u_sym, 'Vars', [x y]);
    fv = matlabFunction(v_sym, 'Vars', [x y]);

    % F factor
    F = exp(-8*pi^2*visc_nu*time);
    
    % Computation of velocity and pressure fields at time t
    for i = 1:N+2
        for j = 1:N+2
            u_an(i,j) = F*fu(ux(i,j), uy(i,j));
            v_an(i,j) = F*fv(vx(i,j), vy(i,j));
            p_field_an(i,j) = -rho*F^2*(cos(2*pi*px(i,j))^2 + cos(2*pi*py(i,j))^2)/2;
        end
    end
    
    [px_an, py_an] = grad(p_field_an, h, N);
end