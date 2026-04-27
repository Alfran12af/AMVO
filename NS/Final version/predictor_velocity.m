function [u_p, v_p, Ru, Rv, dt] = predictor_velocity(u, v, cu, cv, du, dv, h, Ruanterior, Rvanterior, dt, visc_nu)

    % Prellocate variables
    u_p = zeros(size(u));
    v_p = zeros(size(v));
    Ru = zeros(size(u));
    Rv = zeros(size(cv));
    
    for i = 2:N+1
        for j = 2:N+1
            % Known values converted to a constant R
            Ru(i,j) = -(cu(i,j)./h^2) + (du(i,j)./h^2)*visc_nu; 
            Rv(i,j) = -(cv(i,j)./h^2) + (dv(i,j)./h^2)*visc_nu; 
        
            %Velocity predictor
            u_p(i,j) = u(i,j) + dt*(3/2*Ru(i,j) - 1/2*Ruanterior(i,j));
            v_p(i,j) = v(i,j) + dt*(3/2*Rv(i,j) - 1/2*Rvanterior(i,j));
        end
    end
    % Halo Update
    u_p = halo_update(u_p);
    v_p = halo_update(v_p);
end