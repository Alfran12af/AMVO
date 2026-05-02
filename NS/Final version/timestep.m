function [dt] = timestep(N, L, u_field, v_field, visc)
    
    h = L/N;
    
    % Maximum velocities:
    u_max = max(max(abs(u_field(:))));
    v_max = max(max(abs(v_field(:))));
    
    % Minimum convective time limit:
    dtconv = min((h)/u_max, (h)/v_max);
    
    % Diffusive time limit:
    dtdiff = 0.5*(h)^2/visc;
    
    % Minimum time step possible
    f = 0.1; % f = security factor
    dt = f*min(dtconv, dtdiff);

end