function [u_p, v_p, Ru, Rv] = predictor_velocity( ...
    u, v, cu, cv, du, dv, h, N, Ru_prev, Rv_prev, dt, visc_nu)
%
% Computes the predictor velocity using a second-order Adams-Bashforth scheme
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Inputs:
%   u        : horizontal velocity field at time n
%   v        : vertical velocity field at time n
%   cu, cv   : convective terms
%   du, dv   : diffusive terms
%   h        : grid spacing
%   N        : number of cells
%   Ru_prev  : residual at previous time step (u component)
%   Rv_prev  : residual at previous time step (v component)
%   dt       : time step
%   visc_nu  : kinematic viscosity
%
% Outputs:
%   u_p : predictor horizontal velocity field
%   v_p : predictor vertical velocity field
%   Ru  : current residual (u component)
%   Rv  : current residual (v component)
%

    % Allocate memory
    u_p = zeros(N+2, N+2);
    v_p = zeros(N+2, N+2);
    Ru  = zeros(N+2, N+2);
    Rv  = zeros(N+2, N+2);

    % --- Compute residuals and predictor velocity ---
    for j = 2:1:N+1
        for i = 2:1:N+1

            % --- Residuals (Navier--Stokes RHS without pressure) ---
            Ru(i,j) = -cu(i,j) + visc_nu * du(i,j);
            Rv(i,j) = -cv(i,j) + visc_nu * dv(i,j);

            % --- Adams--Bashforth (2nd order) predictor ---
            if all(Ru_prev(:) == 0)
            % --- First time step (Euler) ---
            u_p(i,j) = u(i,j) + dt * Ru(i,j);
            v_p(i,j) = v(i,j) + dt * Rv(i,j);
            else
            % --- Adams-Bashforth (2nd order) ---
            u_p(i,j) = u(i,j) + dt * ( ...
                3/2*Ru(i,j) - 1/2*Ru_prev(i,j) );
        
            v_p(i,j) = v(i,j) + dt * ( ...
                3/2*Rv(i,j) - 1/2*Rv_prev(i,j) );
            end

        end
    end

    % Apply halo update
    u_p = halo_update(u_p);
    v_p = halo_update(v_p);

end