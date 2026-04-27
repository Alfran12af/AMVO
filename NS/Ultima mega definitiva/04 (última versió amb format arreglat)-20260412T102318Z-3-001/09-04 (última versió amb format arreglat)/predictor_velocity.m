function [u_p, v_p, Ru, Rv] = predictor_velocity( ...
    u, v, cu, cv, du, dv, h, N, Ru_prev, Rv_prev, dt, visc_nu)
% PREDICTOR_VELOCITY Computes predictor velocity using Adams–Bashforth (2nd order)
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This function computes the predictor (intermediate) velocity field
%   using an explicit second-order Adams–Bashforth time integration scheme.
%
%   The predictor step corresponds to:
%
%       du/dt = R(u) = -C(u) + νD(u)
%
%   where:
%     - C(u) is the convective term
%     - D(u) is the diffusive term
%
%   The pressure gradient is NOT included in this step. The resulting
%   velocity field is later corrected using the projection method.
%
%   For the first time step, a first-order Euler scheme is used.
%
% Input:
%   u, v     : velocity fields at time step n
%   cu, cv   : convective terms
%   du, dv   : diffusive terms
%   h        : grid spacing (uniform mesh)
%   N        : number of control volumes per direction
%   Ru_prev  : residual at previous time step (u-component)
%   Rv_prev  : residual at previous time step (v-component)
%   dt       : time step
%   visc_nu  : kinematic viscosity
%
% Output:
%   u_p : predictor velocity (u-component)
%   v_p : predictor velocity (v-component)
%   Ru  : current residual (u-component)
%   Rv  : current residual (v-component)
%
% Notes:
%   - Explicit time integration (Adams–Bashforth, 2nd order)
%   - First step uses Euler scheme
%   - Output is not divergence-free (requires projection step)
%

    %% 1. Memory Allocation

    u_p = zeros(N+2, N+2);
    v_p = zeros(N+2, N+2);
    Ru  = zeros(N+2, N+2);
    Rv  = zeros(N+2, N+2);


    %% 2. Residual and Predictor Computation (Interior Nodes)

    for j = 2:N+1
        for i = 2:N+1

            %% 2.1 Residuals (Navier–Stokes RHS without pressure)

            Ru(i,j) = -cu(i,j) + visc_nu * du(i,j);
            Rv(i,j) = -cv(i,j) + visc_nu * dv(i,j);


            %% 2.2 Time Integration (Adams–Bashforth)

            if all(Ru_prev(:) == 0)

                % --- First time step (Euler scheme) ---
                u_p(i,j) = u(i,j) + dt * Ru(i,j);
                v_p(i,j) = v(i,j) + dt * Rv(i,j);

            else

                % --- Adams–Bashforth (2nd order) ---
                u_p(i,j) = u(i,j) + dt * ( ...
                    3/2 * Ru(i,j) - 1/2 * Ru_prev(i,j) );

                v_p(i,j) = v(i,j) + dt * ( ...
                    3/2 * Rv(i,j) - 1/2 * Rv_prev(i,j) );

            end

        end
    end


    %% 3. Halo Update (Periodic Boundary Conditions)

    u_p = halo_update(u_p);
    v_p = halo_update(v_p);

end