function dt = timestep(N, L, u_field, v_field, visc)
% TIMESTEP Computes stable time step based on CFL and diffusion criteria
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This function computes a stable time step for the explicit time
%   integration of the Navier–Stokes equations.
%
%   The time step is constrained by:
%     - Convective (CFL) condition
%     - Diffusive stability condition
%
%   The final time step is selected as the minimum of both constraints,
%   multiplied by a safety factor.
%
% Input:
%   N       : number of control volumes per direction
%   L       : domain length
%   u_field : horizontal velocity field
%   v_field : vertical velocity field
%   visc    : kinematic viscosity
%
% Output:
%   dt : stable time step
%
% Notes:
%   - CFL condition: dt ≤ h / max(|u|, |v|)
%   - Diffusive limit: dt ≤ h² / (2ν)
%   - A safety factor is applied for robustness
%

    %% 1. Grid Spacing

    h = L / N;


    %% 2. Maximum Velocities

    u_max = max(abs(u_field(:)));
    v_max = max(abs(v_field(:)));


    %% 3. Convective Stability Limit (CFL)

    dtconv_x = h / u_max;
    dtconv_y = h / v_max;
    dtconv   = min(dtconv_x, dtconv_y);


    %% 4. Diffusive Stability Limit

    dtdiff = 0.5 * (h^2 / visc);


    %% 5. Final Time Step (Safety Factor)

    f  = 0.1;  % safety factor
    dt = f * min(dtconv, dtdiff);

end