function div = diverg(u, v, h, N)
%
% Computes the numerical divergence of a velocity field
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Inputs:
%   u : horizontal velocity field
%   v : vertical velocity field
%   h : grid spacing
%   N : number of cells
%
% Outputs:
%   div : numerical divergence field (integral form)
%

    % Allocate memory
    div = zeros(N+2, N+2);

    % Evaluate divergence on the mesh
    for j = 2:1:N+1
        for i = 2:1:N+1

            % --- Velocity components ---
            u_p = u(i,j);
            v_p = v(i,j);

            u_w = u(i-1,j);
            v_s = v(i,j-1);

            % --- Divergence (finite volume form) ---
            div(i,j) = h * (u_p - u_w + v_p - v_s);

        end
    end

    % Apply halo update
    div = halo_update(div);

end

