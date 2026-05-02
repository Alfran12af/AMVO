function q_field = vector2field(q_vec, N)
%
% Converts a vector into a scalar field (with halo)
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Inputs:
%   q_vec : vectorized field (size: N^2 x 1)
%   N     : number of cells
%
% Outputs:
%   q_field : scalar field with halo (size: (N+2)x(N+2))
%

    % --- Allocate memory ---
    q_field = zeros(N+2, N+2);

    % --- Reconstruct internal field from vector ---
    q_field(2:N+1, 2:N+1) = reshape(q_vec, [N, N]);

    % --- Apply halo update ---
    q_field = halo_update(q_field);

end