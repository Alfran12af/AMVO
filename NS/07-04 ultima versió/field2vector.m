function b = field2vector(div, N)
%
% Converts a scalar field (with halo) into a vector form
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Inputs:
%   div : scalar field with halo (size: (N+2)x(N+2))
%   N   : number of cells
%
% Outputs:
%   b   : vectorized field without halo (size: N^2 x 1)
%

    % --- Extract internal field (remove halo) ---
    field_internal = div(2:N+1, 2:N+1);

    % --- Convert field into column vector ---
    b = reshape(field_internal, [N*N, 1]);

end