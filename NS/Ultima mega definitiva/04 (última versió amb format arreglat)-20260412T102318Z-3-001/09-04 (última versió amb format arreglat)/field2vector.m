function b = field2vector(div, N)
% FIELD2VECTOR Converts a scalar field (with halo) into vector form
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This function converts a scalar field defined on a structured grid
%   (including halo nodes) into a column vector containing only the
%   interior values.
%
%   The resulting vector is used as the right-hand side of the linear
%   system associated with the Poisson equation:
%
%       A p = b
%
%   where b corresponds to the divergence of the predictor velocity.
%
% Input:
%   div : scalar field with halo (size: (N+2) x (N+2))
%   N   : number of control volumes per direction
%
% Output:
%   b   : column vector of size (N^2 x 1), containing interior values
%
% Notes:
%   - Halo nodes are excluded from the vectorization
%   - Ordering follows MATLAB column-wise convention (reshape)
%

    %% 1. Extract Interior Field (Remove Halo)

    field_internal = div(2:N+1, 2:N+1);


    %% 2. Vectorization (Column-Wise)

    b = reshape(field_internal, [N*N, 1]);

end