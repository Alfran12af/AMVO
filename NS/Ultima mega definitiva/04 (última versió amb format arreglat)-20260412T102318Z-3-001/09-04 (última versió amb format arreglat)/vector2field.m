function q_field = vector2field(q_vec, N)
% VECTOR2FIELD Converts a vector into a scalar field (with halo)
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This function reconstructs a scalar field from its vectorized form,
%   typically obtained after solving a linear system such as the Poisson
%   equation:
%
%       A p = b
%
%   The vector is reshaped into a 2D field corresponding to the interior
%   nodes, and halo nodes are added to enforce periodic boundary conditions.
%
% Input:
%   q_vec : vectorized field (size: N^2 x 1)
%   N     : number of control volumes per direction
%
% Output:
%   q_field : scalar field with halo (size: (N+2) x (N+2))
%
% Notes:
%   - Interior values are reshaped using MATLAB column-wise ordering
%   - Must be consistent with FIELD2VECTOR function
%   - Halo nodes are updated using periodic boundary conditions
%

    %% 1. Memory Allocation

    q_field = zeros(N+2, N+2);


    %% 2. Field Reconstruction (Interior Nodes)

    q_field(2:N+1, 2:N+1) = reshape(q_vec, [N, N]);


    %% 3. Halo Update (Periodic Boundary Conditions)

    q_field = halo_update(q_field);

end