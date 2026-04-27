function A = laplacian_matrix(N)
% LAPLACIAN_MATRIX Builds the discrete Laplacian matrix (2D, periodic BCs)
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This function constructs the discrete Laplacian operator for a 2D
%   structured grid using a finite difference / finite volume formulation.
%
%   The matrix corresponds to the Poisson equation:
%
%       ∇²p = b   →   A p = b
%
%   where A is a sparse matrix of size (N^2 x N^2), including periodic
%   boundary conditions in both spatial directions.
%
%   Each row represents the 5-point stencil:
%
%       -4p(i,j) + p(i+1,j) + p(i-1,j) + p(i,j+1) + p(i,j-1)
%
% Input:
%   N : number of control volumes per direction
%
% Output:
%   A : Laplacian matrix of size (N^2 x N^2)
%
% Notes:
%   - Periodic boundary conditions are enforced via index wrapping
%   - Matrix is singular (null space = constant pressure)
%   - Singularity is removed by modifying one diagonal entry
%

    %% 1. Matrix Initialization

    A = zeros(N*N);


    %% 2. Assembly of Laplacian Operator

    for row = 1:N
        for col = 1:N

            %% 2.1 Linear Index

            idx = (row-1)*N + col;


            %% 2.2 Central Coefficient

            A(idx, idx) = -4;


            %% 2.3 Neighbor Indices (Periodic Wrapping)

            % North neighbor
            if row < N
                idx_north = row*N + col;
            else
                idx_north = col;
            end

            % South neighbor
            if row > 1
                idx_south = (row-2)*N + col;
            else
                idx_south = (N-1)*N + col;
            end

            % East neighbor
            if col < N
                idx_east = (row-1)*N + (col+1);
            else
                idx_east = (row-1)*N + 1;
            end

            % West neighbor
            if col > 1
                idx_west = (row-1)*N + (col-1);
            else
                idx_west = (row-1)*N + N;
            end


            %% 2.4 Assign Neighbor Contributions

            A(idx, idx_north) = 1;
            A(idx, idx_south) = 1;
            A(idx, idx_east)  = 1;
            A(idx, idx_west)  = 1;

        end
    end


    %% 3. Singularity Fix

    % Remove null space (constant pressure mode)
    A(1,1) = -5;

end