function [A] = laplacian_matrix(N)
%
% Builds the discrete Laplacian matrix for a 2D NxN grid
% using periodic boundary conditions (halo-compatible)
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Inputs:
%   N : number of cells per direction
%
% Outputs:
%   A : Laplacian matrix (N^2 x N^2)

    % Initialize matrix
    A = zeros(N*N);

    % Loop over grid cells
    for row = 1:N
        for col = 1:N

            % Linear index of current cell
            idx = (row-1)*N + col;

            % Central coefficient
            A(idx, idx) = -4;

            % Neighbor indices

            % North
            if row < N
                idx_north = row*N + col;
            else
                idx_north = col;
            end

            % South
            if row > 1
                idx_south = (row-2)*N + col;
            else
                idx_south = (N-1)*N + col;
            end

            % East
            if col < N
                idx_east = (row-1)*N + (col+1);
            else
                idx_east = (row-1)*N + 1;
            end

            % West
            if col > 1
                idx_west = (row-1)*N + (col-1);
            else
                idx_west = (row-1)*N + N;
            end

            % Assign neighbor contributions
            A(idx, idx_north) = 1;
            A(idx, idx_south) = 1;
            A(idx, idx_east)  = 1;
            A(idx, idx_west)  = 1;

        end
    end

    % --- Fix singularity ---
    A(1,1) = -5;

end