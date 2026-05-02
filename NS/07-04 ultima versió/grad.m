function [out_x, out_y] = grad(q, h, N)
%
% Computes the numerical gradient of a scalar field
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Inputs:
%   q : scalar field (with halo)
%   h : grid spacing
%   N : number of cells
%
% Outputs:
%   out_x : gradient of the field in x-direction
%   out_y : gradient of the field in y-direction
%

    % Allocate memory
    out_x = zeros(N+2, N+2);
    out_y = zeros(N+2, N+2);

    % --- Compute gradients on the mesh ---
    for j = 2:1:N+1
        for i = 2:1:N+1

            % Gradient in x-direction
            out_x(i,j) = (q(i+1,j) - q(i,j)) / h;

            % Gradient in y-direction
            out_y(i,j) = (q(i,j+1) - q(i,j)) / h;

        end
    end

    % Apply halo update
    out_x = halo_update(out_x);
    out_y = halo_update(out_y);

end