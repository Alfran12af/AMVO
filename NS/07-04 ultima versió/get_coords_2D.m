function [ux, uy, vx, vy, px, py] = get_coords_2D(L, N)
%
% Computes the node coordinates for a 2D staggered grid
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Inputs:
%   L : domain length
%   N : number of cells
%
% Outputs:
%   ux : x-coordinates of u velocity nodes
%   uy : y-coordinates of u velocity nodes
%   vx : x-coordinates of v velocity nodes
%   vy : y-coordinates of v velocity nodes
%   px : x-coordinates of pressure nodes
%   py : y-coordinates of pressure nodes

    % Grid spacing
    D = L/N;

    % Allocate memory
    ux = zeros(N+2, N+2);
    uy = zeros(N+2, N+2);
    vx = zeros(N+2, N+2);
    vy = zeros(N+2, N+2);
    px = zeros(N+2, N+2);
    py = zeros(N+2, N+2);

    % Evaluate coordinates on the mesh
    for j = 2:1:N+1
        for i = 2:1:N+1

            % --- u velocity nodes ---
            ux(i,j) = (i-1)*D;
            uy(i,j) = (j-1)*D - D/2;

            % --- v velocity nodes ---
            vx(i,j) = (i-1)*D - D/2;
            vy(i,j) = (j-1)*D;

            % --- pressure nodes ---
            px(i,j) = (i-1)*D - D/2;
            py(i,j) = (j-1)*D - D/2;

        end
    end

end


