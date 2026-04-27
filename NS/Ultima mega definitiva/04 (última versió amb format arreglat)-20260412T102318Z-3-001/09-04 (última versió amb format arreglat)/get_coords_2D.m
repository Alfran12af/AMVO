function [ux, uy, vx, vy, px, py] = get_coords_2D(L, N)
% GET_COORDS_2D Computes node coordinates for a 2D staggered grid (Arakawa-C)
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This function generates the spatial coordinates of velocity and pressure
%   nodes for a 2D structured staggered grid (Arakawa-C configuration).
%
%   In this arrangement:
%     - u-velocity is defined at vertical cell faces
%     - v-velocity is defined at horizontal cell faces
%     - pressure is defined at cell centers
%
%   The grid is uniform and periodic, consistent with the formulation used
%   throughout the project.
%
% Input:
%   L : domain length
%   N : number of control volumes per direction
%
% Output:
%   ux : x-coordinates of u-velocity nodes
%   uy : y-coordinates of u-velocity nodes
%   vx : x-coordinates of v-velocity nodes
%   vy : y-coordinates of v-velocity nodes
%   px : x-coordinates of pressure nodes (cell centers)
%   py : y-coordinates of pressure nodes (cell centers)
%
% Notes:
%   - Grid includes halo nodes (size: (N+2) x (N+2))
%   - Coordinates correspond to interior nodes only (halo not populated)
%   - Uniform spacing: Δ = L / N
%

    %% 1. Grid Spacing

    D = L / N;


    %% 2. Memory Allocation

    ux = zeros(N+2, N+2);
    uy = zeros(N+2, N+2);
    vx = zeros(N+2, N+2);
    vy = zeros(N+2, N+2);
    px = zeros(N+2, N+2);
    py = zeros(N+2, N+2);


    %% 3. Coordinate Evaluation (Interior Nodes)

    for j = 2:N+1
        for i = 2:N+1

            % --- 3.1 u-velocity nodes (vertical faces) ---
            ux(i,j) = (i-1) * D;
            uy(i,j) = (j-1) * D - D/2;

            % --- 3.2 v-velocity nodes (horizontal faces) ---
            vx(i,j) = (i-1) * D - D/2;
            vy(i,j) = (j-1) * D;

            % --- 3.3 pressure nodes (cell centers) ---
            px(i,j) = (i-1) * D - D/2;
            py(i,j) = (j-1) * D - D/2;

        end
    end

end