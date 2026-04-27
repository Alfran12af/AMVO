function [cu, cv] = convective_2D(u, v, h, N)
%
% Computes the numerical convective terms of a velocity field
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Inputs:
%   u  : horizontal velocity field
%   v  : vertical velocity field
%   h  : grid spacing
%   N  : number of cells
%
% Outputs:
%   cu : numerical convective term for u
%   cv : numerical convective term for v

    % Allocate memory
    cu = zeros(N+2, N+2);
    cv = zeros(N+2, N+2);

    % Evaluate convective terms on the mesh
    for i = 2:N+1
        for j = 2:N+1

            % --- Calculation of the horizontal and vertical face velocities ---
            ue = (u(i+1,j) + u(i,j)) / 2;
            uw = (u(i-1,j) + u(i,j)) / 2;
            un = (u(i,j+1) + u(i,j)) / 2;
            us = (u(i,j-1) + u(i,j)) / 2;

            ve = (v(i+1,j) + v(i,j)) / 2;
            vw = (v(i-1,j) + v(i,j)) / 2;
            vn = (v(i,j+1) + v(i,j)) / 2;
            vs = (v(i,j-1) + v(i,j)) / 2;

            % --- Calculation of the horizontal and vertical face flow terms ---
            Feh = h * (u(i+1,j) + u(i,j)) / 2;
            Fwh = h * (u(i-1,j) + u(i,j)) / 2;
            Fnh = h * (v(i,j) + v(i+1,j)) / 2;
            Fsh = h * (v(i,j-1) + v(i+1,j-1)) / 2;

            Fev = h * (u(i,j+1) + u(i,j)) / 2;
            Fwv = h * (u(i-1,j+1) + u(i-1,j)) / 2;
            Fnv = h * (v(i,j+1) + v(i,j)) / 2;
            Fsv = h * (v(i,j) + v(i,j-1)) / 2;

            % --- Calculation of the convective terms ---
            cu(i,j) = ue*Feh - uw*Fwh + un*Fnh - us*Fsh;
            cv(i,j) = ve*Fev - vw*Fwv + vn*Fnv - vs*Fsv;

        end
    end

    % Normalize by the control volume area
    cu = cu./(h)^2;
    cv = cv./(h)^2;

end