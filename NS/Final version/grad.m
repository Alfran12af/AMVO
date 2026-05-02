function [out_x, out_y] = grad(q, h, N)
    % Calculates the gradients from a field distribution, specifically
    % applied to pressure gradients
    % Input:
    %   q: field distribution
    %   N: mesh size
    %   h: space between nodes
    % Output: 
    %   [out_x]: gradient-x of the field distribution
    %   [out_y]: gradient-y of the field distribution
    

    % Preallocating variables
    out_x = zeros(N+2, N+2);
    out_y = zeros(N+2, N+2);

    % Computing of the gradients
    for i=2:N+1
        for j=2:N+1
            out_x(i,j) = (q(i+1,j) - q(i,j))/h;
            out_y(i,j) = (q(i,j+1) - q(i,j))/h;
        end
    end

    % Halo update
    out_x = halo_update(out_x);
    out_y = halo_update(out_y);
end
