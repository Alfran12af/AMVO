function [q_field] = vector2field(q_vec, N)
    % Transforms a vector into a field distribution
    % Input:
    %   q_vec: vector
    %   N: mesh size
    % Output: 
    %   q_field: field distribution
    
    
    % Preallocating variables
    q_field = zeros(N+2);
    
    % Calculation of the field distribution
    for i = 2:N+1
        for j = 2:N+1
             % Current control volume number of the grid 
             ncv = ((i-1) - 1)*N + j - 1;
             % Assigning the value of ncv to the field
             q_field(i,j) = q_vec(ncv);
        end
    end

    % Halo update
    q_field = halo_update(q_field);
end
