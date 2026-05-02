function print_field(u)
%
% Prints a field following the notation used in the course slides
%
% Author: Manel Soria (2025)
%
% Inputs:
%   u : field to print
%
% Outputs:
%   (none)

    % Mesh size
    N = size(u,1) - 2;

    % --- Print field values ---
    for j = N+2:-1:1

        fprintf('j=%2d ', j);

        for i = 1:N+2
            fprintf(' %+8.3e', u(i,j));
        end

        fprintf('\n');

    end

end