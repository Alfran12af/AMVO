function print_field(u)
% PRINT_FIELD Prints a 2D field following course notation (with halo)
%
% Author: Manel Soria (2025)
%
% Description:
%   This function prints a 2D scalar or vector field to the console using
%   the indexing convention adopted in the course slides.
%
%   The field is displayed row by row, starting from the top boundary
%   (highest j index) down to the bottom, matching the physical layout
%   of the domain.
%
%   Values are printed using scientific notation for clarity, which is
%   particularly useful in CFD applications where magnitudes can vary
%   significantly.
%
% Input:
%   u : field to print (including halo nodes)
%
% Output:
%   (none)
%
% Notes:
%   - Field size is assumed to be (N+2) x (N+2), including halo nodes
%   - Output follows (i,j) indexing with j decreasing vertically
%   - Useful for debugging and verification of field values
%

    %% 1. Mesh Size

    N = size(u,1) - 2;


    %% 2. Field Printing

    for j = N+2:-1:1

        fprintf('j=%2d ', j);

        for i = 1:N+2
            fprintf(' %+8.3e', u(i,j));
        end

        fprintf('\n');

    end

end