function F = halo_update(F)
%
% Updates the halo of a field using periodic boundaries
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Inputs:
%   F : field to update
%
% Outputs:
%   F : field with the halo updated

    % Mesh size
    N = size(F,1) - 2;

    % --- Horizontal halo update ---
    F(1,:)   = F(N+1,:); % Row 1 becomes equal to row N+1
    F(N+2,:) = F(2,:);   % Row N+2 becomes equal to row 2

    % --- Vertical halo update ---
    F(:,1)   = F(:,N+1); % Column 1 becomes equal to column N+1
    F(:,N+2) = F(:,2);   % Column N+2 becomes equal to column 2

end