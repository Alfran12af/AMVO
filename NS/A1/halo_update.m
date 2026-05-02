function F = halo_update(F)

% halo_update function updates the halo of a field introduced
% input:
%   F = field to update
% output:
%   F = field with the halo updated

N = size(F,1)-2; % Defining the Mesh size.

F(1,:)=F(N+1,:); % Column 1 becomes equal to the column N+1
F(N+2,:)=F(2,:); % Column N+2 becomes equal to the column 2

F(:,1)=F(:,N+1); % Row 1 becomes equal to the row N+1
F(:,N+2)=F(:,2); % Row N+2 becomes equal to the row 2
