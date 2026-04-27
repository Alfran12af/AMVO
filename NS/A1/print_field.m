function print_field( u )
% print_field( u ) prints a field with the notation described in
% slides
% u: field to print
% Written by: Manel Soria 2025
% Example of use: print_field(up);
N=size(u,1)-2; % mesh size
for j=N+2:-1:1
fprintf('j=%2d ',j);
for i=1:N+2
fprintf(' %+8.3e' ,u(i,j));
end
fprintf('\n');
end
end