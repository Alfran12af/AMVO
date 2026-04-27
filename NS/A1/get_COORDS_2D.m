function [ux, uy, vx, vy, px, py] = get_COORDS_2D(L, N)
% Gets the node distribution for a 2D symmetrical case
% INPUT
%
%
%

D = L/N;
ux = zeros(N+2, N+2);
uy = zeros(N+2, N+2);
vx = zeros(N+2, N+2);
vy = zeros(N+2, N+2);
px = zeros(N+2, N+2);
py = zeros(N+2, N+2);

for j=2:1:N+1

    for i=2:1:N+1
        ux(i,j) = (i-1)*D;
        vx(i,j) = (i-1)*D - D/2;
        px(i,j) = (i-1)*D - D/2;

        uy(i,j) = (j-1)*D - D/2;
        vy(i,j) = (j-1)*D;
        py(i,j) = (j-1)*D - D/2;
    end



end



