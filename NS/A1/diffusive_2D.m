function [dif_u, dif_v] = diffusive_2D(N, u, v, h)
%
%
% Dado que es una malla uniforme donde Delta_x y Delta_y son iguales y de
% valor Delta, se anulan estos valores junto con el denominador de la
% derivada
dif_u = zeros(N+2, N+2);
dif_v = zeros(N+2, N+2);

for j = 2:1:N+1
    for i = 2:1:N+1
        % East
        dudx_e = u(i+1,j) - u(i,j);
        dvdx_e = v(i+1,j) - v(i,j);

        % West
        dudx_w = u(i,j) - u(i-1,j);
        dvdx_w = v(i,j) - v(i-1,j);

        % North
        dudy_n = u(i,j+1) - u(i,j);
        dvdy_n = v(i,j+1) - v(i,j);


        % South
        dudy_s = u(i,j) - u(i,j-1);
        dvdy_s = v(i,j) - v(i,j-1);

        % Sum equals int_V of \nabla · grad(velocity) dV 
        dif_u(i,j) = dudx_e + dudy_n - dudx_w - dudy_s;
        dif_v(i,j) = dvdx_e + dvdy_n - dvdx_w - dvdy_s;


    end
end

dif_u = dif_u./h^2;
dif_v = dif_v./h^2;

end

