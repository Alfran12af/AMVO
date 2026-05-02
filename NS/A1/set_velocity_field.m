function [u_field, v_field] = set_velocity_field(N, u, v, ux, uy, vx, vy)

%
%
%
% 
syms x y

u_fun = matlabFunction(u, 'Vars', [x y]);
v_fun = matlabFunction(v, 'Vars', [x y]);

u_field = u_fun(ux, uy);
v_field = v_fun(vx, vy);

u_field = halo_update(u_field);
v_field = halo_update(v_field);


end