function [conv_u_an, conv_v_an, diff_u_an, diff_v_an] = analytic_terms(u_sym, v_sym, N, xu, yu, xv, yv)
%
% Computes analytically the convective and diffusive terms of a velocity field
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez 
%
% Inputs:
%   u_sym  : symbolic expression of horizontal velocity
%   v_sym  : symbolic expression of vertical velocity
%   N      : number of cells
%   xu,yu  : coordinates where u is defined
%   xv,yv  : coordinates where v is defined
%
% Outputs:
%   conv_u_an : analytical convective term for u
%   conv_v_an : analytical convective term for v
%   diff_u_an : analytical diffusive term for u
%   diff_v_an : analytical diffusive term for v

    % Allocate memory
    conv_u_an = zeros(N+2, N+2);
    conv_v_an = zeros(N+2, N+2);
    diff_u_an = zeros(N+2, N+2);
    diff_v_an = zeros(N+2, N+2);

    syms x y

    % --- Convective terms ---
    Cu_sym = diff(u_sym^2, x) + diff(u_sym*v_sym, y);
    Cv_sym = diff(v_sym*u_sym, x) + diff(v_sym^2, y);

    % --- Diffusive terms (Laplacian) ---
    Du_sym = diff(u_sym, x, 2) + diff(u_sym, y, 2);
    Dv_sym = diff(v_sym, x, 2) + diff(v_sym, y, 2);

    % Convert symbolic expressions to callable functions
    Cu_function = matlabFunction(Cu_sym, 'Vars', [x y]);
    Cv_function = matlabFunction(Cv_sym, 'Vars', [x y]);
    Du_function = matlabFunction(Du_sym, 'Vars', [x y]);
    Dv_function = matlabFunction(Dv_sym, 'Vars', [x y]);

    % Evaluate fields on the mesh
    for i = 2:N+1
        for j = 2:N+1

            conv_u_an(i,j) = Cu_function(xu(i,j), yu(i,j));
            conv_v_an(i,j) = Cv_function(xv(i,j), yv(i,j));

            diff_u_an(i,j) = Du_function(xu(i,j), yu(i,j));
            diff_v_an(i,j) = Dv_function(xv(i,j), yv(i,j));

        end
    end
end