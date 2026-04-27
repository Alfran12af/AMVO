function [conv_u_an, conv_v_an, diff_u_an, diff_v_an] = analytic_terms( ...
    u_sym, v_sym, N, xu, yu, xv, yv)
% ANALYTIC_TERMS Computes analytical convective and diffusive terms (MMS)
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This function computes the analytical convective and diffusive terms
%   of a 2D velocity field using symbolic differentiation.
%
%   The convective terms are evaluated in conservative (divergence) form,
%   consistent with the finite control-volume (FCV) discretization used
%   in Part A:
%
%       ∇·(u_i u) = ∂(u_i u)/∂x + ∂(u_i v)/∂y
%
%   The diffusive terms correspond to the Laplacian of each velocity
%   component:
%
%       ∇²u = ∂²u/∂x² + ∂²u/∂y²
%
%   These analytical expressions are used for verification through the
%   Method of Manufactured Solutions (MMS).
%
% Input:
%   u_sym : symbolic expression of horizontal velocity
%   v_sym : symbolic expression of vertical velocity
%   N     : number of control volumes per direction
%   xu,yu : coordinates of u-velocity nodes (staggered grid)
%   xv,yv : coordinates of v-velocity nodes
%
% Output:
%   conv_u_an : analytical convective term for u-component
%   conv_v_an : analytical convective term for v-component
%   diff_u_an : analytical diffusive term for u-component
%   diff_v_an : analytical diffusive term for v-component
%
% Notes:
%   - Convective term is computed in divergence form (important for MMS)
%   - Fields are evaluated on staggered grid locations
%   - Periodic boundary conditions are enforced via halo update
%

    %% 1. Memory Allocation

    conv_u_an = zeros(N+2, N+2);
    conv_v_an = zeros(N+2, N+2);
    diff_u_an = zeros(N+2, N+2);
    diff_v_an = zeros(N+2, N+2);


    %% 2. Symbolic Definitions

    syms x y


    %% 3. Analytical Expressions

    % --- 3.1 Convective terms (conservative form) ---
    Cu_sym = diff(u_sym * u_sym, x) + diff(u_sym * v_sym, y);
    Cv_sym = diff(v_sym * u_sym, x) + diff(v_sym * v_sym, y);

    % --- 3.2 Diffusive terms (Laplacian) ---
    Du_sym = diff(u_sym, x, 2) + diff(u_sym, y, 2);
    Dv_sym = diff(v_sym, x, 2) + diff(v_sym, y, 2);


    %% 4. Symbolic → Numerical Functions

    Cu_function = matlabFunction(Cu_sym, 'Vars', [x y]);
    Cv_function = matlabFunction(Cv_sym, 'Vars', [x y]);
    Du_function = matlabFunction(Du_sym, 'Vars', [x y]);
    Dv_function = matlabFunction(Dv_sym, 'Vars', [x y]);


    %% 5. Field Evaluation (Interior Nodes)

    for i = 2:N+1
        for j = 2:N+1

            % --- 5.1 Convective terms ---
            conv_u_an(i,j) = Cu_function(xu(i,j), yu(i,j));
            conv_v_an(i,j) = Cv_function(xv(i,j), yv(i,j));

            % --- 5.2 Diffusive terms ---
            diff_u_an(i,j) = Du_function(xu(i,j), yu(i,j));
            diff_v_an(i,j) = Dv_function(xv(i,j), yv(i,j));

        end
    end


    %% 6. Halo Update (Periodic Boundary Conditions)

    conv_u_an = halo_update(conv_u_an);
    conv_v_an = halo_update(conv_v_an);
    diff_u_an = halo_update(diff_u_an);
    diff_v_an = halo_update(diff_v_an);

end