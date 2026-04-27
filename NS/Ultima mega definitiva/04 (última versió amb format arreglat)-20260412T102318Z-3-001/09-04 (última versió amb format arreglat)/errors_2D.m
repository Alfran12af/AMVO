function [error_cu, error_cv, error_du, error_dv] = errors_2D( ...
    N, cu, cv, du, dv, cua, cva, dua, dva)
% ERRORS_2D Computes maximum errors between numerical and analytical terms
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Description:
%   This function computes the maximum absolute error between numerical
%   and analytical convective and diffusive terms, as part of the
%   Method of Manufactured Solutions (MMS) verification.
%
%   The error is evaluated using the infinity norm (maximum norm) over
%   the interior nodes of the domain.
%
% Input:
%   N   : number of control volumes per direction
%   cu  : numerical convective term for u-component
%   cv  : numerical convective term for v-component
%   du  : numerical diffusive term for u-component
%   dv  : numerical diffusive term for v-component
%   cua : analytical convective term for u-component
%   cva : analytical convective term for v-component
%   dua : analytical diffusive term for u-component
%   dva : analytical diffusive term for v-component
%
% Output:
%   error_cu : maximum error (∞-norm) of convective term for u
%   error_cv : maximum error (∞-norm) of convective term for v
%   error_du : maximum error (∞-norm) of diffusive term for u
%   error_dv : maximum error (∞-norm) of diffusive term for v
%
% Notes:
%   - Halo nodes are excluded from error computation
%   - Used for mesh convergence analysis (expected O(h^2))
%

    %% 1. Error Computation (Interior Nodes)

    error_cu = max(max(abs(cu(2:N+1,2:N+1) - cua(2:N+1,2:N+1))));
    error_cv = max(max(abs(cv(2:N+1,2:N+1) - cva(2:N+1,2:N+1))));
    error_du = max(max(abs(du(2:N+1,2:N+1) - dua(2:N+1,2:N+1))));
    error_dv = max(max(abs(dv(2:N+1,2:N+1) - dva(2:N+1,2:N+1))));

end