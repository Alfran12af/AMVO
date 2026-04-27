function [error_cu, error_cv, error_du, error_dv] = errors_2D(N, cu, cv, du, dv, cua, cva, dua, dva)
%
% Computes the maximum errors between numerical and analytical terms
%
% Authors: Joan Duro | Axl Francia | Pol Jimenez
%
% Inputs:
%   cu  : numerical convective term for u
%   cv  : numerical convective term for v
%   du  : numerical diffusive term for u
%   dv  : numerical diffusive term for v
%   cua : analytical convective term for u
%   cva : analytical convective term for v
%   dua : analytical diffusive term for u
%   dva : analytical diffusive term for v
%
% Outputs:
%   error_cu : maximum error of convective term for u
%   error_cv : maximum error of convective term for v
%   error_du : maximum error of diffusive term for u
%   error_dv : maximum error of diffusive term for v

    % --- Calculation of the maximum errors ---
    error_cu = max(max(abs(cu(2:N+1,2:N+1) - cua(2:N+1,2:N+1))));
    error_cv = max(max(abs(cv(2:N+1,2:N+1) - cva(2:N+1,2:N+1))));
    error_du = max(max(abs(du(2:N+1,2:N+1) - dua(2:N+1,2:N+1))));
    error_dv = max(max(abs(dv(2:N+1,2:N+1) - dva(2:N+1,2:N+1))));

end