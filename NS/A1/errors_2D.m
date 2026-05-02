function [error_cu, error_cv, error_du, error_dv] = errors_2D(cu, cv, du, dv, cua, cva, dua, dva)

% Calculation of the maximum errors
    error_cu = max(max(abs(cu - cua))) * 100;
    error_cv = max(max(abs(cv - cva))) * 100;
    error_du = max(max(abs(du - dua))) * 100;
    error_dv = max(max(abs(dv - dva))) * 100;


end