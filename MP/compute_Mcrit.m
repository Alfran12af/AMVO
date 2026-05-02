function Mcrit = compute_Mcrit(M_inf, Cp_KT, Cp_star)

% -------------------------------------------------------------------------
% COMPUTE CRITICAL MACH NUMBER
%
% Busca la intersecció Cp_KT = Cp_star per cada alpha
%
% Inputs:
%   M_inf   -> vector Mach
%   Cp_KT   -> (Mach x alpha)
%   Cp_star -> (Mach x 1)
%
% Output:
%   Mcrit   -> (alpha x 1)
% -------------------------------------------------------------------------

[ n_M, n_alpha ] = size(Cp_KT);

Mcrit = zeros(n_alpha,1);

for j = 1:n_alpha

    % diferència entre corbes
    diff = abs(Cp_KT(:,j) - Cp_star);

    % índex del mínim
    [~, idx] = min(diff);

    Mcrit(j) = M_inf(idx);

end

end