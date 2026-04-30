function [Cp_KT, Cp_star] = KarmanTsien(N, alpha, M_inf, gamma, cp_values)

% -------------------------------------------------------------------------
% KARMAN-TSIEN COMPRESSIBILITY CORRECTION
%
% cp_values -> (alpha x panels)
% Cp_KT     -> (Mach x alpha)
% Cp_star   -> (Mach x 1)
% -------------------------------------------------------------------------

n_alpha = length(alpha);
n_M = length(M_inf);

Cp_KT = zeros(n_M, n_alpha);
Cp_star = zeros(n_M,1);

%% ------------------------------------------------------------------------
% Cp* (condició sònica)

for i = 1:n_M

    M = M_inf(i);

    Cp_star(i) = (2/(gamma*M^2)) * ...
        ( ((1 + (gamma-1)/2*M^2)/(1 + (gamma-1)/2))^(gamma/(gamma-1)) - 1 );

end

%% ------------------------------------------------------------------------
% Cp amb Karman-Tsien

for j = 1:n_alpha

    % 🔴 TEORIA: el punt crític és el Cp mínim
    Cp0 = min(cp_values(j,:));

    for i = 1:n_M

        M = M_inf(i);

        beta = sqrt(1 - M^2);

        Cp_KT(i,j) = (Cp0 / beta) / ...
            (1 + (gamma-1)/2 * M^2 * (Cp0/(2*beta)));

    end
end

end