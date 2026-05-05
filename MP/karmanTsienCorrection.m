function [Cp_KT, Cp_star] = karmanTsienCorrection(N, alpha_deg, M_inf, gamma, cp_values)

% -------------------------------------------------------------------------
% KARMAN-TSIEN COMPRESSIBILITY CORRECTION
%
% Computes:
%   - Compressibility-corrected minimum pressure coefficient (Cp)
%   - Critical (sonic) pressure coefficient Cp*
%
% Inputs:
%   N           -> number of panels (not explicitly used, kept for consistency)
%   alpha_deg   -> angle of attack (degrees)
%   M_inf       -> Mach number vector
%   gamma       -> ratio of specific heats
%   cp_values   -> incompressible Cp (alpha x panels)
%
% Outputs:
%   Cp_KT       -> corrected Cp (Mach x alpha)
%   Cp_star     -> critical Cp* (Mach x 1)
% -------------------------------------------------------------------------

n_alpha = length(alpha_deg);
n_M = length(M_inf);

Cp_KT   = zeros(n_M, n_alpha);
Cp_star = zeros(n_M,1);

%% ------------------------------------------------------------------------
% Sonic Cp*

for i = 1:n_M

    M = M_inf(i);

    Cp_star(i) = (2/(gamma*M^2)) * ...
        ( ((1 + (gamma-1)/2*M^2)/(1 + (gamma-1)/2))^(gamma/(gamma-1)) - 1 );

end

%% ------------------------------------------------------------------------
% Karman-Tsien correction (using minimum Cp)

for j = 1:n_alpha

    Cp0 = min(cp_values(j,:)); % critical Cp at each alpha

    for i = 1:n_M

        M = M_inf(i);
        beta = sqrt(1 - M^2);

        Cp_KT(i,j) = Cp0 / (beta + (M^2/(1+beta))*(Cp0/2));

    end

end

end