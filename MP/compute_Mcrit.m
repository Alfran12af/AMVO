function Mcrit = compute_Mcrit(M_inf, Cp_KT, Cp_star)

% -------------------------------------------------------------------------
% COMPUTE CRITICAL MACH NUMBER
%
% Computes the critical Mach number (Mcrit) by finding the intersection
% between the Karman-Tsien corrected pressure coefficient (Cp_KT)
% and the critical pressure coefficient (Cp*).
%
% The intersection is obtained using linear interpolation between points
% where a sign change occurs.
%
% Inputs:
%   M_inf   -> vector of Mach numbers
%   Cp_KT   -> corrected Cp (Mach x alpha)
%   Cp_star -> critical Cp* (Mach x 1)
%
% Output:
%   Mcrit   -> critical Mach number for each angle of attack (alpha x 1)
% -------------------------------------------------------------------------

[ n_M, n_alpha ] = size(Cp_KT);

Mcrit = zeros(n_alpha,1);

for j = 1:n_alpha

    % Difference between both curves
    diff = Cp_KT(:,j) - Cp_star;

    % --- Find sign change (intersection interval)
    idx = find(diff(1:end-1).*diff(2:end) < 0, 1);

    if ~isempty(idx)
        
        % --- Linear interpolation
        M1 = M_inf(idx);
        M2 = M_inf(idx+1);
        
        f1 = diff(idx);
        f2 = diff(idx+1);
        
        Mcrit(j) = M1 - f1 * (M2 - M1) / (f2 - f1);
        
    else
        % --- Fallback: closest point if no sign change is detected
        [~, idx_min] = min(abs(diff));
        Mcrit(j) = M_inf(idx_min);
    end

end

end