function Cl_Mach = compute_Cl_KT(c, Mach_list, N, nv, pl, alpha_deg, cp)

% -------------------------------------------------------------------------
% COMPUTE Cl vs Mach (Karman-Tsien correction)
%
% Computes the lift coefficient for different Mach numbers using the
% Karman-Tsien compressibility correction applied to Cp.
%
% Inputs:
%   c          -> chord length
%   Mach_list  -> vector of Mach numbers
%   N          -> number of panels
%   nv         -> normal vectors (Nx2)
%   pl         -> panel lengths (Nx1)
%   alpha_deg  -> angle of attack (degrees)
%   cp         -> incompressible Cp (Nx1)
%
% Output:
%   Cl_Mach    -> lift coefficient for each Mach number
% -------------------------------------------------------------------------

cp = cp(:); % ensure column vector

% --- Convert angle to radians
alpha = deg2rad(alpha_deg);

% --- Flow direction unit vector (perpendicular to freestream)
Kinf = [-sin(alpha), cos(alpha)];

Cl_Mach = zeros(length(Mach_list),1);

for j = 1:length(Mach_list)
    
    M = Mach_list(j);
    beta = sqrt(1 - M^2);
    
    Cl = 0;
    
    for i = 1:N
        
        Cp0 = cp(i);
        
        % --- Karman-Tsien correction
        Cp_KT = Cp0 / (beta + (M^2/(1+beta))*(Cp0/2));
        
        % --- Lift contribution from panel i
        Cl = Cl - Cp_KT * pl(i) * dot(nv(i,:), Kinf) / c;
        
    end
    
    Cl_Mach(j) = Cl;
    
end

end