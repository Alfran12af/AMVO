function [x, z, xc, zc, nv, tv, cp, Cl, Cm14, xcp, pl, vortex] = ...
    constantstrength_tandem(data, N, c1, c2, d, alpha, delta)

% -------------------------------------------------------------------------
% CONSTANT STRENGTH VORTEX PANEL METHOD - TANDEM AIRFOIL
%
% Computes the aerodynamic coefficients of a two-element airfoil
% (main element + aft element) using a constant-strength vortex panel method.
%
% Inputs:
%   data   -> airfoil coordinates (Nx3: index, x, z)
%   N      -> number of panels per element
%   c1     -> main airfoil chord
%   c2     -> aft airfoil chord
%   d      -> gap between elements
%   alpha  -> angle of attack (radians)
%   delta  -> deflection angle of aft element (radians)
%
% Outputs:
%   x, z   -> combined geometry coordinates
%   xc, zc -> control points
%   nv     -> normal vectors
%   tv     -> tangent vectors
%   cp     -> pressure coefficient distribution
%   Cl     -> lift coefficient
%   Cm14   -> moment coefficient about quarter chord
%   xcp    -> center of pressure
%   pl     -> panel lengths
%   vortex -> vortex strength distribution
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% GEOMETRY DEFINITION

% Base airfoil (normalized)
x_base = data(:,2);
z_base = data(:,3);

% --- Main element
x1 = c1 * x_base;
z1 = c1 * z_base;

% --- Aft element (rotation + translation)
R = [cos(delta) -sin(delta); sin(delta) cos(delta)];
coords2 = [x_base z_base] * R';

x2 = c2 * coords2(:,1) + (c1 + d);
z2 = c2 * coords2(:,2);

% --- Combine both elements
x = [x1; x2];
z = [z1; z2];

NT = 2*N; % total panels

ivel = zeros(NT,NT,2);

%% ------------------------------------------------------------------------
% CONTROL POINTS

xc = zeros(NT,1);
zc = zeros(NT,1);

% Main element
for i = 1:N
    xc(i) = (x(i) + x(i+1)) / 2;
    zc(i) = (z(i) + z(i+1)) / 2;
end

% Aft element
for i = 1:N
    idx = i + N;
    xc(idx) = (x(idx) + x(idx+1)) / 2;
    zc(idx) = (z(idx) + z(idx+1)) / 2;
end

%% ------------------------------------------------------------------------
% PANEL GEOMETRY

[nv, tv, pl] = nortanlen(NT, x, z);

Qinf = [cos(alpha), sin(alpha)];

%% ------------------------------------------------------------------------
% INFLUENCE MATRIX

a = zeros(NT, NT);
b = zeros(NT, 1);

for i = 1:NT
    
    % RHS
    b(i) = -dot(Qinf, tv(i,:));
    
    for j = 1:NT
        
        if i == j
            a(i,j) = -0.5;
        else
            
            % Panel j local reference frame
            sinalphaj = -tv(j,2);
            cosalphaj =  tv(j,1);
            plj = pl(j);

            % Transform control point i into panel j coordinates
            xloc = (xc(i)-x(j))*cosalphaj - (zc(i)-z(j))*sinalphaj;
            zloc = (xc(i)-x(j))*sinalphaj + (zc(i)-z(j))*cosalphaj;

            % Geometric terms
            [r1, r2, theta1, theta2] = rtheta(plj, xloc, zloc);

            % Induced velocity (local)
            wij = (1/(4*pi)) * log(r2^2 / r1^2);
            uij = (theta2 - theta1) / (2*pi);

            % Convert to global coordinates
            u = uij*cosalphaj + wij*sinalphaj;
            w = -uij*sinalphaj + wij*cosalphaj;

            % Influence coefficient
            a(i,j) = dot([u, w], tv(i,:));
        end
    end
end

%% ------------------------------------------------------------------------
% KUTTA CONDITIONS (one per element)

% Main element trailing edge
k1 = round(N/4);
a(k1,:) = 0;
a(k1,1) = 1;
a(k1,N) = 1;
b(k1) = 0;

% Aft element trailing edge
k2 = N + round(N/4);
a(k2,:) = 0;
a(k2,N+1) = 1;
a(k2,NT) = 1;
b(k2) = 0;

%% ------------------------------------------------------------------------
% SOLVE SYSTEM

vortex = a \ b;

% Smooth Kutta panels
vortex(k1) = (vortex(k1-1) + vortex(k1+1)) / 2;
vortex(k2) = (vortex(k2-1) + vortex(k2+1)) / 2;

%% ------------------------------------------------------------------------
% POST-PROCESSING

Cl   = 0;
Cm0  = 0;
Cm14 = 0;
cp   = zeros(NT,1);

c_tot = c1 + c2 + d;

for i = 1:NT
    
    % Lift contribution
    Cl = Cl + vortex(i) * pl(i) / (norm(Qinf) * c_tot);
    
    % Pressure coefficient
    cp(i) = 1 - (vortex(i)/norm(Qinf))^2;

    % Moments
    Cm14 = Cm14 + cp(i) * ((xc(i)-0.25) * pl(i) / c_tot^2);
    Cm0  = Cm0  + cp(i) * ( xc(i)        * pl(i) / c_tot^2);
    
end

Cl = 2 * Cl;

% Center of pressure
xcp = c_tot * (1/4 - Cm0 / (Cl * cos(alpha)));

end