function [x, z, xc, zc, nv, tv, cp, Cl, Cm14, xcp, pl, vortex] = ...
    constantstrength_project(data, N, c, alpha)

% -------------------------------------------------------------------------
% CONSTANT STRENGTH VORTEX PANEL METHOD
%
% Computes aerodynamic coefficients of an airfoil using a vortex panel method.
%
% Inputs:
%   data  -> airfoil coordinates (Nx3: index, x, z)
%   N     -> number of panels
%   c     -> chord length
%   alpha -> angle of attack (radians)
%
% Outputs:
%   cp    -> pressure coefficient distribution (Nx1)
%   Cl    -> lift coefficient
%   Cm14  -> moment coefficient about quarter chord
%   xcp   -> center of pressure
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% GEOMETRY

x = data(:,2);
z = data(:,3);

% Control points
xc = (x(1:end-1) + x(2:end)) / 2;
zc = (z(1:end-1) + z(2:end)) / 2;

% Panel geometry
[nv, tv, pl] = nortanlen(N, x, z);

% Freestream
Qinf = [cos(alpha), sin(alpha)];
Q = norm(Qinf);

%% ------------------------------------------------------------------------
% INITIALIZATION

a = zeros(N,N);          % Influence matrix
b = zeros(N,1);          % RHS
ivel = zeros(N,N,2);     % Induced velocity storage

%% ------------------------------------------------------------------------
% BUILD INFLUENCE MATRIX

for i = 1:N

    % No-penetration condition
    b(i) = -dot(Qinf, tv(i,:));

    for j = 1:N

        if i == j
            a(i,j) = -0.5; % self influence
        else
            
            % Panel orientation
            cos_a = tv(j,1);
            sin_a = -tv(j,2);

            % Transform control point to local panel coordinates
            dx = xc(i) - x(j);
            dz = zc(i) - z(j);

            x_local =  dx*cos_a - dz*sin_a;
            z_local =  dx*sin_a + dz*cos_a;

            % Geometric integrals
            [r1, r2, theta1, theta2] = rtheta(pl(j), x_local, z_local);

            % Induced velocity (local)
            w = (1/(4*pi)) * log(r2^2 / r1^2);
            u = (theta2 - theta1) / (2*pi);

            % Transform back to global coordinates
            vx =  u*cos_a + w*sin_a;
            vz = -u*sin_a + w*cos_a;

            ivel(i,j,1) = vx;
            ivel(i,j,2) = vz;

            % Influence coefficient
            a(i,j) = dot([vx, vz], tv(i,:));
        end
    end
end

%% ------------------------------------------------------------------------
% KUTTA CONDITION

idx = round(N/4);

a(idx,:) = 0;
b(idx) = 0;

a(idx,1) = 1;
a(idx,N) = 1;

%% ------------------------------------------------------------------------
% SOLVE SYSTEM

vortex = a \ b;

% Smooth Kutta panel
vortex(idx) = 0.5 * (vortex(idx-1) + vortex(idx+1));

%% ------------------------------------------------------------------------
% POST-PROCESSING

cp   = zeros(N,1);
Cl   = 0;
Cm14 = 0;
Cm0  = 0;

for i = 1:N

    % Lift contribution
    Cl = Cl + vortex(i) * pl(i) / (Q * c);

    % Induced velocity
    vind = [0,0];
    for j = 1:N
        vind = vind + vortex(j) * squeeze(ivel(i,j,:))';
    end

    V = Qinf + vind;

    % Pressure coefficient (simplified vortex relation)
    cp(i) = 1 - (vortex(i)/Q)^2;

    % Geometry increments
    dx = x(i+1) - x(i);
    dz = z(i+1) - z(i);

    % Moments
    Cm14 = Cm14 + cp(i)*((xc(i)-0.25)*dx/c^2 + zc(i)*dz/c^2);
    Cm0  = Cm0  + cp(i)*( xc(i)*dx/c^2        + zc(i)*dz/c^2);

end

Cl = 2 * Cl;

%% ------------------------------------------------------------------------
% CENTER OF PRESSURE

xcp = c * (1/4 - Cm0/(Cl*cos(alpha)));

end