function [x, z, xc, zc, nv, tv, cp, Cl, Cd, Cm14, xcp, pl, vortex] = ...
    vortexPanel2D(data, N, c, alpha)

% -------------------------------------------------------------------------
% 2D VORTEX PANEL METHOD (CONSTANT STRENGTH)
%
% Computes aerodynamic coefficients of an airfoil using a constant-strength
% vortex panel method (incompressible potential flow).
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
%   Cd    -> drag coefficient (≈ 0 in potential flow)
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
[nv, tv, pl] = panelGeometry2D(N, x, z);

% Freestream
Qinf = [cos(alpha), sin(alpha)];
Q = norm(Qinf);

%% ------------------------------------------------------------------------
% INITIALIZATION

a = zeros(N,N);
b = zeros(N,1);
ivel = zeros(N,N,2);

%% ------------------------------------------------------------------------
% BUILD INFLUENCE MATRIX

for i = 1:N

    b(i) = -dot(Qinf, tv(i,:));

    for j = 1:N

        if i == j
            a(i,j) = -0.5;
        else
            
            cos_a = tv(j,1);
            sin_a = -tv(j,2);

            dx = xc(i) - x(j);
            dz = zc(i) - z(j);

            x_local =  dx*cos_a - dz*sin_a;
            z_local =  dx*sin_a + dz*cos_a;

            [r1, r2, theta1, theta2] = panelInfluenceGeometry(pl(j), x_local, z_local);

            w = (1/(4*pi)) * log(r2^2 / r1^2);
            u = (theta2 - theta1) / (2*pi);

            vx =  u*cos_a + w*sin_a;
            vz = -u*sin_a + w*cos_a;

            ivel(i,j,1) = vx;
            ivel(i,j,2) = vz;

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

% Smooth trailing edge
vortex(idx) = 0.5 * (vortex(idx-1) + vortex(idx+1));

%% ------------------------------------------------------------------------
% POST-PROCESSING

cp   = zeros(N,1);
Cl   = 0;
Cd   = 0;
Cm14 = 0;
Cm0  = 0;

% Unit vectors
iinf = [cos(alpha), sin(alpha)];        % flow direction
kinf = [-sin(alpha), cos(alpha)];       % lift direction

for i = 1:N

    % ------------------------------
    % Velocity at control point
    vind = [0,0];
    for j = 1:N
        vind = vind + vortex(j) * squeeze(ivel(i,j,:))';
    end

    V = Qinf + vind;

    % ------------------------------
    % Pressure coefficient
    cp(i) = 1 - (norm(V)/Q)^2;

    % ------------------------------
    % Lift & Drag from Cp
    Cl = Cl - cp(i) * pl(i) * dot(nv(i,:), kinf) / c;
    Cd = Cd - cp(i) * pl(i) * dot(nv(i,:), iinf) / c;

    % ------------------------------
    % Moment coefficients (consistent formulation)
    Cm14 = Cm14 - cp(i) * pl(i) * ...
           ( (xc(i)-0.25)*nv(i,2) - zc(i)*nv(i,1) ) / c^2;

    Cm0  = Cm0 - cp(i) * pl(i) * ...
           ( xc(i)*nv(i,2) - zc(i)*nv(i,1) ) / c^2;

end

%% ------------------------------------------------------------------------
% CENTER OF PRESSURE

xcp = c * (1/4 - Cm0/(Cl*cos(alpha)));

end