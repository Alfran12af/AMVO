function [x, z, xc, zc, nv, tv, cp, Cl, Cd, Cm14, xcp, pl, vortex] = ...
    vortexPanelTandem(data, N, c1, c2, d, alpha, delta)

% -------------------------------------------------------------------------
% 2D VORTEX PANEL METHOD - TANDEM AIRFOIL
%
% Computes aerodynamic coefficients of a two-element airfoil
% using a constant-strength vortex panel method.
%
% Outputs include Cl, Cd (pressure drag ~0), Cp, and moments.
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% GEOMETRY

x_base = data(:,2);
z_base = data(:,3);

% Main element
x1 = c1 * x_base;
z1 = c1 * z_base;

% Aft element (rotation + translation)
R = [cos(delta) -sin(delta); sin(delta) cos(delta)];
coords2 = [x_base z_base] * R';

x2 = c2 * coords2(:,1) + (c1 + d);
z2 = c2 * coords2(:,2);

% Combine
x = [x1; x2];
z = [z1; z2];

NT = 2*N;

%% ------------------------------------------------------------------------
% CONTROL POINTS

xc = zeros(NT,1);
zc = zeros(NT,1);

for i = 1:NT
    xc(i) = (x(i) + x(i+1)) / 2;
    zc(i) = (z(i) + z(i+1)) / 2;
end

%% ------------------------------------------------------------------------
% PANEL GEOMETRY

[nv, tv, pl] = panelGeometry2D(NT, x, z);

Qinf = [cos(alpha), sin(alpha)];
Q = norm(Qinf);

%% ------------------------------------------------------------------------
% INFLUENCE MATRIX

a = zeros(NT, NT);
b = zeros(NT, 1);
ivel = zeros(NT, NT, 2);

for i = 1:NT
    
    b(i) = -dot(Qinf, tv(i,:));
    
    for j = 1:NT
        
        if i == j
            a(i,j) = -0.5;
        else
            
            cos_a = tv(j,1);
            sin_a = -tv(j,2);

            dx = xc(i) - x(j);
            dz = zc(i) - z(j);

            xloc = dx*cos_a - dz*sin_a;
            zloc = dx*sin_a + dz*cos_a;

            [r1, r2, theta1, theta2] = panelInfluenceGeometry(pl(j), xloc, zloc);

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
% KUTTA CONDITIONS

k1 = round(N/4);
k2 = N + round(N/4);

% Main airfoil
a(k1,:) = 0;
a(k1,1) = 1;
a(k1,N) = 1;
b(k1) = 0;

% Aft airfoil
a(k2,:) = 0;
a(k2,N+1) = 1;
a(k2,NT) = 1;
b(k2) = 0;

%% ------------------------------------------------------------------------
% SOLVE

vortex = a \ b;

vortex(k1) = 0.5*(vortex(k1-1) + vortex(k1+1));
vortex(k2) = 0.5*(vortex(k2-1) + vortex(k2+1));

%% ------------------------------------------------------------------------
% POST-PROCESSING

cp   = zeros(NT,1);
Cl   = 0;
Cd   = 0;
Cm14 = 0;
Cm0  = 0;

c_tot = c1 + c2 + d;

% Unit vectors
iinf = [cos(alpha), sin(alpha)];
kinf = [-sin(alpha), cos(alpha)];

for i = 1:NT
    
    % ------------------------------
    % Velocity
    vind = [0,0];
    for j = 1:NT
        vind = vind + vortex(j) * squeeze(ivel(i,j,:))';
    end

    V = Qinf + vind;

    % ------------------------------
    % Cp (CORRECT)
    cp(i) = 1 - (norm(V)/Q)^2;

    % ------------------------------
    % Lift & Drag
    Cl = Cl - cp(i) * pl(i) * dot(nv(i,:), kinf) / c_tot;
    Cd = Cd - cp(i) * pl(i) * dot(nv(i,:), iinf) / c_tot;

    % ------------------------------
    % Moments
    Cm14 = Cm14 - cp(i) * pl(i) * ...
           ( (xc(i)-0.25)*nv(i,2) - zc(i)*nv(i,1) ) / c_tot^2;

    Cm0  = Cm0 - cp(i) * pl(i) * ...
           ( xc(i)*nv(i,2) - zc(i)*nv(i,1) ) / c_tot^2;

end

%% ------------------------------------------------------------------------
% CENTER OF PRESSURE

xcp = c_tot * (1/4 - Cm0/(Cl*cos(alpha)));

end