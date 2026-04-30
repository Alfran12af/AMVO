function [x,z,xc,zc,nv,tv,cp,Cl,Cm14,xcp,pl,vortex] = constantstrength_project(data,N,c,alpha)

% -------------------------------------------------------------------------
% CONSTANT STRENGTH VORTEX PANEL METHOD
%
% Inputs:
%   data  -> airfoil coordinates
%   N     -> number of panels
%   c     -> chord
%   alpha -> angle of attack (rad)
%
% Outputs:
%   cp    -> pressure coefficient
%   Cl    -> lift coefficient
%   Cm14  -> moment coefficient at quarter chord
%   xcp   -> center of pressure
% -------------------------------------------------------------------------

%% GEOMETRY

x = data(:,2);
z = data(:,3);

% Control points
xc = (x(1:end-1) + x(2:end)) / 2;
zc = (z(1:end-1) + z(2:end)) / 2;

% Normals, tangents and panel lengths
function [nv,tv,pl] = nortanlen(N,x,z)
% Outputs
%   nv: normal vectors at each control point
%   tv: tangent vectors at each control point
%   pl: panel lengths

nv=zeros(N,2); 
tv=zeros(N,2); 
pl=zeros(N,1); 
nx=zeros(N,1); nz=zeros(N,1);
tx=zeros(N,1); tz=zeros(N,1);
 
for i=1:N
    pl(i)=sqrt((x(i+1)-x(i))^2+(z(i+1)-z(i))^2); % panel length
    nx(i)=(z(i)-z(i+1))/pl(i); % normal vector x components (sin alpha)
    nz(i)=(x(i+1) - x(i))/pl(i); % normal vector z components (cos alpha)
    nv(i,:)=[nx(i), nz(i)];
    
    tx(i)= nz(i); % tangent vector x components (cos alpha)
    tz(i)=-nx(i); % tangent vector z components (-sin alpha)
    tv(i,:)=[tx(i), tz(i)];

end
end[nv, tv, pl] = nortanlen(N, x, z);

% Freestream velocity
Qinf = [cos(alpha), sin(alpha)];
Q = norm(Qinf);

%% INITIALIZATION

a = zeros(N,N);          % Influence matrix
b = zeros(N,1);          % RHS vector
ivel = zeros(N,N,2);     % Induced velocities

%% BUILD INFLUENCE MATRIX

for i = 1:N

    % RHS (no penetration condition)
    b(i) = -dot(Qinf, tv(i,:));

    for j = 1:N

        if i == j
            % Self influence
            a(i,j) = -1/2;

        else
            % Panel orientation
            cos_a = tv(j,1);
            sin_a = -tv(j,2);

            % Local coordinates
            dx = xc(i) - x(j);
            dz = zc(i) - z(j);

            x_local =  dx*cos_a - dz*sin_a;
            z_local =  dx*sin_a + dz*cos_a;

            % Geometric terms
            [r1, r2, theta1, theta2] = rtheta(pl(j), x_local, z_local);

            % Induced velocity in local frame
            w = (1/(4*pi)) * log(r2^2 / r1^2);
            u = (theta2 - theta1) / (2*pi);

            % Transform to global frame
            vx =  u*cos_a + w*sin_a;
            vz = -u*sin_a + w*cos_a;

            ivel(i,j,1) = vx;
            ivel(i,j,2) = vz;

            % Influence coefficient
            a(i,j) = dot([vx, vz], tv(i,:));
        end
    end
end

%% KUTTA CONDITION

idx = round(N/4);

a(idx,:) = 0;
b(idx) = 0;
a(idx,1) = 1;
a(idx,N) = 1;

%% SOLVE SYSTEM

vortex = a \ b;

% Smooth eliminated panel
vortex(idx) = 0.5 * (vortex(idx-1) + vortex(idx+1));

%% COMPUTE Cp, Cl, Cm

cp = zeros(N,1);
Cl = 0;
Cm14 = 0;
Cm0 = 0;

for i = 1:N

    % Lift contribution
    Cl = Cl + vortex(i)*pl(i)/(Q*c);

    % Velocity at control point
    vind = [0,0];
    for j = 1:N
        vind = vind + vortex(j) * squeeze(ivel(i,j,:))';
    end

    v = Qinf + vind;

    % Pressure coefficient
    cp(i) = 1 - (vortex(i)/Q)^2;

    % Moment coefficients
    dx = x(i+1) - x(i);
    dz = z(i+1) - z(i);

    Cm14 = Cm14 + cp(i)*((xc(i)-0.25)*dx/c^2 + zc(i)*dz/c^2);
    Cm0  = Cm0  + cp(i)*( xc(i)*dx/c^2        + zc(i)*dz/c^2);

end

Cl = 2 * Cl;

%% CENTER OF PRESSURE

xcp = c * (1/4 - Cm0/(Cl*cos(alpha)));

end