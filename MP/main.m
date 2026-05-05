clc; clear; close all;
tic

%% ------------------------------------------------------------------------
% LOAD GEOMETRY

addpath('/Users/joandurobayona/Documents/Matlab/AMVO/MP/HQ_300')
addpath('/Users/joandurobayona/Documents/Matlab/AMVO/MP/NACA0012')

chord = 1;
alpha_deg = 0:2:8;
alpha = deg2rad(alpha_deg);
gamma = 1.4;

%% ------------------------------------------------------------------------
% CONVERGENCE

files = {'HQ300_16.txt','HQ300_32.txt','HQ300_64.txt','HQ300_128.txt','HQ300_256.txt','HQ300_512.txt'};
files2 = {'NACA_0012_N_32.txt','NACA_0012_N_64.txt','NACA_0012_N_128.txt','NACA_0012_N_256.txt','NACA_0012_N_512.txt'};
HQ300_Cl_conv = zeros(length(files),1);
NACA0012_Cl_conv = zeros(length(files2),1);

alpha_test = deg2rad(4);

for k = 1:length(files)

    data = load(files{k});
    Nloc = size(data,1)-1;

    [~,~,~,~,~,~,~,HQ300_Cl,~,~,~,~,~] = ...
    vortexPanel2D(data,Nloc,chord,alpha_test);

    HQ300_Cl_conv(k) = HQ300_Cl;

end

for j = 1:length(files2)

    data2 = load(files2{j});
    Nloc2 = size(data2,1)-1;

    [~,~,~,~,~,~,~,NACA0012_Cl,~,~,~,~,~] = ...
    vortexPanel2D(data2,Nloc2,chord,alpha_test);

    NACA0012_Cl_conv(j) = NACA0012_Cl;

end


%% ------------------------------------------------------------------------
% PART 1: HQ300 AIRFOIL

HQ300 = load('HQ300_512.txt');

N = size(HQ300,1)-1;


for i=1:length(alpha)

    [x,z,xc,zc,nv,tv,cp,Cl,Cd,Cm14,xcp,pl,vortex] = ...
    vortexPanel2D(HQ300, N, chord, alpha(i));

    cp_values(i,:) = cp;
    cl_values(i) = Cl;
    cd_values(i) = Cd;
    cm_values(i) = Cm14;
    xcp_values(i) = xcp;
    cp_min(i) = min(cp);

end

% --- regressions
regr_cl = polyfit(alpha, cl_values, 1);
cl_alpha = regr_cl(1);
cl_0 = regr_cl(2);

regr_cm = polyfit(alpha, cm_values, 1);
cm_alpha = regr_cm(1);

xac = 0.25 - cm_alpha / cl_alpha;

fprintf('\n=====================================================\n');
fprintf('PART 1: HQ300 AIRFOIL\n');
fprintf('=====================================================\n');

fprintf('\nAlpha (deg):     ');
fprintf('%8.1f', alpha_deg);

fprintf('\nCl:              ');
fprintf('%8.4f', cl_values);

fprintf('\nCm(c/4):         ');
fprintf('%8.4f', cm_values);

fprintf('\n\nAerodynamic center (x/c) = %.4f\n', xac);

%% ------------------------------------------------------------------------
% PART 2: CRITICAL MACH

alpha_Mcrit = deg2rad([0 2 4]);

M_inf = linspace(0.3,0.8,100);

[Cp_KT, Cp_star] = karmanTsienCorrection(N, alpha_Mcrit, M_inf, gamma, cp_values(1:3,:));

Mcrit = criticalMachNumber(M_inf, Cp_KT, Cp_star);

fprintf('\n=====================================================\n');
fprintf('PART 2: CRITICAL MACH (Karman-Tsien)\n');
fprintf('=====================================================\n');

fprintf('\nAlpha (deg)    Mcrit\n');

for i = 1:length(Mcrit)
    fprintf('%10.1f    %6.4f\n', alpha_Mcrit(i)*180/pi, Mcrit(i));
end

%% ------------------------------------------------------------------------
% PART 3: Cl vs Mach (alpha = 4º)

alpha_target = deg2rad(4);
idx_alpha = find(alpha_deg==4);

Mach_list = [Mcrit(idx_alpha)-0.15,...
             Mcrit(idx_alpha)-0.10,...
             Mcrit(idx_alpha)-0.05,...
             Mcrit(idx_alpha)];

[~,~,~,~,nv_ref,~,cp_ref,~,~,~,~,pl_ref,~] = ...
    vortexPanel2D(HQ300, N, chord, alpha_target);

Cl_Mach = liftVsMachKT(chord, Mach_list, N, nv_ref, pl_ref, rad2deg(alpha_target), cp_ref);

fprintf('\n=====================================================\n');
fprintf('PART 3: Cl vs Mach (alpha = 4 deg)\n');
fprintf('=====================================================\n');

fprintf('\nMach        Cl\n');

for i = 1:length(Mach_list)
    fprintf('%6.3f    %6.4f\n', Mach_list(i), Cl_Mach(i));
end

%% ------------------------------------------------------------------------
% PART 4: TANDEM AIRFOIL

NACA0012 = load('NACA_0012_N_512.txt');

alpha2 = deg2rad(alpha_deg);

c1 = 0.64;
c2 = 0.34;
d  = 0.02;

for i = 1:length(alpha2)

    [~,~,~,~,~,~,~,Cl,~,Cm14,~,~,~] = ...
    vortexPanelTandem(NACA0012, N, c1, c2, d, alpha2(i), 0);

    cl_tandem(i) = Cl;
    cm_tandem(i) = Cm14;

end

fprintf('\n=====================================================\n');
fprintf('PART 4: TANDEM AIRFOIL (delta = 0)\n');
fprintf('=====================================================\n');

fprintf('\nAlpha (deg):     ');
fprintf('%8.1f', alpha_deg);

fprintf('\nCl (tandem):     ');
fprintf('%8.4f', cl_tandem);

fprintf('\nCm(c/4):         ');
fprintf('%8.4f', cm_tandem);

%% ------------------------------------------------------------------------
% PART 5: DEFLECTION EFFECT

delta_deg = 0:4:16;
delta = deg2rad(delta_deg);
alpha_fixed = deg2rad(4);

for i = 1:length(delta)

    [~,~,~,~,~,~,~,Cl,~,Cm14,~,~,~] = ...
    vortexPanelTandem(NACA0012, N, c1, c2, d, alpha_fixed, delta(i));

    cl_delta(i) = Cl;
    cm_delta(i) = Cm14;

end

fprintf('\n=====================================================\n');
fprintf('PART 5: DEFLECTION EFFECT (alpha = 4 deg)\n');
fprintf('=====================================================\n');

fprintf('\nDelta (deg):     ');
fprintf('%8.1f', delta_deg);

fprintf('\nCl:              ');
fprintf('%8.4f', cl_delta);

fprintf('\nCm(c/4):         ');
fprintf('%8.4f', cm_delta);

% --------------------------------------------------
% 1. Cp distribution
figure;
for i = 1:length(alpha)
    subplot(2,3,i)
    plot(xc, cp_values(i,:),'LineWidth',2)
    set(gca,'YDir','reverse')
    xlabel('x/c'); ylabel('Cp')
    title(['\alpha = ' num2str(alpha_deg(i)) 'º'])
    grid on; axis square
end


% --------------------------------------------------
% 2. Polar: Cl vs Cd
figure;
plot(cd_values, cl_values,'-o','LineWidth',2)
xlabel('C_d'); ylabel('C_l')
title('Polar curve (HQ300)')
grid on


% --------------------------------------------------
% 3. Cl, Cd, Cm vs alpha
figure;
subplot(1,3,1)
plot(alpha_deg, cl_values,'-o','LineWidth',2)
xlabel('\alpha (deg)'); ylabel('C_l')
grid on

subplot(1,3,2)
plot(alpha_deg, cd_values,'-o','LineWidth',2)
xlabel('\alpha (deg)'); ylabel('C_d')
grid on

subplot(1,3,3)
plot(alpha_deg, cm_values,'-o','LineWidth',2)
xlabel('\alpha (deg)'); ylabel('C_m')
grid on


% --------------------------------------------------
% 4. Cp_min vs alpha
figure;
plot(alpha_deg, cp_min,'-o','LineWidth',2)
xlabel('\alpha (deg)'); ylabel('C_{p,min}')
title('Minimum Cp vs alpha')
grid on


% --------------------------------------------------
% 5. Vortex distribution (gamma)
figure;
f = 0.05;
idx = find(alpha_deg==4);

for i = 1:N
    vec = f*vortex(i)*nv(i,:);
    if vortex(i)<0
        quiver(xc(i),zc(i),-vec(1),-vec(2),'r'); hold on
    else
        quiver(xc(i),zc(i),vec(1),vec(2),'b'); hold on
    end
end
plot(x,z,'k')
axis equal
title('\Gamma distribution (\alpha = 4º)')


% --------------------------------------------------
% 6. Convergence HQ300
figure;
plot([16 32 64 128 256 512], HQ300_Cl_conv,'-o','LineWidth',2)
xlabel('Panels'); ylabel('C_l')
title('Convergence HQ300')
grid on


% --------------------------------------------------
% 7. Cp_min vs Mach (IMPORTANT)
figure;
plot(M_inf, Cp_KT(:,2),'LineWidth',2); hold on
plot(M_inf, Cp_star,'--','LineWidth',2)
set(gca,'YDir','reverse')
xlabel('Mach'); ylabel('Cp')
title('Cp_{min} vs Mach')
legend('Cp_{KT}','Cp^*')
grid on

% --------------------------------------------------
% 8. Polar Tandem
figure;
plot(cd_values, cl_values,'-o','LineWidth',2)
xlabel('C_d'); ylabel('C_l')
title('Polar Tandem')
grid on


% --------------------------------------------------
% 9. Cl, Cm vs alpha (tandem)
figure;
subplot(1,2,1)
plot(alpha_deg, cl_tandem,'-o','LineWidth',2)
xlabel('\alpha'); ylabel('C_l')
title('Tandem Lift')

subplot(1,2,2)
plot(alpha_deg, cm_tandem,'-o','LineWidth',2)
xlabel('\alpha'); ylabel('C_m')
title('Tandem Moment')


% --------------------------------------------------
% 10. Convergence NACA0012
figure;
plot([32 64 128 256 512], NACA0012_Cl_conv,'-o','LineWidth',2)
xlabel('Panels'); ylabel('C_l')
title('Convergence NACA0012')
grid on


% --------------------------------------------------
% 11. Cl, Cm vs delta
figure;
subplot(1,2,1)
plot(delta_deg, cl_delta,'-o','LineWidth',2)
xlabel('\delta'); ylabel('C_l')
title('Lift vs deflection')

subplot(1,2,2)
plot(delta_deg, cm_delta,'-o','LineWidth',2)
xlabel('\delta'); ylabel('C_m')
title('Moment vs deflection')

toc