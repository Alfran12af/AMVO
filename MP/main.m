clc; clear; close all;
tic

%% ------------------------------------------------------------------------
% LOAD GEOMETRY

addpath('/Users/joandurobayona/Downloads/code/HQ_300')

HQ300 = load('HQ300_128.txt');
% NACA0012 = load('NACA_0012_N_128.txt');

chord = 1;
alpha_deg = 0:2:8;
alpha = deg2rad(alpha_deg);
gamma = 1.4;

%% ------------------------------------------------------------------------
% PART 1: HQ300 AIRFOIL (Cl, Cm vs alpha + aerodynamic center)

N = size(HQ300,1)-1;


for i=1:length(alpha)

    [x,z,xc,zc,nv,tv,cp,Cl,Cm14,xcp,pl,vortex] = ...
        constantstrength_project(HQ300, N, chord, alpha(i));

    cp_values(i,:) = cp;
    cl_values(i) = Cl;
    cm_values(i) = Cm14;
    xcp_values(i) = xcp;

end

% --- Linear regression (Cl vs alpha)
regr_cl = polyfit(alpha, cl_values, 1);
cl_alpha = regr_cl(1);
cl_0 = regr_cl(2);

% --- Linear regression (Cm vs alpha)
regr_cm = polyfit(alpha, cm_values, 1);
cm_alpha = regr_cm(1);

% --- Aerodynamic center
xac = 0.25 - cm_alpha / cl_alpha;

fprintf('Aerodynamic center (x/c) = %.4f\n', xac);

%% ------------------------------------------------------------------------
% PART 2: CRITICAL MACH (KÁRMÁN-TSIEN)

M_inf = linspace(0.3,0.8,100);

% 🔴 NOVA FUNCIÓ (substitueix Laitone)
[Cp_KT, Cp_star] = KarmanTsien(N, alpha, M_inf, gamma, cp_values);

% 🔴 NOVA FUNCIÓ
Mcrit = compute_Mcrit(M_inf, Cp_KT, Cp_star);

%% ------------------------------------------------------------------------
% % PART 3: Cl vs Mach (alpha = 4º)
% 
% alpha_target = deg2rad(4);
% 
% % 🔴 NOVA FUNCIÓ
% idx_alpha = find(alpha_deg==4);
% 
% Mach_list = [Mcrit(idx_alpha)-0.15,...
%              Mcrit(idx_alpha)-0.10,...
%              Mcrit(idx_alpha)-0.05,...
%              Mcrit(idx_alpha)];
% 
% % Recompute for alpha = 4º
% [~,~,~,~,nv_ref,~,cp_ref,~,~,~,pl_ref,~] = ...
%     constantstrength_project(HQ300, N, chord, alpha_target);
% 
% Cl_Mach = compute_Cl_KT(chord, Mach_list, N, nv_ref, gamma, pl_ref, alpha_target, cp_ref);

%% ------------------------------------------------------------------------
% % PART 4: TWO-ELEMENT AIRFOIL (NACA 0012 TANDEM)
% 
% alpha_deg2 = 0:2:8;
% alpha2 = deg2rad(alpha_deg2);
% 
% c1 = 0.64;
% c2 = 0.34;
% d  = 0.02;
% 
% for i = 1:length(alpha2)
% 
%     % 🔴 REUTILITZEM però adaptat
%     [x,z,xc,zc,nv,tv,cp,Cl,Cm14,xcp,pl,vortex] = ...
%         constantstrength_tandem(NACA0012, N, c1, c2, d, alpha2(i), 0);
% 
%     cl_tandem(i) = Cl;
%     cm_tandem(i) = Cm14;
% 
% end

%% ------------------------------------------------------------------------
% % PART 5: EFFECT OF DEFLECTION (delta_e)
% 
% delta_deg = 0:4:16;
% delta = deg2rad(delta_deg);
% 
% alpha_fixed = deg2rad(4);
% 
% for i = 1:length(delta)
% 
%     [x,z,xc,zc,nv,tv,cp,Cl,Cm14,xcp,pl,vortex] = ...
%         constantstrength_tandem(NACA0012, N, c1, c2, d, alpha_fixed, delta(i));
% 
%     cl_delta(i) = Cl;
%     cm_delta(i) = Cm14;
% 
% end

%% ------------------------------------------------------------------------
% PLOTS

figure;
for i = 1:length(alpha)
    subplot(2,3,i)
    plot(xc, cp_values(i,:),'LineWidth',2)
    set(gca,'YDir','reverse')
    xlabel('x/c'); ylabel('Cp')
    title(['\alpha = ' num2str(alpha_deg(i)) 'º'])
    grid on; axis square
end

figure;

subplot(1,2,1)
plot(alpha_deg, cl_values,'-o','LineWidth',2)
xlabel('\alpha (º)'); ylabel('C_l')
title('Lift coefficient')
grid on; axis square

subplot(1,2,2)
plot(alpha_deg, cm_values,'-o','LineWidth',2)
xlabel('\alpha (º)'); ylabel('C_{m_{1/4}}')
title('Moment coefficient')
grid on; axis square

figure;
f = 0.05;

idx = 3; % alpha = 4º
[~,~,xc,zc,nv,~,~,~,~,~,~,vortex] = ...
    constantstrength_project(HQ300, N, chord, alpha(idx));

for i=1:N
    vec = f*vortex(i)*nv(i,:);
    if vortex(i)<0
        quiver(xc(i),zc(i),-vec(1),-vec(2),'r'); hold on
    else
        quiver(xc(i),zc(i),vec(1),vec(2),'b'); hold on
    end
end
plot(HQ300(:,2), HQ300(:,3),'k')
axis equal
title('\alpha = 4º vortex distribution')

figure;
plot(alpha_deg(1:length(Mcrit)), Mcrit,'-o','LineWidth',2)
xlabel('\alpha (º)')
ylabel('M_{crit}')
title('Critical Mach vs \alpha')
grid on

figure;
plot(Mach_list, Cl_Mach,'-o','LineWidth',2)
xlabel('M_{\infty}')
ylabel('C_l')
title('Lift vs Mach')
grid on

%% ------------------------------------------------------------------------
% CONVERGENCE STUDY (OPTIONAL)

files = {'HQ300_16.txt','HQ300_32.txt','HQ300_64.txt','HQ300_128.txt','HQ300_256.txt','HQ300_512.txt'};
Cl_conv = zeros(length(files),1);

alpha_test = deg2rad(4);

for k = 1:length(files)

    data = load(files{k});
    Nloc = size(data,1)-1;

    [~,~,~,~,~,~,~,Cl,~,~,~,~] = ...
        constantstrength_project(data,Nloc,chord,alpha_test);

    Cl_conv(k) = Cl;

end

figure;
plot([16 32 64 128 256 512], Cl_conv,'-o','LineWidth',2)
xlabel('Number of panels')
ylabel('C_l')
title('Convergence study (HQ300)')
grid on

toc