%% PREAMBLE
close all;
clear variables;
clc;

set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);

%% Solve
k = linspace(0, pi, 100);

hkp_pa = (3.*sin(k))./(2+cos(k));
hkp_cd = sin(k);

%% Plot
figure(1); clf; hold on; grid on
xlabel('hk')
ylabel('$hk^\prime$')
plot(k, hkp_pa)
plot(k, hkp_cd)
legend({'Pade', 'Central Difference'}, 'Location','northwest')
