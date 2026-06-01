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

%% Plotting

[zR, zI] = meshgrid(linspace(-4, 1, 800), linspace(-4, 4, 800));
Z = zR + 1i*zI;

sigma_euler = abs(1 + Z);
sigma_rk2   = abs(1 + Z + Z.^2/2);
sigma_rk3   = abs(1 + Z + Z.^2/2 + Z.^3/6);
sigma_rk4   = abs(1 + Z + Z.^2/2 + Z.^3/6 + Z.^4/24);
sigma_ab2   = abs((Z .* (Z - 1)) ./ (0.5*(3*Z - 1)));  

figure; hold on; grid on; axis equal;

contour(zR, zI, sigma_euler, [1 1], 'k',  'LineWidth', 2);
contour(zR, zI, sigma_rk2,   [1 1], 'b',  'LineWidth', 2);
contour(zR, zI, sigma_rk3,   [1 1], 'm',  'LineWidth', 2);
contour(zR, zI, sigma_rk4,   [1 1], 'r',  'LineWidth', 2);

theta = linspace(0, 2*pi, 2000);
eit   = exp(1i*theta);
z_ab2 = (eit .* (eit - 1)) ./ (0.5*(3*eit - 1));
plot(real(z_ab2), imag(z_ab2), 'g', 'LineWidth', 2);

xline(0,'k--','LineWidth',0.5); yline(0,'k--','LineWidth',0.5);
xlabel('$\lambda_R \Delta t$'); ylabel('$\lambda_I \Delta t$');
title('Stability Boundaries');
legend('Euler','RK2','RK3','RK4','AB2','Location','northwest');
xlim([-4 1]); ylim([-4 4]);