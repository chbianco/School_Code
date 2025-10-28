%% Preamble
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

%% Exact vs Expansion

x = linspace(0,1, 1000);
T = besseli(1, 2*0.6.*sqrt(x))./(besseli(1, 2*0.6));

T_est = 1 + 0.6^2 .*(x-1);

figure(1); hold on;
plot(x, T);
plot(x, T_est);
grid on
xlabel('$\tilde{x}$')
ylabel('$\tilde{T}$')
legend({'Exact Solution', 'Expansion Solution'}, 'Location','best')
hold off

%% Exact vs Boundary

T_bound = exp(-0.6.*(1-x));

figure(2); hold on;
plot(x, T);
plot(x, T_bound);
grid on
xlabel('$\tilde{x}$')
ylabel('$\tilde{T}$')
legend({'Exact Solution', 'Inner Solution'}, 'Location','best')
hold off