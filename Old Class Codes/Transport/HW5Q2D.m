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
x = linspace(0, 1, 1000);

c1 = (exp(1) - exp(1.*x))./(exp(1) -1);
c2 = (exp(30) - exp(30.*x))./(exp(30) -1);

figure(1); hold on;

plot(x, c1);
plot(x, c2);

xlabel('$x\L$')
ylabel('$C/C_0$')
grid on 
legend({'P=1', 'P=30'})
hold off

