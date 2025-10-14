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
x = linspace(0, 10, 1000);

theta = (0.1.*exp(x))./(-0.9.*exp(-x)+1-0.1+0.1.*exp(x));

figure(1); hold on;

plot(x, theta);
xlabel('$Pe$')
ylabel('$\theta$')
grid on 


