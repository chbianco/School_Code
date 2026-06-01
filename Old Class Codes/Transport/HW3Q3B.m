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

cplus = exp(-4.*exp(-x));
cminus = exp(4.*exp(-x));

figure(1); hold on;

plot(x, cplus);
plot(x,cminus);

xlabel('$x/\xi$')
ylabel('$C/C_0$')
legend({'$K^+$', '$Cl^-$'}, 'Location', 'best')
grid on 


