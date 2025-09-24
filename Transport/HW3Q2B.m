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
r = linspace(0,1,100);
lambdas = [0, 0.1, 1];

figure(1); hold on; %Initialize figure for plotting

for n = 1:3
v = r./4 - lambdas(n)/2 - 1/4;

plot(r,-v, 'DisplayName', sprintf('Slip Length = %d', lambdas(n)));


end

xlabel('$r^* = r/a$')
ylabel('$u^*$')
legend('Location', 'best')
grid on 
