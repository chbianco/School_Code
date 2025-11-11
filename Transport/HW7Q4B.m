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

%% Plotting
x = linspace(-1,1,1000);
ts = [0.01, 0.05, 1];

figure(1); hold on; %Initialize figure for plotting

for n = 1:3
T = 0.5*(1+erf(x./(2*sqrt(ts(n)))));
plot(x,T, 'DisplayName', ['$t$' '=' num2str(ts(n))]);


end

xlabel('$x$')
ylabel('$T^*$')
legend('Location', 'best')
grid on 
hold off

