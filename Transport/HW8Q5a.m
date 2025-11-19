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
r = linspace(-1,1,1000);
c = 1./sqrt(r.^2 + 1);

figure(1);
plot(r, c);
xlabel('$r$')
ylabel('$C$')
grid on 
hold off

%% Plotting for b 

flux = 1./(2*pi.*(1+r.^2).^(3/2));

figure(2);
plot(r, flux);
xlabel('$r$')
ylabel('Surface Flux')
grid on 
hold off
