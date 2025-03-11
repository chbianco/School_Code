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
set(groot, 'defaultFigureColor', 'white');

s = linspace(0, 10^-4, 1000);
v = linspace(0,30,1000);

 mu = 1.83e-5; % Dynamic viscosity, in Pa-s
 rho = 1000; % Density of the particle, kg/m^3
 r = 5*10^-5;
 u = 20;

 Res = rho.*u.*s./mu;
 Rev = rho.*v.*r./mu;

figure(1)
plot(s, Res);
xlabel('Particle Radius (m)')
ylabel('Re')

figure(2)
plot(v, Rev)
xlabel('Initial Particle Velocity (m/s)')
ylabel('Re')