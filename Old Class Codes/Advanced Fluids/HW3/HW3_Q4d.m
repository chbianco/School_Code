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

%% Plotting

s = linspace(0, 10^-4, 1000); %Vector of particle sizes
mu = 1.8*10^-5; % Dynamic viscosity of air, in Pa-s
rho = 1000; % Density of the particle, kg/m^3
u_vec = [10, 20 30]; %Particle velocities, in m/s

figure(1); hold on; %Initialize figure for plotting


for i = 1:length(u_vec)
    
    u = u_vec(i);
    Re = rho.*u.*s./mu;
    plot(s, Re, 'DisplayName', sprintf('$u_0 = %d$ m/s', u));

end

xlabel('Particle Radius (m)')
ylabel('Re')
xlim([0 10^-4])
legend('Location', 'best')
grid on 
