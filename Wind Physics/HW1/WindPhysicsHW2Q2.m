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

%% Dry air
z = linspace(0, 2000, 1000);
vpt_d = (20 - 0.01.*z)./((1-2.25577*10^(-5).*z).^(5.25588));

figure(1)
xlabel('Virtual Potential Temperature')
ylabel('z (m)')
grid on 
hold on 

plot(vpt_d, z)

hold off

%% Wet air
vpt_w = ((20 - 0.01.*z)./((1-2.25577*10^(-5).*z).^(5.25588))).*(1 + (1/1000).*(20-10.*(z/2000)));
figure(2)
xlabel('Virtual Potential Temperature')
ylabel('z (m)')
grid on 
hold on 

plot(vpt_w, z)

hold off

%% Comparison 
temp = (20-0.01).*z; 

figure(3)
hold on
xlabel('Temperature Comparison')
ylabel('z (m)')
grid on 

plot(vpt_d, z)
plot(vpt_w, z)
plot(temp, z)

legend({ 'VPT without Humidity', 'VPT with Humidity', 'Raw Temp'}, 'Location', 'best')
hold off


