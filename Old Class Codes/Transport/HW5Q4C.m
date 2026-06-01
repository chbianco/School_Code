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

cin = 1-exp(-50.*x);
cout = 1./(0.2.*(1-x)+1);
ccomp = 1./(0.2.*(1-x)+1) -exp(-50.*x)/(1.2);

figure(1); hold on;

plot(x, cin);
plot(x, cout);
plot(x, ccomp);


xlabel('$x\L$')
ylabel('$C/C_0$')
grid on 
legend({'$C_{in}/C_0$', '$C_{out}/C_0$', '$C_{composite}/C_0$'})
hold off

