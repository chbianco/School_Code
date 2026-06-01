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

%% Solving the ODE
%Exact
eps = 0.1;
Km = 1;
gam = 0.5;

t_span = [0, 0.3];

%y(1) is s, y(2) is c
y0 = [1; 0]; 

ode = @(t, y) [-y(1) + (Km + y(1) - gam)*y(2);
    (-(y(1) + Km)*y(2) + y(1))/eps];


[t, y] = ode45(ode, t_span, y0);

s = y(:,1);
c = y(:,2);

%Approximate
ode_a = @(t, s) -gam * s / (Km + s);

[t_a, s_a] = ode45(ode_a, t_span, 1);
c_a = s_a ./ (Km + s_a);

%% Plotting
%Numerical Solution
figure(1); hold on;
plot(t, s);
plot(t, c);
grid on
xlabel('$\tilde{t}$')
ylabel('Non-Dim Concentration')
legend({'s', 'c'}, 'Location','best')
title('Numerical Solution')
hold off

%Approximate solution
figure(2); hold on;
plot(t_a, s_a);
plot(t_a, c_a);
grid on
xlabel('$\tilde{t}$')
ylabel('Non-Dim Concentration')
ylim([0 1]);
legend({'s', 'c'}, 'Location','best')
title('Approximate Solution (3.4)')
hold off
