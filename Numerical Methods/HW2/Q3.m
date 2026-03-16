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

%% Solving
h = 1.2;
t = 0:h:15;
v = zeros(1, length(t));

f = @(t, V) -3.*t.*V./(1+t) + (2.*(1+t).^3).*exp(-t);
alp = @(t) 3.*t./(1+t); bet = @(t) (2.*(1+t).^3).*exp(-t);

exact = exp(-t).*(t + 1).^3;
v(1) = 1;

v_fe = v;
v_be = v;
v_trp = v;
v_rk2 = v;
v_rk4 = v; 

for j = 2:length(t)
%Forward Euler 
v_fe(j) = v_fe(j-1) +h.*f(t(j-1), v_fe(j-1));

%Backward/explicit Euler
v_be(j) = (v_be(j-1) + h.*((2.*(1+t(j-1)).^3).*exp(-t(j-1))))./(1 + h.*3.*t(j-1)./(1+t(j-1)));

%Trapezoidal
v_trp(j) = (v_trp(j-1) + (h/2)*(f(t(j-1), v_trp(j-1)) + bet(t(j))))/(1+(h/2)*alp(t(j)));

%RK2
foo = v_rk2(j-1) + (h/2)*f(t(j-1),v_rk2(j-1));
v_rk2(j) = v_rk2(j-1) + h*f(t(j-1) + h/2, foo);

%RK4
k1 = h*f(t(j-1), v_rk4(j-1));
k2 = h*f(t(j-1) + h/2, v_rk4(j-1) + 0.5*k1);
k3 = h*f(t(j-1) + h/2, v_rk4(j-1) + 0.5*k2);
k4 = h*f(t(j-1) + h, v_rk4(j-1) + k3);
v_rk4(j) = v_rk4(j-1) + (1/6)*k1 + (1/3)*(k2+k3)+(1/6)*k4;

end

%% Plotting
figure(1); hold on; grid on;
ylim([-1 5])
plot(t, exact)
plot(t, v_fe)
plot(t, v_be)
plot(t, v_trp)
plot(t, v_rk2)
plot(t, v_rk4)
title('h = 1.2')
legend({'Exact Solution', 'Explicit Euler', 'Backward Euler', 'Trapezoidal', 'RK2', 'RK4'})