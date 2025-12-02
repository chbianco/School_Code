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


%% 2H
t = linspace(0,1,1000);
 
k_max = 1000;
theta_h = zeros(1,length(t));

for k = 1:k_max
theta_k = (2.*(sin(pi*k) - pi*k)./(pi^2*k^2)).*exp(-pi^2*k^2.*t).*sin(pi*k*0.5);
theta_h = theta_h + theta_k;
end
%% 2H Plotting

theta= 0.5 + theta_h;

theta_sim1 = erfc(0.5./(2.*sqrt(t)));

figure(1);
plot(t, theta);
xlabel('$t$')
ylabel('$\theta(0.5,t)$')
grid on 

figure(2);
plot(t, theta_sim1);
xlabel('$t$')
ylabel('$\theta(0.5,t)$')
title('Similarity Solution')
grid on 



%% 2I
x = linspace(0,1,1000);
k_max = 5;
theta_h = zeros(1,length(t));

for k = 1:k_max
theta_k = (2.*(sin(pi*k) - pi*k)./(pi^2*k^2)).*exp(-pi^2*k^2.*0.1).*sin(pi*k.*x);
theta_h = theta_h + theta_k;
end

theta_sim = erfc(x./(2*sqrt(0.1)));
theta = 1 - x + theta_h;

figure(3); hold on;
plot(x, theta);
plot(x, theta_sim)
xlabel('$x$')
ylabel('$\theta(x,0.1)$')
legend({'Series solution', 'Similarity (short time) solution'})
grid on 
hold off

%% 2J

F_0t = zeros(1, length(t));
F_1t = zeros(1, length(t));

for k = 1:k_max
F_0t_k = sqrt(2).*(sin(pi*k)-pi*k).*exp(-(pi^2)*(k^2).*t);
F_1t_k =  sqrt(2).*(sin(pi*k)-pi*k).*exp(-(pi^2)*(k^2).*t).*cos(pi*k);

F_0t = F_0t + F_0t_k;
F_1t = F_1t + F_1t_k;
end

F_0 = 1 - F_0t;
F_1 = 1 - F_1t;


figure(4); hold on;
plot(t, F_0);
plot(t, F_1)
xlabel('$t$')
ylabel('Mass Flux')
legend({'x = 0', 'x = 1'})
grid on 
hold off

%% 3B
j_max = 1000;
ts = linspace(0, 6, 1000);
xs = linspace(0, 1, 1000);

c0 = zeros(1, length(xs));
c01 = zeros(1, length(xs));
c03 = zeros(1, length(xs));
c07 = zeros(1, length(xs));
c3 = zeros(1, length(xs));

for j = 1:j_max
c0_j = sqrt(2).*cos((j + 1/2).*pi.*xs).*(sqrt(2)/2)*(j*pi + pi/2)*sin(j*pi + pi/2).*exp((-j^2).*(pi^2).*0).*(exp(-0).*(0.*exp(0) + 1) -1);
c01_j = sqrt(2).*cos((j + 1/2).*pi.*xs).*(sqrt(2)/2)*(j*pi + pi/2)*sin(j*pi + pi/2).*exp((-j^2).*(pi^2).*0.1).*(exp(-0.1).*(0.1.*exp(0.1) + 1) -1);
c03_j = sqrt(2).*cos((j + 1/2).*pi.*xs).*(sqrt(2)/2)*(j*pi + pi/2)*sin(j*pi + pi/2).*exp((-j^2).*(pi^2).*0.3).*(exp(-0.3).*(0.3.*exp(0.3) + 1) -1);
c07_j = sqrt(2).*cos((j + 1/2).*pi.*xs).*(sqrt(2)/2)*(j*pi + pi/2)*sin(j*pi + pi/2).*exp((-j^2).*(pi^2).*0.7).*(exp(-0.7).*(0.7.*exp(0.7) + 1) -1);
c3_j = sqrt(2).*cos((j + 1/2).*pi.*xs).*(sqrt(2)/2)*(j*pi + pi/2)*sin(j*pi + pi/2).*exp((-j^2).*(pi^2).*3).*(exp(-3).*(3.*exp(3) + 1) -1);

c0 = c0 + c0_j;
c01 = c01 + c01_j;
c03 = c03 + c03_j;
c07 = c07 + c07_j;
c3 = c3 + c3_j;
end

figure(5); hold on;
plot(xs, 1 + c0);
plot(xs, 1 + c01);
plot(xs, 1 + c03);
plot(xs, 1 + c07);
plot(xs, 1 + c3);
xlabel('$x$')
ylabel('Sheet Concentration')
legend({'t = 0', 't = 0.1', 't = 0.3', 't = 0.7', 't = 3'})
grid on 
hold off