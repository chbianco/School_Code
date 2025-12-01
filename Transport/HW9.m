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

figure(1);
plot(t, theta);
xlabel('$t$')
ylabel('$\theta(0.5,t)$')
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
figure(2); hold on;
plot(x, theta);
plot(x, theta_sim)
xlabel('$x$')
ylabel('$\theta(x,0.1)$')
legend({'Series solution', 'Similarity (short time) solution'})
grid on 
hold off