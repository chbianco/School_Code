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

t_vec = [0, 0.1, 0.3, 0.7, 3];

figure(5); hold on;
xlabel('$x$')
ylabel('Sheet Concentration')
grid on 

for n = 1:length(t_vec)
C = zeros(1, length(ts));
t = t_vec(n);

for j = 1:j_max
    Cj = sqrt(2).*cos((j + 0.5).*xs.*pi).*...
        exp(-(pi^2)*((j+0.5)^2)*t) * ((2/(2*pi*j + pi))*2*sin(j*pi + pi/2) + ...
       sqrt(2)*(j*pi+pi/2)*sin(j*pi + pi/2)*(2 * exp(-t)...
       * (((4 * (pi^2) * (j^2) + 4 * (pi^2) * j + (pi^2) - 4) * exp(t) -...
       4 * (pi^2) * (j^2) - 4 * (pi^2) * j - (pi^2)) * exp((pi^2) * (j^2) *...
       t + (pi^2) * j * t + ((pi^2) * t) / 4) + 4 * exp(t))) /...
       (pi^2 * (16 * (pi^2) * j^4 + 32 * (pi^2) * (j^3) + (24 * (pi^2) - 16) *...
       (j^2) + (8 * pi^2 - 16) * j + (pi^2) - 4)));

    C = C + Cj;
end

plot(xs, C, 'DisplayName', ['$t$' '=' num2str(t_vec(n))])
end
legend
hold off
%% LEss bad 
f = 1/2 - 0.5*exp(-ts);
figure(6); hold on;
xlabel('$t$')
ylabel('f(t)')
grid on 
plot(f, ts);
hold off;
% 
% figure(7)
% xlabel('$t$')
% ylabel('C(0,t)')
% grid on 

%% 5 B
Ra = linspace(0,16,1000);

U1 = sqrt(Ra -1);
U2 = -sqrt(Ra -1);
y1 = sqrt(Ra -1)./Ra;
y2 = -sqrt(Ra -1)./Ra;
x1 = 1./Ra;
x2 = -1./Ra;

figure(7); hold on;
plot(Ra, U1)
plot(Ra, U2)
plot(Ra, y1)
plot(Ra, y2)
plot(Ra, x1)
plot(Ra, x2)
ylim([-1 1])
legend({'U1', 'U2', 'y1', 'y2', 'x1', 'x2'})


hold off;