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

%% Solve
% Parameters
lam = 2.74;   
om = 5.48;   
b0 = pi/2; 
xi = 0.96;
eta = 0.24;

%RHS 
function dudt = rhs(u, lam, om, b0, xi, eta)
    al=u(1); dal=u(2); ga=u(3); dga=u(4);
    C = cos(b0)*cos(al-ga) - sin(b0)*sin(al-ga);
    S = sin(b0)*cos(al-ga) + cos(b0)*sin(al-ga);
    D = 1 - xi*eta*C^2;
    rhs_a = -lam^2*sin(al) - eta*S*dga^2;  
    rhs_g = -om^2*sin(ga)  + xi*S*dal^2;    
    ddal = (rhs_a - eta*C*rhs_g) / D;      
    ddga = (rhs_g - xi*C*rhs_a)  / D;    
    dudt = [dal; ddal; dga; ddga];
end

%RK4
function [t, U] = rk4_solve(f, t0, tf, u0, h)
    t = (t0:h:tf)';
    N = length(t);
    U = zeros(N, length(u0));
    U(1,:) = u0(:)';
    for n = 1:N-1
        k1 = f(U(n,:)');
        k2 = f(U(n,:)' + h/2*k1);
        k3 = f(U(n,:)' + h/2*k2);
        k4 = f(U(n,:)' + h*k3);
        U(n+1,:) = U(n,:) + (h/6)*(k1+2*k2+2*k3+k4)';
    end
end

f = @(u) rhs(u, lam, om, b0, xi, eta);
u0 = [pi/12; pi; 0; 0];
[t, U] = rk4_solve(f, 0, 100, u0, 0.005);

%% Plot 
figure;
subplot(2,1,1);
plot(t, U(:,1)*180/pi, 'b'); hold on;
plot(t, U(:,3)*180/pi, 'r');
xlabel('t (s)'); ylabel('Angle (deg)');
legend('$\alpha$', '$\gamma$'); grid on;

subplot(2,1,2);
plot(t, U(:,2), 'b'); hold on;
plot(t, U(:,4), 'r');
xlabel('t (s)'); ylabel('Angular velocity (rad/s)');
legend('$\dot\alpha$', '$\dot\gamma$'); grid on;


%% Chaotic ICs
u0_chaos1 = [pi/2; 5; 0; 0];
u0_chaos2 = [pi/2; 5.025; 0; 0];  

[t2, U2a] = rk4_solve(f, 0, 50, u0_chaos1, 0.0005);
[t2, U2b] = rk4_solve(f, 0, 50, u0_chaos2, 0.0005);

figure;
subplot(2,1,1);
plot(t2, U2a(:,1)*180/pi, 'b'); hold on;
plot(t2, U2b(:,1)*180/pi, 'r--');
xlabel('t (s)'); ylabel('$\alpha$ (deg)');
legend('$\dot{\alpha}_0=5$','$\dot{\alpha}_0=5.025$');
title('$\alpha$'); grid on;

subplot(2,1,2);
plot(t2, U2a(:,3)*180/pi, 'b'); hold on;
plot(t2, U2b(:,3)*180/pi, 'r--');
xlabel('t (s)'); ylabel('$\gamma$ (deg)');
legend('$\dot{\alpha}_0=5$','$\dot{\alpha}_0=5.025$');
title('$\gamma$'); grid on;


