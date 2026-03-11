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
set(groot, 'defaultLineLineWidth', 1);

%% RK4 Solver
function [t, U] = rk4_solve(f, t0, tf, u0, h)
    t = (t0:h:tf)';
    N = length(t);
    U = zeros(N, numel(u0));
    U(1,:) = u0(:)';
    for n = 1:N-1
        k1 = f(U(n,:)');
        k2 = f(U(n,:)' + h/2*k1);
        k3 = f(U(n,:)' + h/2*k2);
        k4 = f(U(n,:)' + h*k3);
        U(n+1,:) = U(n,:) + (h/6)*(k1 + 2*k2 + 2*k3 + k4)';
    end
end

sig = 10; b = 8/3;

%% Part a
r = 20;
h = 0.005;
f20 = @(u) [sig*(u(2)-u(1));  r*u(1)-u(2)-u(1)*u(3);  u(1)*u(2)-b*u(3)];
[t, U20] = rk4_solve(f20, 0, 25, [1;1;1], h);
x=U20(:,1); y=U20(:,2); z=U20(:,3);

figure('Name','Part (a): r=20');
subplot(2,3,1); plot(x,y,'b'); xlabel('x'); ylabel('y'); title('xy plane'); grid on;
subplot(2,3,2); plot(x,z,'b'); xlabel('x'); ylabel('z'); title('xz plane'); grid on;
subplot(2,3,3); plot(y,z,'b'); xlabel('y'); ylabel('z'); title('yz plane'); grid on;
subplot(2,3,4); plot(t,x,'r'); xlabel('t'); ylabel('x'); grid on;
subplot(2,3,5); plot(t,y,'g'); xlabel('t'); ylabel('y'); grid on;
subplot(2,3,6); plot(t,z,'b'); xlabel('t'); ylabel('z'); grid on;
sgtitle('r=20');

%% Part b
r = 28;
h = 0.005;
f28 = @(u) [sig*(u(2)-u(1));  r*u(1)-u(2)-u(1)*u(3);  u(1)*u(2)-b*u(3)];
[t, U28] = rk4_solve(f28, 0, 25, [1;1;1], h);
x=U28(:,1); y=U28(:,2); z=U28(:,3);

figure('Name','Part (a): r=20');
subplot(2,3,1); plot(x,y,'b'); xlabel('x'); ylabel('y'); title('xy plane'); grid on;
subplot(2,3,2); plot(x,z,'b'); xlabel('x'); ylabel('z'); title('xz plane'); grid on;
subplot(2,3,3); plot(y,z,'b'); xlabel('y'); ylabel('z'); title('yz plane'); grid on;
subplot(2,3,4); plot(t,x,'r'); xlabel('t'); ylabel('x'); grid on;
subplot(2,3,5); plot(t,y,'g'); xlabel('t'); ylabel('y'); grid on;
subplot(2,3,6); plot(t,z,'b'); xlabel('t'); ylabel('z'); grid on;
sgtitle('r=28');

figure('Name','Part (b): r=28 3D');
plot3(z,y,x,'b','LineWidth',0.3); 
xlabel('z'); ylabel('y'); zlabel('x');
title('r=28: 3D trajectory'); grid on; view(30,30);

%% Part c
[t, Ua] = rk4_solve(f28, 0, 25, [6;6;6], h);
[t, Ub] = rk4_solve(f28, 0, 25, [6;6.01;6], h);

figure('Name','Part (c): sensitivity');
subplot(3,1,1);
plot(t, Ua(:,1), 'b', t, Ub(:,1), 'r--');
xlabel('t'); ylabel('x'); legend('(6,6,6)','(6,6.01,6)'); grid on;
title('r=28');

subplot(3,1,2);
plot(t, Ua(:,2), 'b', t, Ub(:,2), 'r--', 'LineWidth', 1);
xlabel('t'); ylabel('y'); grid on;

subplot(3,1,3);
plot(t, Ua(:,3), 'b', t, Ub(:,3), 'r--', 'LineWidth', 1);
xlabel('t'); ylabel('z'); grid on;