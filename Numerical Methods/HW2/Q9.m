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

%% Setup
%Parameters
a = @(x) -(x+3)./(x+1);
b = @(x)  (x+3)./(x+1).^2;
f = @(x)  2*(x+1) + 3*b(x);

L  = 2;
TA = 5;   % T(0)
TB = 4;   % T(2)

%RK4
function [x, U] = rk4_bvp(rhs, x0, xf, u0, h)
    x = (x0:h:xf)';
    N = length(x);
    U = zeros(N, numel(u0));
    U(1,:) = u0(:)';
    for n = 1:N-1
        k1 = rhs(x(n),       U(n,:)');
        k2 = rhs(x(n)+h/2,   U(n,:)' + h/2*k1);
        k3 = rhs(x(n)+h/2,   U(n,:)' + h/2*k2);
        k4 = rhs(x(n)+h,     U(n,:)' + h*k3);
        U(n+1,:) = U(n,:) + (h/6)*(k1+2*k2+2*k3+k4)';
    end
end

h = 0.01;

rhs_full = @(x, u) [u(2); f(x) - a(x)*u(2) - b(x)*u(1)];
rhs_hom  = @(x, u) [u(2); -a(x)*u(2) - b(x)*u(1)];

%% aii
% Shoot u: full eqn, u(0)=TA, u'(0)=0
[xg, U_u] = rk4_bvp(rhs_full, 0, L, [TA; 0], h);
% Shoot v: homogeneous, v(0)=0, v'(0)=1
[~,  U_v] = rk4_bvp(rhs_hom,  0, L, [0;  1], h);

% Linear combination
s_dir = (TB - U_u(end,1)) / U_v(end,1);
T_shoot_dir = U_u(:,1) + s_dir * U_v(:,1);

%% aii
s_ins = -U_u(end,2) / U_v(end,2);
T_shoot_ins = U_u(:,1) + s_ins * U_v(:,1);

%% Part b
N  = 21;
dx = L/(N-1);
xj = linspace(0, L, N)';
aj = a(xj);  bj = b(xj);  fj = f(xj);

%Matrix
A  = zeros(N, N);
rhs_fd = zeros(N, 1);

%BCs
A(1,1) = 1;  rhs_fd(1) = TA;
A(N,N) = 1;  rhs_fd(N) = TB;

% Interior points
for j = 2:N-1
    Lj = 1/dx^2 - aj(j)/(2*dx);
    Dj = -2/dx^2 + bj(j);
    Rj = 1/dx^2 + aj(j)/(2*dx);
    A(j, j-1) = Lj;
    A(j, j)   = Dj;
    A(j, j+1) = Rj;
    rhs_fd(j) = fj(j);
end

T_fd = A \ rhs_fd;

%% Plots
figure;
plot(xg, T_shoot_dir, 'b-', 'DisplayName', 'Shooting Method'); hold on;
plot(xj, T_fd, 'ro',  'MarkerSize', 6, 'DisplayName', 'Finite Difference');
xlabel('x'); ylabel('T(x)');
title('Temperature distribution: T(0)=5, T(2)=4');
legend; grid on;

figure;
plot(xg, T_shoot_ins, 'b-');
xlabel('x'); ylabel('T(x)');
title(sprintf('Insulated end (dT/dx=0 at x=2): T(2) = %.4f', T_shoot_ins(end)));
grid on;
fprintf('Part (a)(iii): Temperature at insulated end T(L) = %.6f\n', T_shoot_ins(end));
