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
k1 = 2e3; k2 = 1e-3; k3 = 10;
u0 = [1; 5e-5; 0; 0];   % [CS, CE, CES, CP]

function dudt = rhs(u, k1, k2, k3)
    CS=u(1); CE=u(2); CES=u(3);
    dudt = [-k1*CS*CE + k2*CES; -k1*CS*CE + (k2+k3)*CES; k1*CS*CE - (k2+k3)*CES; k3*CES];
end

function J = jac(u, k1, k2, k3)
    CS=u(1); CE=u(2);
    J = [-k1*CE,  -k1*CS,   k2,      0;
         -k1*CE,  -k1*CS,   k2+k3,   0;
          k1*CE,   k1*CS,  -(k2+k3), 0;
          0,       0,        k3,      0];
end

%RK4
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
        U(n+1,:) = U(n,:) + (h/6)*(k1+2*k2+2*k3+k4)';
    end
end

%Linearlized trapezoidal
function [t, U] = lin_trap_solve(f, jac_f, t0, tf, u0, h)
    t = (t0:h:tf)';
    N = length(t);
    m = numel(u0);
    U = zeros(N, m);
    U(1,:) = u0(:)';
    for n = 1:N-1
        un = U(n,:)';
        fn = f(un);
        Jn = jac_f(un);
        A  = eye(m) - (h/2)*Jn;
        du = A \ (h * fn);
        U(n+1,:) = (un + du)';
    end
end

f = @(u) rhs(u, k1, k2, k3);
jf = @(u) jac(u, k1, k2, k3);

f_ode = @(t, u) rhs(u, k1, k2, k3); 

%% Part a
%ai
tic;
[t_rk, U_rk] = rk4_solve(f, 0, 1, u0, 1e-5);
t_rk4 = toc;
fprintf('RK4 time: %.3f s, steps: %d\n', t_rk4, length(t_rk));

%aii
tic;
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t_s, U_s] = ode23s(f_ode, [0 1], u0, opts);
t_stiff = toc;
fprintf('ode23s time: %.3f s, steps: %d\n', t_stiff, length(t_s));

%Plots
figure;
vars   = {'$C_S$', '$C_E$', '$C_{ES}$', '$C_P$'};
for i = 1:4
    subplot(2,2,i);
    loglog(t_rk(2:end), U_rk(2:end,i), 'b-',  'LineWidth', 1.2); hold on;
    loglog(t_s(2:end),  U_s(2:end,i),  'r--', 'LineWidth', 1.2);
    xlabel('t'); ylabel(vars{i});
    legend('RK4','ode23s','Location','best');
    grid on;
    title(vars{i});
end
sgtitle('Enzyme Concentrations');



%% Part b
tic;
[t_lt, U_lt] = lin_trap_solve(f, jf, 0, 1, u0, 1e-5);
t_lintrap = toc;
fprintf('Trapezoidal time: %.3f s, steps: %d\n', t_lintrap, length(t_lt));

figure;
vars   = {'$C_S$', '$C_E$', '$C_{ES}$', '$C_P$'};
for i = 1:4
    subplot(2,2,i);
    loglog(t_lt(2:end), U_lt(2:end,i), 'b-'); hold on;
    xlabel('t'); ylabel(vars{i});
    grid on;
    title(vars{i});
end
sgtitle('Enzyme Concentrations (Linearized Trapezoidal)');