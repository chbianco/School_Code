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
% RHS
function dF = blasius_rhs(eta, F)
    f1=F(1); f2=F(2); f3=F(3);
    dF = [-0.5*f1*f3; f1; f2];         
end

%RK4 
function [eta, F] = rk4_blasius(f1_0, eta_end, deta)
    eta = (0:deta:eta_end)';
    N   = length(eta);
    F   = zeros(N, 3);
    F(1,:) = [f1_0, 0, 0];   
    for n = 1:N-1
        k1 = blasius_rhs(eta(n),          F(n,:)');
        k2 = blasius_rhs(eta(n)+deta/2,   F(n,:)' + deta/2*k1);
        k3 = blasius_rhs(eta(n)+deta/2,   F(n,:)' + deta/2*k2);
        k4 = blasius_rhs(eta(n)+deta,     F(n,:)' + deta*k3);
        F(n+1,:) = F(n,:) + (deta/6)*(k1+2*k2+2*k3+k4)';
    end
end

%% Solving
eta_inf = 10;
deta    = 0.01;

%Guesses
s0 = 1.0; [~, F0] = rk4_blasius(s0, eta_inf, deta); r0 = F0(end,2) - 1;
s1 = 0.5; [~, F1] = rk4_blasius(s1, eta_inf, deta); r1 = F1(end,2) - 1;

for iter = 1:20
    s2 = s1 - r1*(s1-s0)/(r1-r0);
    [~, F2] = rk4_blasius(s2, eta_inf, deta);
    r2 = F2(end,2) - 1;
    fprintf('iter %2d: f''''(0) = %.10f, f''(inf)-1 = %.2e\n', iter, s2, r2);
    if abs(r2) < 1e-10, break; end
    s0=s1; r0=r1; s1=s2; r1=r2;
end

fprintf('\nConverged: f''''(0) = %.6f\n', s2);

%Solution
[eta, F] = rk4_blasius(s2, eta_inf, deta);
f2 = F(:,2);   
f3 = F(:,3);

%% Plots
U = 1;
nu = 1.5e-5;
Rex_vals = [100, 500, 1000];
colors   = {'b','r','g'};

figure('Name','u/U profile');
hold on; grid on;
for k = 1:3
    Rex = Rex_vals(k);
    x = Rex * nu / U;  
    y = eta * x / sqrt(Rex);
    plot(f2, y, colors{k}, 'DisplayName', sprintf('$Re_x = %d$', Rex));
end
xlabel('u/U'); ylabel('y (m)');
title('Horizontal velocity profile'); legend; 

figure('Name','v/U profile');
hold on; grid on;
for k = 1:3
    Rex = Rex_vals(k);
    x = Rex * nu / U;
    y = eta * x / sqrt(Rex);
    vU = (eta .* f2 - f3) / (2*sqrt(Rex));
    plot(vU, y, colors{k}, 'DisplayName', sprintf('$Re_x = %d$', Rex));
end
xlabel('v/U'); ylabel('y (m)');
title('Vertical velocity profile'); legend;