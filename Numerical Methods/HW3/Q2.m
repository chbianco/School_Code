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


%% Parameters
alpha = 0;
c     = 1; 
L     = 1.0;
t_end = 1.0;

%% Grid
N  = 100;
x  = linspace(0, L, N+1)';
dx = x(2) - x(1);

%% Initial condition
u0 = exp(-200*(x - 0.25).^2);
dt_max_RK4 = 2*sqrt(2) * dx / c;

%% FE
dt_FE   = 0.9 * dt_max_RK4;
nsteps_FE = round(t_end / dt_FE);
u_FE    = u0;
t_FE    = 0;
snap_FE = zeros(N+1, 3);
snap_FE(:,1) = u0;
snap_idx_FE  = 2;
t_snaps      = [0.5, 1.0];
diverge_time = NaN;

for step = 1:nsteps_FE
    t_FE = step * dt_FE;
    k1   = rhs_central(u_FE, c, alpha, dx, N);
    u_FE = u_FE + dt_FE * k1;

    if any(~isfinite(u_FE)) || max(abs(u_FE)) > 1e6
        diverge_time = t_FE;
        fprintf('  -> DIVERGED at t = %.4f (step %d/%d)\n', ...
                t_FE, step, nsteps_FE);
        break;
    end

    for s = 1:2
        if abs(t_FE - t_snaps(s)) < 0.5*dt_FE && snap_idx_FE <= 3
            snap_FE(:, snap_idx_FE) = u_FE;
            snap_idx_FE = snap_idx_FE + 1;
        end
    end
end
if isnan(diverge_time)
    snap_FE(:,3) = u_FE;
    fprintf('  -> Completed without divergence');
end

%% RK4
dt_RK4   = 0.9 * dt_max_RK4;
nsteps_RK4 = round(t_end / dt_RK4);
u_RK4    = u0;
snap_RK4 = zeros(N+1, 3);
snap_RK4(:,1) = u0;
snap_idx_RK4  = 2;

for step = 1:nsteps_RK4
    t_RK4 = step * dt_RK4;
    k1 = rhs_central(u_RK4, c, alpha, dx, N);
    k2 = rhs_central(u_RK4 + dt_RK4/2*k1, c, alpha, dx, N);
    k3 = rhs_central(u_RK4 + dt_RK4/2*k2, c, alpha, dx, N);
    k4 = rhs_central(u_RK4 + dt_RK4*k3,   c, alpha, dx, N);
    u_RK4 = u_RK4 + (dt_RK4/6)*(k1 + 2*k2 + 2*k3 + k4);

    for s = 1:2
        if abs(t_RK4 - t_snaps(s)) < 0.5*dt_RK4 && snap_idx_RK4 <= 3
            snap_RK4(:, snap_idx_RK4) = u_RK4;
            snap_idx_RK4 = snap_idx_RK4 + 1;
        end
    end
end
snap_RK4(:,3) = u_RK4;

%% Plots
t_labels = {'t = 0', 't = 0.5', 't = 1'};
t_plot   = [0, 0.5, 1.0];
clr = [0.13 0.47 0.71;   % blue
       0.84 0.15 0.16;   % red
       0.17 0.63 0.17];  % green

figure('Position', [80 80 1200 620]);

for k = 1:3
    subplot(2, 3, k);
    plot(x, snap_RK4(:,k), '-', 'Color', clr(k,:), 'LineWidth', 2.2);
    xlabel('x', 'FontSize', 11);
    ylabel('u', 'FontSize', 11);
    title(sprintf('RK4 — %s', t_labels{k}), 'FontSize', 11);
    xlim([0 1]); ylim([-0.15 1.15]);
    grid on; box on;
end

for k = 1:3
    subplot(2, 3, 3+k);
    if ~isnan(diverge_time) && t_plot(k) >= diverge_time
        % Solution had already diverged at this snapshot time — plot nothing
        text(0.5, 0.5, sprintf('Diverged\nbefore t = %.1f', t_plot(k)), ...
            'HorizontalAlignment', 'center', 'FontSize', 11, ...
            'Color', [0.7 0 0], 'FontWeight', 'bold');
    else
        plot(x, snap_FE(:,k), '-', 'Color', clr(k,:), 'LineWidth', 2.2);
    end
    xlabel('x', 'FontSize', 11);
    ylabel('u', 'FontSize', 11);
    title(sprintf('Fwd Euler — %s', t_labels{k}), 'FontSize', 11);
    xlim([0 1]); ylim([-0.15 1.15]);
    grid on; box on;
end

%% RHS Operator 
function dudt = rhs_central(u, c, alpha, dx, N)
    dudt = zeros(N+1, 1);
    dudt(2:N) = -c * (u(3:N+1) - u(1:N-1)) / (2*dx) ...
               + alpha * (u(3:N+1) - 2*u(2:N) + u(1:N-1)) / dx^2;
    dudt(N+1) = -c * (u(N+1) - u(N)) / dx;
end