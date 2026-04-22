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
M  = 20;
N  = 20;
dt = 0.05;
t_end = 1.0;

%% Grid
x  = linspace(-1, 1, M+1);
y  = linspace(-1, 1, N+1);
dx = x(2) - x(1);
dy = y(2) - y(1);
[X, Y] = meshgrid(x, y);
X = X'; Y = Y';

%% Source term and exact steady state
Q  = 2*(2 - X.^2 - Y.^2);
phi_ss = (X.^2 - 1) .* (Y.^2 - 1);

%% Initial condition
phi = zeros(M+1, N+1);

%% Tridiagonal coefficients
rx = dt / (2 * dx^2);
ry = dt / (2 * dy^2);

Mx = M - 1;  
Ny = N - 1;

% tridiag(a, b, c, f): a=main diagonal, b=sub-diagonal, c=super-diagonal
ax_main  = (1 + 2*rx) * ones(Mx, 1);
ax_sub   = -rx * ones(Mx, 1);
ax_super = -rx * ones(Mx, 1);

ay_main  = (1 + 2*ry) * ones(Ny, 1);
ay_sub   = -ry * ones(Ny, 1);
ay_super = -ry * ones(Ny, 1);

%% Storage for snapshots at t = 0, 0.25, 1.0
snap_times  = [0.0, 0.25, 1.0];
snapshots   = cell(length(snap_times), 1);
snap_labels = {'t = 0.0', 't = 0.25', 't = 1.0'};
snapshots{1} = phi;   % t = 0

%% Time integration
t = 0.0;
nsteps = round(t_end / dt);

for step = 1:nsteps
    t = t + dt;
    xi = zeros(M+1, N+1);
    xi(2:M, 2:N) = phi(2:M, 2:N) ...
        + ry * (phi(2:M, 3:N+1) - 2*phi(2:M, 2:N) + phi(2:M, 1:N-1));

    % Step 2: r = (I + rx*Ax) * xi + dt*q
    rhs = zeros(M+1, N+1);
    rhs(2:M, 2:N) = xi(2:M, 2:N) ...
        + rx * (xi(3:M+1, 2:N) - 2*xi(2:M, 2:N) + xi(1:M-1, 2:N)) ...
        + dt * Q(2:M, 2:N);

    eta = zeros(M+1, N+1);
    for j = 2:N
        d = rhs(2:M, j);
        eta(2:M, j) = tridiag(ax_main, ax_sub, ax_super, d);
    end

    phi_new = zeros(M+1, N+1);
    for i = 2:M
        d = eta(i, 2:N)';
        phi_new(i, 2:N) = tridiag(ay_main, ay_sub, ay_super, d)';
    end
    phi = phi_new;

    for k = 1:length(snap_times)
        if abs(t - snap_times(k)) < 0.5*dt
            snapshots{k} = phi;
        end
    end
end

snapshots{end} = phi;

%% Error vs exact steady state
err = max(max(abs(phi - phi_ss)));
fprintf('Max error at t=1 vs exact steady state: %.4e\n', err);
fprintf('(~%.2f%% relative error)\n', err / max(max(abs(phi_ss))) * 100);

%% Plotting
figure('Position', [100 100 1200 380]);
for k = 1:3
    subplot(1, 3, k);
    surf(X, Y, snapshots{k}, 'EdgeColor', 'none');
    xlabel('x'); ylabel('y'); zlabel('$\phi$');
    title(snap_labels{k}, 'FontSize', 13);
    xlim([-1 1]); ylim([-1 1]); zlim([0 1.05]);
    colormap(parula);
    view([-35 25]);
    box on; grid on;
end


%% Large dt convergence
fprintf('\nConvergence study with large time steps (dt = 1.0):\n');
dt_big   = 1.0;
rx2      = dt_big / (2*dx^2);
ry2      = dt_big / (2*dy^2);

ax2_main  = (1 + 2*rx2) * ones(Mx, 1);
ax2_sub   = -rx2 * ones(Mx, 1);
ax2_super = -rx2 * ones(Mx, 1);
ay2_main  = (1 + 2*ry2) * ones(Ny, 1);
ay2_sub   = -ry2 * ones(Ny, 1);
ay2_super = -ry2 * ones(Ny, 1);

phi2 = zeros(M+1, N+1);
for ns = 1:4
    xi2 = zeros(M+1, N+1);
    xi2(2:M,2:N) = phi2(2:M,2:N) ...
        + ry2*(phi2(2:M,3:N+1) - 2*phi2(2:M,2:N) + phi2(2:M,1:N-1));

    rhs2 = zeros(M+1, N+1);
    rhs2(2:M,2:N) = xi2(2:M,2:N) ...
        + rx2*(xi2(3:M+1,2:N) - 2*xi2(2:M,2:N) + xi2(1:M-1,2:N)) ...
        + dt_big * Q(2:M,2:N);

    eta2 = zeros(M+1, N+1);
    for j = 2:N
        eta2(2:M,j) = tridiag(ax2_main, ax2_sub, ax2_super, rhs2(2:M,j));
    end

    pn = zeros(M+1, N+1);
    for i = 2:M
        pn(i,2:N) = tridiag(ay2_main, ay2_sub, ay2_super, eta2(i,2:N)')';
    end
    phi2 = pn;

    e = max(max(abs(phi2 - phi_ss)));
    fprintf('  t = %d (step %d):  max error = %.4e\n', ns, ns, e);
end