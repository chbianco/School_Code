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

%% Parameters
alpha = 0.01;
L = 1.0;
t_end = 1.0;
N = 200;

%% Grid
x = linspace(0, L, N+1)';
dx = x(2) - x(1);
xi = x(2:N);
Ni = N - 1;

%% Time step: satisfy both CFL and diffusion stability constraints
dt_diff = dx^2 / (2*alpha);
dt_adv  = dx   / 1.0;
dt = 0.4 * min(dt_diff, dt_adv);

%% Operators
e = ones(Ni, 1);
D1 = spdiags([-e, zeros(Ni,1), e], [-1, 0, 1], Ni, Ni) / (2*dx);
D2 = spdiags([e, -2*e, e], [-1, 0, 1], Ni, Ni) / dx^2;
I = speye(Ni);

%% Initial condition
u = sin(pi * xi);

%% Storage for snapshots
snap_times = [0, 0.25, 0.5, 0.75, 1.0];
Ns = length(snap_times);
snaps = zeros(N+1, Ns);
snaps(:,1) = [0; u; 0];   % t = 0
snap_idx = 2;
t = 0.0;
nsteps = round(t_end / dt);

%% Time loop
for step = 1:nsteps
    t = t + dt;
    L_op = alpha*D2 - spdiags(u, 0, Ni, Ni)*D1;
    % LHS and RHS matrices
    LHS = I - (dt/2) * L_op;
    RHS = I + (dt/2) * L_op;
    rhs = RHS * u;

    u = LHS \ rhs;

    % Save snapshots
    for k = snap_idx:Ns
        if abs(t - snap_times(k)) < 0.5*dt
            snaps(:,k) = [0; u; 0];
            snap_idx = k + 1;
            break;
        end
    end
end
snaps(:,Ns) = [0; u; 0];

%% Plot
clr = [0.13 0.47 0.71;
       0.17 0.63 0.17;
       0.84 0.15 0.16;
       0.84 0.50 0.05;
       0.50 0.18 0.56];
t_labels = {'$t = 0$', '$t = 0.25$', '$t = 0.5$', '$t = 0.75$', '$t = 1$'};
figure
hold on
for k = 1:Ns
    plot(x, snaps(:,k), '-', 'Color', clr(k,:));
end
hold off
xlabel('$x$');
ylabel('$u$');
xlim([0 1]);
ylim([-0.05 1.1]);
grid on;
legend(t_labels, 'Location', 'northwest');