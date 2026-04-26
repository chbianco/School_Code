%% Preamble
close all; clc; clear all;
set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);

%% Setup
nu   = 0.01;
CFL  = 0.4;
T_out = [1.0, 2.0, 3.0];
N_list = [32, 64, 128, 256];

%% Helper: nonlinear term with anti-aliasing

function Nhat = nonlinear(u_hat, N, k, anti_alias)
    if anti_alias
        M = floor(3*N/2);
        u_hat_pad = zeros(M,1);
        u_hat_pad(1:N/2)       = u_hat(1:N/2);
        u_hat_pad(M-N/2+1:M)   = u_hat(N/2+1:N);
        u_pad = real(ifft(u_hat_pad)) * (M/N);
        U2_M  = fft(u_pad.^2 / 2);
        U2 = zeros(N,1);
        U2(1:N/2)     = U2_M(1:N/2);
        U2(N/2+1:N)   = U2_M(M-N/2+1:M);
        Nhat = -1i * k .* U2;
    else
        u = real(ifft(u_hat));
        Nhat = -1i * k .* fft(u.^2 / 2);
    end
end

%% Solver
function results = solve_burgers(N, T_out, nu, CFL, anti_alias)
    x = linspace(-pi, pi, N+1)'; x = x(1:end-1);
    h = 2*pi/N;
    k = [0:N/2-1, 0, -N/2+1:-1]';
 
    u = -sin(x);
    u_hat = fft(u);
    L = -nu * k.^2;
 
    results = struct();
    t  = 0;
    dt0 = 0.001;
 
    %Forward Euler step to get NL_{n-1}
    NL_old = nonlinear(u_hat, N, k, anti_alias);
    u_hat  = ((1 + L*dt0/2) .* u_hat + dt0*NL_old) ./ (1 - L*dt0/2);
    t      = dt0;
    NL_old2 = NL_old;
    NL_old  = nonlinear(u_hat, N, k, anti_alias);
 
    T_max = max(T_out);
 
    while t < T_max - 1e-12
        umax = max(abs(real(ifft(u_hat)))) + 1e-10;
        dt   = min([CFL*h/umax, 0.005, T_max-t]);
        for T = sort(T_out)
            if t < T && T <= t+dt+1e-12
                dt = T - t; break;
            end
        end
 
        NL_rhs  = 1.5*NL_old - 0.5*NL_old2;
        lhs     = 1 - L*dt/2;
        u_hat   = ((1 + L*dt/2).*u_hat + dt*NL_rhs) ./ lhs;
 
        NL_old2 = NL_old;
        NL_old  = nonlinear(u_hat, N, k, anti_alias);
        t = t + dt;
 
        for j = 1:length(T_out)
            T = T_out(j);
            fname = sprintf('t%d', round(T*10));
            if abs(t-T) < 1e-7 && ~isfield(results, fname)
                results.(fname).u = real(ifft(u_hat));
                results.(fname).x = x;
                label = 'AA'; if ~anti_alias; label = 'no-AA'; end
                fprintf('  N=%3d %s: t=%.4f  max|u|=%.4f\n', N, label, t, max(abs(real(ifft(u_hat)))));
            end
        end
    end
end

%% part a
results_noaa = struct();
for N = N_list
    res = solve_burgers(N, T_out, nu, CFL, false);
    results_noaa.(sprintf('N%d',N)) = res;
end
 
figure('Position',[50 50 1400 450]);
cols = {'r','y','b','g'};
for ti = 1:3
    T = T_out(ti);
    subplot(1,3,ti);
    for ni = 1:length(N_list)
        N = N_list(ni);
        res = results_noaa.(sprintf('N%d',N));
        fname = sprintf('t%d', round(T*10));
        plot(res.(fname).x, res.(fname).u, cols{ni}, 'LineWidth',1.8, ...
             'DisplayName', sprintf('N=%d',N)); hold on;
    end
    title(sprintf('t = %.0f', T)); xlabel('x'); ylabel('u(x,t)');
    legend('Location','best'); grid on; hold off;
end
sgtitle('Burgers  $\nu= 0.01$,  $u_0=-\sin x$  (no anti-aliasing, AB2+CN)');

%% part b
res_aa_32  = solve_burgers(32,  [2.0], nu, CFL, true);
res_aa_256 = solve_burgers(256, [2.0], nu, CFL, true);
 
figure('Position',[50 550 1100 450]);
for pi_ = 1:2
    N = [32, 256]; N = N(pi_);
    subplot(1,2,pi_);
    res_noaa_plot = results_noaa.(sprintf('N%d',N));
    if pi_ == 1; res_aa_plot = res_aa_32; else; res_aa_plot = res_aa_256; end
    plot(res_aa_plot.t20.x,   res_aa_plot.t20.u,   'b-',  'LineWidth',2,   'DisplayName','Anti-aliased (3/2 rule)');
    hold on;
    plot(res_noaa_plot.t20.x, res_noaa_plot.t20.u, 'r--', 'LineWidth',1.5, 'DisplayName','No anti-aliasing');
    title(sprintf('t=2,  N=%d', N)); xlabel('x'); ylabel('u(x,t)');
    legend('Location','best'); grid on; hold off;
end
sgtitle('Burgers at t=2: Anti-aliased vs No Anti-aliasing');