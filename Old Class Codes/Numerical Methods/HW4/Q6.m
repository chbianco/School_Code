%% Preamble
close all; clc; clear all;
set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);

%% Modified wavenumber plots 

kh = linspace(-pi, pi, 1000);

kp_CDS2 = sin(kh);
kp_UDS1 = (1 - exp(-1j*kh)) / 1j;
kp_QUICK2 = (3 - 4*exp(-1j*kh) + exp(-2j*kh)) / (2j);
kp_FOURIER = kh;

tol = 1e-12;
safe_kh = kh;
safe_kh(abs(kh) < tol) = tol;

c_ratio_CDS2   = real(kp_CDS2)./ safe_kh;
c_ratio_UDS1   = real(kp_UDS1)./ safe_kh;
c_ratio_QUICK2 = real(kp_QUICK2)./ safe_kh;
c_ratio_FOURIER = ones(size(kh));

ci_ratio_CDS2   = imag(kp_CDS2)./ safe_kh;
ci_ratio_UDS1   = imag(kp_UDS1)./ safe_kh;
ci_ratio_QUICK2 = imag(kp_QUICK2)./ safe_kh;
ci_ratio_FOURIER = zeros(size(kh));

%Re(c'/c) vs k*dx
figure
plot(kh, c_ratio_CDS2,    'b-',  'LineWidth', 2,   'DisplayName', 'CDS2');   hold on;
plot(kh, c_ratio_UDS1,    'r--', 'LineWidth', 2,   'DisplayName', 'UDS1');
plot(kh, c_ratio_QUICK2,  'm-.', 'LineWidth', 2,   'DisplayName', 'QUICK2');
plot(kh, c_ratio_FOURIER, 'g:',  'LineWidth', 2.5, 'DisplayName', 'FOURIER (exact)');
yline(1, 'k:', 'LineWidth', 1.2, 'DisplayName', 'Exact = 1');
yline(0, 'k-', 'LineWidth', 0.5);

xlabel('$k\Delta x$',      'Interpreter','latex', 'FontSize', 14);
ylabel('$\mathrm{Re}(c''/c)$', 'Interpreter','latex', 'FontSize', 14);
title('Modified phase speed: $\mathrm{Re}(c''/c)$ vs $k\Delta x$', ...
      'Interpreter','latex', 'FontSize', 13);
legend('Location','southwest', 'FontSize', 11);
grid on;
xlim([-pi pi]);
ylim([-0.2 1.2]);
xticks([-pi, -pi/2, 0, pi/2, pi]);
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});

%Im(c'/c) vs k*dx
figure
plot(kh, ci_ratio_CDS2,'b-', 'DisplayName', 'CDS2');   hold on;
plot(kh, ci_ratio_UDS1, 'r--', 'DisplayName', 'UDS1');
plot(kh, ci_ratio_QUICK2,'m-.','DisplayName', 'QUICK2');
plot(kh, ci_ratio_FOURIER, 'g:',  'LineWidth', 2.5, 'DisplayName', 'FOURIER (exact)');
yline(0, 'k-', 'LineWidth', 0.8);

xlabel('$k\Delta x$',      'Interpreter','latex', 'FontSize', 14);
ylabel('$\mathrm{Im}(c''/c)$', 'Interpreter','latex', 'FontSize', 14);
title('Dissipation: $\mathrm{Im}(c''/c)$ vs $k\Delta x$', ...
      'Interpreter','latex', 'FontSize', 13);
legend('Location','southwest', 'FontSize', 11);
grid on;
xlim([-pi pi]);
xticks([-pi, -pi/2, 0, pi/2, pi]);
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});


%% Setup
c = 2*pi;
N = 200;
x = linspace(-pi, pi, N+1)'; x = x(1:end-1);
dx = 2*pi/N;
T = 2*pi/c;
t_end = 10*T;
CFL = 0.1;
dt  = CFL*dx/c;
Nsteps = ceil(t_end/dt);
dt = t_end/Nsteps;
u0 = exp(-4*x.^2);

%Exact soln
exact = @(t) exp(-4*( mod(x - c*t + pi, 2*pi) - pi ).^2);

%% Helpers
% Spatial derivatives
function du = cds2(u, dx)
    du = (circshift(u,-1) - circshift(u,1)) / (2*dx);
end
 
function du = uds1(u, dx)
    du = (u - circshift(u,1)) / dx;
end
 
function du = quick2(u, dx)
    du = (3*u - 4*circshift(u,1) + circshift(u,2)) / (2*dx);
end
 
function du = fourier_deriv(u)
    N = length(u);
    k = [0:N/2-1, 0, -N/2+1:-1]';
    du = real(ifft(1j*k.*fft(u)));
end

%RHS
rhs_CDS2    = @(u) -c * cds2(u, dx);
rhs_UDS1    = @(u) -c * uds1(u, dx);
rhs_QUICK2  = @(u) -c * quick2(u, dx);
rhs_FOURIER = @(u) -c * fourier_deriv(u);

% RK4
function u_new = rk4(u, dt, rhs)
    k1 = rhs(u);
    k2 = rhs(u + dt/2*k1);
    k3 = rhs(u + dt/2*k2);
    k4 = rhs(u + dt*k3);
    u_new = u + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

%% Solving
% Time integration
u_CDS2    = u0;
u_UDS1    = u0;
u_QUICK2  = u0;
u_FOURIER = u0;
 
E_CDS2    = zeros(Nsteps+1,1); E_CDS2(1)    = sum(u0.^2)*dx;
E_UDS1    = zeros(Nsteps+1,1); E_UDS1(1)    = sum(u0.^2)*dx;
E_QUICK2  = zeros(Nsteps+1,1); E_QUICK2(1)  = sum(u0.^2)*dx;
E_FOURIER = zeros(Nsteps+1,1); E_FOURIER(1) = sum(u0.^2)*dx;
t_hist    = linspace(0, t_end, Nsteps+1)';
 
fprintf('Integrating: N=%d, dt=%.5f, Nsteps=%d\n', N, dt, Nsteps);
for step = 1:Nsteps
    u_CDS2    = rk4(u_CDS2,    dt, rhs_CDS2);
    u_UDS1    = rk4(u_UDS1,    dt, rhs_UDS1);
    u_QUICK2  = rk4(u_QUICK2,  dt, rhs_QUICK2);
    u_FOURIER = rk4(u_FOURIER, dt, rhs_FOURIER);
    E_CDS2(step+1) = sum(u_CDS2.^2)*dx;
    E_UDS1(step+1)  = sum(u_UDS1.^2)*dx;
    E_QUICK2(step+1)  = sum(u_QUICK2.^2)*dx;
    E_FOURIER(step+1) = sum(u_FOURIER.^2)*dx;
end
 
u_ex = exact(t_end);
fprintf('\nMax errors at t=%.1f:\n', t_end);
fprintf('  CDS2:    %.3e\n', max(abs(u_CDS2    - u_ex)));
fprintf('  UDS1:    %.3e\n', max(abs(u_UDS1    - u_ex)));
fprintf('  QUICK2:  %.3e\n', max(abs(u_QUICK2  - u_ex)));
fprintf('  FOURIER: %.3e\n', max(abs(u_FOURIER - u_ex)));

%% Figures
% h
figure;
plot(x, u_ex,      'k-',  'LineWidth',2.5, 'DisplayName','Exact'); hold on;
plot(x, u_CDS2,    'b-',  'LineWidth',1.8, 'DisplayName','CDS2');
plot(x, u_UDS1,    'r--', 'LineWidth',1.8, 'DisplayName','UDS1');
plot(x, u_QUICK2,  'm-.', 'LineWidth',1.8, 'DisplayName','QUICK2');
plot(x, u_FOURIER, 'g:',  'LineWidth',2,   'DisplayName','FOURIER');
xlabel('x'); ylabel('u(x,t)');
legend('Location','best'); grid on; hold off;

% i
figure;
t_norm = t_hist / T;
plot(t_norm, E_CDS2,'b-',  'LineWidth',1.8, 'DisplayName','CDS2');   hold on;
plot(t_norm, E_UDS1,'r--', 'LineWidth',1.8, 'DisplayName','UDS1');
plot(t_norm, E_QUICK2,'m-.', 'LineWidth',1.8, 'DisplayName','QUICK2');
plot(t_norm, E_FOURIER,'g:', 'DisplayName','FOURIER');
yline(E_FOURIER(1), 'k:', 'LineWidth',1, 'DisplayName','Exact (constant)');
xlabel('t / T  (flow-through times)'); ylabel('$\int u^2 dx$');
title('Global energy vs time'); legend('Location','best'); grid on; hold off;

%% Part j
Nc  = 50;
dxc = 2*pi/Nc;
 
function C = build_C(Nc, dxc, scheme, c)
    C = zeros(Nc,Nc);
    for j = 1:Nc
        im1 = mod(j-2,Nc)+1; im2 = mod(j-3,Nc)+1; ip1 = mod(j,Nc)+1;
        if strcmp(scheme,'CDS2')
            C(j,ip1) =  1/(2*dxc); C(j,im1) = -1/(2*dxc);
        elseif strcmp(scheme,'UDS1')
            C(j,j)   =  1/dxc;     C(j,im1) = -1/dxc;
        elseif strcmp(scheme,'QUICK2')
            C(j,j)   =  3/(2*dxc); C(j,im1) = -4/(2*dxc); C(j,im2) = 1/(2*dxc);
        end
    end
    C = -c*C;
end
 
figure;
scheme_names = {'CDS2','UDS1','QUICK2'};
cols = {'b','r','m'};
for si = 1:3
    s = scheme_names{si};
    C = build_C(Nc, dxc, s, c);
    S = C + C';
    eigs_S = eig(S);
    skew_norm = max(max(abs(S)));
    fprintf('\n%s: max|C+C^T| = %.2e\n', s, skew_norm);
    fprintf('  Re(eig) range: [%.3f, %.3f]\n', min(real(eigs_S)), max(real(eigs_S)));
 
    subplot(1,3,si);
    plot(real(eigs_S), imag(eigs_S), [cols{si},'o'], 'MarkerSize',5, 'MarkerFaceColor',cols{si});
    hold on; xline(0,'k--','LineWidth',0.8); hold off;
    title(sprintf('%s: $\\mathrm{eig}(C+C^T)$,  $\\max|S|=%.1e$', s, skew_norm), ...
          'Interpreter','latex', 'FontSize', 11);
    xlabel('$\mathrm{Re}(\lambda)$', 'Interpreter','latex', 'FontSize', 12);
    ylabel('$\mathrm{Im}(\lambda)$', 'Interpreter','latex', 'FontSize', 12);
    grid on;
end
sgtitle('Eigenvalue spectra of $C + C^T$,  $N=50$', 'Interpreter','latex', 'FontSize', 13);