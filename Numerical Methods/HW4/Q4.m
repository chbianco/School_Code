%% Preamble
close all; clc; clear all;
set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);

%% Grid and fuctions
N  = 32;
x  = (0:N-1)' * (2*pi/N);
 
y1 = sin(3*x) + 0.1*sin(16*x);
y2 = cos(3*x) + 0.1*cos(16*x);

%% Aliased
prod_naive = y1 .* y2;
P_naive = fft(prod_naive) / N;

k_vals = [0:N/2-1, -N/2:-1];
for i = 1:N
    if abs(P_naive(i)) > 1e-10
        amp = 2*abs(P_naive(i));
        fprintf('  k=%3d:  amplitude = %.4f\n', k_vals(i), amp);
    end
end

%% Anti-aliased
M = 3*N/2;
xM = (0:M-1)' * (2*pi/M);
 
y1M = sin(3*xM) + 0.1*sin(16*xM);
y2M = cos(3*xM) + 0.1*cos(16*xM);
 
prod_M = y1M .* y2M;
 
% FFT on M grid, then truncate to N modes
PM = fft(prod_M) / M;
 
P_aa = zeros(N, 1);
P_aa(1:N/2)   = PM(1:N/2);        % k = 0 .. N/2-1
P_aa(N/2+1:N) = PM(M-N/2+1:M);   % k = -N/2 .. -1
 
for i = 1:N
    if abs(P_aa(i)) > 1e-10
        amp = 2*abs(P_aa(i));
        fprintf('  k=%3d:  amplitude = %.4f\n', k_vals(i), amp);
    end
end
 
prod_aa = real(ifft(P_aa * N));

%% Fine grid for real comparison
Nf = 1000;
xf = linspace(0, 2*pi, Nf+1); xf = xf(1:end-1)';
prod_true = (sin(3*xf)+0.1*sin(16*xf)) .* (cos(3*xf)+0.1*cos(16*xf));

%% Plots
figure 

% Physical space
subplot(1,3,1);
plot(xf, prod_true, 'g-', 'LineWidth',2, 'DisplayName','True product'); hold on;
plot(x,  prod_naive,'r--o','LineWidth',1.5,'MarkerSize',5,'DisplayName','Naive (aliased)');
plot(x,  prod_aa,   'b:s', 'LineWidth',1.5,'MarkerSize',5,'DisplayName','Anti-aliased (3/2)');
xlabel('x'); ylabel('$y_1 y_2$');
title('Physical space');
legend('Location','best'); grid on;
 
% Spectra
subplot(1,3,2);
k_plot = k_vals;
stem(k_plot, 2*abs(P_naive), 'r', 'filled', 'DisplayName','Naive'); hold on;
stem(k_plot+0.4, 2*abs(P_aa), 'b', 'filled', 'DisplayName','Anti-aliased');
xlabel('k'); ylabel('$2|c_k|$');
title('Fourier spectra');
legend('Location','best'); grid on;
xlim([-17 17]);
 
% Product on M=48 grid vs true
subplot(1,3,3);
plot(xf, prod_true, 'g-', 'LineWidth',2, 'DisplayName','True product'); hold on;
plot(xM, prod_M, 'ko-', 'LineWidth',1.5,'MarkerSize',4,'DisplayName',sprintf('M=%d grid',M));
xlabel('x'); ylabel('$y_1 y_2$');
title(sprintf('Product on 3/2-padded M=%d grid', M));
legend('Location','best'); grid on;
 