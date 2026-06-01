%% Preamble
close all; clc; clear all;
set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);

%% Fourier Coefficients via FFT
 
f = @(x) cos(2*x) + (1/2)*cos(4*x) + (1/6)*cos(12*x);
 
N_values = [8, 16, 32, 64];
 
fprintf('%-5s  %-6s  %-12s  %-12s\n', 'N', 'k', 'Re(f_hat_k)', 'Im(f_hat_k)');
fprintf('%s\n', repmat('-',1,45));
 
for N = N_values
    x = (0:N-1) * (2*pi/N);
    fx = f(x);
 
    F = fft(fx) / N;
 
    fprintf('\nN = %d\n', N);
    for k = 0:N-1
        fprintf('  k=%-3d  %+.6f  %+.6fi\n', k, real(F(k+1)), imag(F(k+1)));
    end
end

%% N = 32 verification
x = (0:N-1) * (2*pi/N);
fx = f(x);
 
F = fft(fx) / N;
fx_reconstructed = N * ifft(F);
 
% Verify reconstruction
max_error = max(abs(fx_reconstructed - fx));
fprintf('\n--- N = 32 Verification ---\n');
fprintf('Max |f_reconstructed - f_original| = %.2e\n', max_error);
if max_error < 1e-10
    fprintf('PASSED: Reconstructed signal is identical to original (within machine precision).\n');
else
    fprintf('FAILED: Max error exceeds threshold.\n');
end

%% Plots
k = 0:N-1;
abs_F = abs(F);
 
figure('Name','Fourier Coefficients (N=32)', 'Position',[100 100 800 500]);
stem(k, abs_F, 'filled', 'Color', [0.1 0.4 0.8], 'MarkerSize', 7, 'LineWidth', 1.5);
hold on;
 
% Annotate expected non-zero frequencies
expected_k  = [2, N-2, 4, N-4, 12, N-12];
expected_amp = [0.5, 0.5, 0.25, 0.25, 1/12, 1/12];
stem(expected_k, expected_amp, 'r^', 'MarkerSize', 9, 'LineWidth', 1.5, 'DisplayName','Expected');
 
xlabel('k', 'FontSize', 13);
ylabel('$|\hat{f}_k|$',          'FontSize', 13);
legend('Computed $|{\hat f}_k|$', 'Expected peaks', 'Location','northeast');
grid on;
xlim([-1 N]);
xticks(0:2:N-1);