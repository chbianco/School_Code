%% Preamble
close all; clc; clear all;
set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);

%% Solving
cases(1).fun    = @(x) sin(3*x) + 3*cos(6*x);
cases(1).dfun  = @(x) 3*cos(3*x) - 18*sin(6*x);
cases(1).name  = 'sin(3x)+3cos(6x)';
 
cases(2).fun   = @(x) 6*x - x.^2;
cases(2).dfun  = @(x) 6 - 2*x;
cases(2).name  = '6x - x^2';
 
N_values = [16, 32];
 
figure('Position',[50 50 1400 700]);
pnum = 1;
 
for c = 1:2
    fun  = cases(c).fun;
    dfun = cases(c).dfun;
 
    for N = N_values
        x  = (0:N-1)' * (2*pi/N);   % column vector
        h  = 2*pi/N;
        fx    = fun(x);
        exact = dfun(x);
 
        %--- FFT derivative ---
        % Multiply FFT by i*k, then IFFT
        F  = fft(fx);
        k  = [0:N/2-1, 0, -N/2+1:-1]';  % wavenumbers (zero out Nyquist mode)
        dF = (1i * k) .* F;
        df_fft = real(ifft(dF));
 
        %--- 2nd-order central finite differences (periodic) ---
        ip1 = mod((0:N-1)+1, N) + 1;  % index i+1 (periodic)
        im1 = mod((0:N-1)-1, N) + 1;  % index i-1 (periodic)
        df_fd = (fx(ip1) - fx(im1)) / (2*h);
 
        err_fft = abs(df_fft - exact);
        err_fd  = abs(df_fd  - exact);
 
        fprintf('Case %d (%s)  N=%2d:  FFT max err = %.3e,  FD2 max err = %.3e\n', ...
            c, cases(c).name, N, max(err_fft), max(err_fd));
 
        %--- Derivative plot ---
        subplot(2,4,pnum);
        plot(x, exact,  'g-',  'LineWidth',2,   'DisplayName','Exact');  hold on;
        plot(x, df_fft, 'b--', 'LineWidth',1.5, 'DisplayName','FFT');
        plot(x, df_fd,  'r:',  'LineWidth',1.5, 'DisplayName','FD2');
        title(sprintf('Deriv: case %d, N=%d', c, N));
        xlabel('x'); 
        ylabel("f'(x)")
        legend('Location','best'); grid on; hold off;
 
        %--- Error plot ---
        subplot(2,4,pnum+1);  % error plot in next column
        semilogy(x, err_fft+1e-17, 'b-', 'LineWidth',1.8, ...
            'DisplayName', sprintf('FFT max=%.2e',max(err_fft)));
        hold on;
        semilogy(x, err_fd+1e-17,  'r-', 'LineWidth',1.8, ...
            'DisplayName', sprintf('FD2 max=%.2e',max(err_fd)));
        title(sprintf('Error: case %d, N=%d', c, N));
        xlabel('x'); ylabel('|error|'); legend('Location','best'); grid on; hold off;
 
        pnum = pnum + 2;
    end
end
sgtitle('FFT vs 2nd-Order FD Differentiation');