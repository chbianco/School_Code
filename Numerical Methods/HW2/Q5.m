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

%% Solving
h = 0.2;
tspan = 0:h:100;
N = length(tspan);

exact = @(y0, t) -log(exp(-y0) + exp(-t) - 1);
f = @(y, t)  exp(y - t);

y0 = -1;

%implicit Euler 
y_fi = zeros(1,N); y_fi(1) = y0;
for n = 1:N-1
    tn1 = tspan(n+1);
    y_iter = y_fi(n);
    for k = 1:50
        F  = y_iter - y_fi(n) - h*exp(y_iter - tn1);
        dF = 1 - h*exp(y_iter - tn1);
        dy = -F/dF;
        y_iter = y_iter + dy;
        if abs(dy) < 1e-12, break; end
    end
    y_fi(n+1) = y_iter;
end

% Linearized implicit Euler 
y_li = zeros(1,N); y_li(1) = y0;
for n = 1:N-1
    fn = exp(y_li(n) - tspan(n));
    y_li(n+1) = y_li(n) + h*fn*(1 - h) / (1 - h*fn);
end

y_ex = exact(y0, tspan);

%% Plotting
figure;
subplot(1,2,1);
plot(tspan, y_ex,  'k-',  'LineWidth', 2, 'DisplayName', 'Exact'); hold on;
plot(tspan, y_fi,  'b--', 'LineWidth', 1.5, 'DisplayName', 'Fully Implicit');
plot(tspan, y_li,  'r:',  'LineWidth', 1.5, 'DisplayName', 'Linearized Implicit');
xlabel('t'); ylabel('y'); title('$y_0 = -10^{-5}$');
legend('Location','best'); grid on;

subplot(1,2,2);
semilogy(tspan, abs(y_fi - y_ex), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Fully Implicit'); hold on;
semilogy(tspan, abs(y_li - y_ex), 'r:',  'LineWidth', 1.5, 'DisplayName', 'Linearized');
xlabel('t'); ylabel('$|error|$'); title('Errors, $y_0 = -10^{-5}$');
legend('Location','best'); grid on;
sgtitle('Problem 4.12(c)');