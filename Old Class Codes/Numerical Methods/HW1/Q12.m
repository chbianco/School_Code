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
N = 200;
x = linspace(0, pi, N);
h = x(2) - x(1);

w = ones(N,1)*h;
w(1)   = h/2;
w(end) = h/2;

K = zeros(N,N);

for i = 1:N
    for j = 1:N
        K(i,j) = 3*(0.5*sin(3*x(i)) - x(i)^2 * x(j));
    end
end

A = K .* w';      
M = eye(N) - A;

f = pi * x.^2;

phi_num = M \ f';

phi_exact = sin(3*x);

%% Plot
figure; hold on;
grid on
plot(x,phi_num,'b');
plot(x,phi_exact,'r--');
legend({'Numerical','Exact'}, 'Location','southeast');
xlabel('x'); ylabel('$\phi(x)$');
