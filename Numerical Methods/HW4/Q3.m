%% Preamble
close all; clc; clear all;
set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);

%% Part a
N = 128;
h = 2*pi/N;
x = linspace(-pi, pi, N+1); x = x(1:end-1); 
[X, Y] = meshgrid(x, x);
f = exp(-(X.^2 + Y.^2)) - 1/(4*pi);

figure('Name','(a) Source term');
surf(X, Y, f, 'EdgeColor','none'); colorbar;
xlabel('x'); ylabel('y'); zlabel('f');
title('Source term $f = e^{-(x^2+y^2)} - 1/(4\pi)$');
shading interp; view([-35 30]); grid on;

%% part c
k = [0:N/2-1, 0, -N/2+1:-1]; 
[KX, KY] = meshgrid(k, k);
K2 = KX.^2 + KY.^2;
 
F_hat = fft2(f);
 
Phi_hat = zeros(N, N);
mask = K2 > 0;
Phi_hat(mask) = F_hat(mask) ./ (-K2(mask));
Phi_hat(1,1) = 0;
 
phi = real(ifft2(Phi_hat));

%% Plot
figure('Name','(c) Solution phi');

surf(X, Y, phi, 'EdgeColor','none'); colorbar;
xlabel('x'); ylabel('y'); zlabel('$\phi$');
shading interp; view([-35 30]); grid on;
 