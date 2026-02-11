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

%% Solve part a
f = @(y) 1-y.^8;
N = 32;
a=0.98;
j = 0:1:N;

z = -1 + 2.*j/N;
x = (1/a).*tanh(z.*atanh(a));
h = abs(z(1) - z(2));
gp = @(y) a./(atanh(a).*(1-(a.^2).*y.^2));

exact = -8.*x.^7;
cd = zeros(1, N+1);

for k = 2:N
    cd(k)=gp(x(k))*(f(x(k+1))-f(x(k-1)))/(2*h);
end
%Boundaries 
cd(1) = gp(x(1)) * (-3*f(x(1)) + 4*f(x(2)) - f(x(3))) / (2*h);
cd(N+1) = gp(x(N+1)) * (3*f(x(N+1)) - 4*f(x(N)) + f(x(N-1))) / (2*h);

%% Solve part b

zb = pi.*j/N;
xb = cos(zb);
hb = abs(zb(1) - zb(2));
gp_b = @(y) -1./(sqrt(1-y.^2));
cd_b = zeros(1, N+1);

for k = 2:N
    cd_b(k)=gp_b(xb(k))*(f(xb(k+1))-f(xb(k-1)))/(2*hb);
end
%Boundaries 
cd_b(1) = gp_b(xb(1)) * (-3*f(xb(1)) + 4*f(xb(2)) - f(xb(3))) / (2*hb);
cd_b(N+1) = gp_b(xb(N+1)) * (3*f(xb(N+1)) - 4*f(xb(N)) + f(xb(N-1))) / (2*hb);



%% Plot
figure(1); clf; hold on; grid on
xlabel('$x$')
ylabel('$f^\prime$')
title('Part a scheme')
plot(x, cd)
plot(x, exact)
legend({'Central Difference', 'Exact'}, 'Location','best')

figure(2); clf; hold on; grid on
xlabel('$x$')
ylabel('$f^\prime$')
title('Part b scheme')
plot(xb, cd_b)
plot(x, exact)
legend({'Central Difference', 'Exact'}, 'Location','best')