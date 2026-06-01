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

%% Solve
N = 24; 
x = linspace(0,1,N+1);
h = x(2)-x(1);
x = x + h/2;
x = x(1: end-1);

I = eye(24,24);
e = ones(N-1,1);
% Create A matrix
A = diag(10/12*ones(N,1)) + diag((1/12)*e,1) + diag((1/12)*e,-1);
A(1,1)=1;
A(1,2) = -11/23;
A(end, end) = 1;
A(end, end-1) = -11/23;
%Create B matrix
B = diag(-2*ones(N,1)) + diag(e,1) + diag(e,-1);
B(1,1)=-36/23;
B(1,2) = 48/23;
B(1,3) = -12/23;
B(end, end) = -36/23;
B(end, end-1) = 48/23;
B(end,end-2) = -12/23;
B = (1/h^2)*B;

y = (A\B + I)\(x.^3)';

%% Plotting

exact = -6.*x + x.^3 + 3*(-1 + 2*cos(1)).*cos(x)*csc(1) + 6.*sin(x);

figure(1); hold on
plot(x, exact, 'r')
plot(x,y, 'b--')
legend({'Numerical', 'Analytical'}, 'Location','best')
