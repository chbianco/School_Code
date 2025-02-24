function [y, f] = cheb_example(eps, n)

% Example of the use of Chebyshev polynomials to solve a
% differential equation.
%
% eps is the small parameter in the equation 
% n   is the number of Chebyshev modes
%
% Updated, Kenny Breuer Feb 2025

[y, t, t1, t2, t3, t4] = setup_cheb(n);  

% Initialize Matrices
B = zeros(n+1,n+1);

% Fill in with equations, in this case:
%
% eps * y'' + (1 + eps) y' + y = 0 
%
%
% We are setting up the equation B a = c, where a are the coefficients
B = eps*t2 + (1 + eps)*t1 + t;

% Boundary conditions
% f = 1 @ y =  1
% f = 0 @ y = -1
B(1,  :)   =  t(1,  :);
B(n+1,:)   =  t(n+1,:);

c      = zeros(n+1,1);
c(1)   = 1;
c(n+1) = 0;

% Solve the linear system
a = B\c;

% Evaluate the soution on the collocation points
f = t*a;

plot(y, f, '-o');
set(gca, 'FontSize', 16)
xlabel('x')
ylabel('y')
grid on
