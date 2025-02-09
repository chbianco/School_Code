function yp = rayleigh(y,f)

global alpha c U Upp

%Solve for tanh profile
% U = tanh(y);
% Upp = -2*tanh(y)*sech(y)^2;

%Solve for sech profile
U = sech(y);
Upp = sech(y)*tanh(y)^2 - sech(y)^3;

yp = [f(2); f(1)*(Upp/(U-c) + alpha^2) ];