function yp = rayleigh(y,f)

global alpha c U Upp

U = tanh(y);
Upp = -2*tanh(y)*sech(y)^2;

yp = [f(2); f(1)*(Upp/(U-c) + alpha^2) ];