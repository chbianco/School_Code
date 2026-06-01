function [S, xeval] = adaptsimp(f,a,b,tol,xeval)

m = (a+b)/2;

fa = f(a);
fm = f(m);
fb = f(b);

xeval = [xeval a m b];

S = simpson_local(fa,fm,fb,a,b);

[S, xeval] = recurse(f,a,b,fa,fm,fb,S,tol,xeval);

end