function [S, xeval] = recurse(f,a,b,fa,fm,fb,S,tol,xeval)

m = (a+b)/2;
lm = (a+m)/2;
rm = (m+b)/2;

flm = f(lm);
frm = f(rm);

xeval = [xeval lm rm];

Sleft  = simpson_local(fa,flm,fm,a,m);
Sright = simpson_local(fm,frm,fb,m,b);

if abs(Sleft + Sright - S) < 15*tol
    S = Sleft + Sright + (Sleft + Sright - S)/15;
else
    [Sleft,  xeval] = recurse(f,a,m,fa,flm,fm,Sleft,tol/2,xeval);
    [Sright, xeval] = recurse(f,m,b,fm,frm,fb,Sright,tol/2,xeval);
    S = Sleft + Sright;
end

end