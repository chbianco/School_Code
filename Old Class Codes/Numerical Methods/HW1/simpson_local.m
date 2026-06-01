function S = simpson_local(fa,fm,fb,a,b)

h = (b-a)/2;
S = h/3*(fa + 4*fm + fb);

end