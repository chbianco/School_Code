A  = [4 1 2 0; 1 3 0 1; 2 0 5 1; 0 1 1 2];
x = [1 1 1 1]';

for i = 1:1000
    x_new = A*x;
    x = x_new/norm(x_new);
end 

eig_max = norm(A*x)

eig_val = max(abs(eig(A)))

