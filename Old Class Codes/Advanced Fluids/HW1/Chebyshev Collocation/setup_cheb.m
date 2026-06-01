function [y, t, t1, t2, t3, t4, t5, t6] = setup_cheb(n)

% Setup Chebyshev collocation points, polynomials and derivatives

% Chebyshev polynomials, generated on collocation points
y  = cos((0:n)*pi/n)';

t  = zeros(n+1,n+1);
t1 = zeros(n+1,n+1);
t2 = zeros(n+1,n+1);
t3 = zeros(n+1,n+1);
t4 = zeros(n+1,n+1);
t5 = zeros(n+1,n+1);
t6 = zeros(n+1,n+1);

for k=1:n+1,
  for j=1:n+1,
    t(j,k) = cos((j-1)*(k-1)*pi/n);
  end
end

% First derivative:
for j = 1:n+1,
  t1(j,1) = 0;
  t1(j,2) =    t(j,1);
  t1(j,3) = 4.*t(j,2);
  for k = 4:n+1,
    t1(j,k)=2*(k-1)*t(j,k-1) + (k-1)*t1(j,k-2)/(k-3);
  end
end

% Second Derivative:
for j = 1:n+1,
  t2(j,1) = 0;
  t2(j,2) =    t1(j,1);
  t2(j,3) = 4.*t1(j,2);
  for k = 4:n+1,
    t2(j,k)=2*(k-1)*t1(j,k-1) +(k-1) *t2(j,k-2)/(k-3);
  end
end

% Third Derivative:
for j = 1:n+1,
  t3(j,1) = 0;
  t3(j,2) =    t2(j,1);
  t3(j,3) = 4.*t2(j,2);
  for k = 4:n+1,
    t3(j,k)=2*(k-1)*t2(j,k-1)+(k-1)*t3(j,k-2)/(k-3);
  end
end

% Fourth Derivative:
for j = 1:n+1,
  t4(j,1) = 0;
  t4(j,2) =    t3(j,1);
  t4(j,3) = 4.*t3(j,2);
  for k = 4:n+1,
    t4(j,k)=2*(k-1)*t3(j,k-1)+(k-1)*t4(j,k-2)/(k-3);
  end
end

% Fifth Derivative:
for j = 1:n+1,
  t5(j,1) = 0;
  t5(j,2) =    t4(j,1);
  t5(j,3) = 4.*t4(j,2);
  for k = 4:n+1,
    t5(j,k)=2*(k-1)*t4(j,k-1)+(k-1)*t5(j,k-2)/(k-3);
  end
end

% Sixth Derivative:
for j = 1:n+1,
  t6(j,1) = 0;
  t6(j,2) =    t5(j,1);
  t6(j,3) = 4.*t5(j,2);
  for k = 4:n+1,
    t6(j,k)=2*(k-1)*t5(j,k-1)+(k-1)*t6(j,k-2)/(k-3);
  end
end

