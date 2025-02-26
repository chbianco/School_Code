% CHANGE THIS
n = 100;
Re = 20000;

% Set up Chebyshev collocation points and differentiation matrices
[y, t, t1, t2, t3, t4] = setup_cheb(n);

% Define the base flow profile U(y) = 1 - y^2 and its derivatives
U = 1 - y.^2;
d2U = -2 * ones(size(y)); % Second derivative of U

% Initialize matrices
A = zeros(n+1, n+1);
B = zeros(n+1, n+1);

% Define the base flow profile U(y) and its derivatives (for interior points)
U_int = U(2:n);
d2U_int = d2U(2:n);

% Convert to diagonal matrices
U_diag = diag(U_int);
d2U_diag = diag(d2U_int);

% Ensure identity matrix has correct size
I = eye(n+1);

% Construct the Orr-Sommerfeld matrix (interior points)
A(2:n, 2:n) = (t4(2:n,2:n) - 2 * t2(2:n,2:n) + I(2:n,2:n)) / Re ...
               - U_diag * (t2(2:n,2:n) - I(2:n,2:n)) - d2U_diag;

B(2:n, 2:n) = t2(2:n, 2:n) - I(2:n, 2:n);

% Apply boundary conditions
A(1, :) = t(1, :);  % u = 0 at y = 1
A(n+1, :) = t(n+1, :); % u = 0 at y = -1
A(2, :) = t1(1, :);  % du/dy = 0 at y = 1
A(n, :) = t1(n+1, :);   % du/dy = 0 at y = -1

B(1, :) = 0;  % No contribution to the eigenvalue problem
B(2, :) = 0;
B(n, :) = 0;
B(n+1, :) = 0;

% Solve the eigenvalue problem
[eigvecs, eigvals_matrix] = eig(A, B);
eigvals = diag(eigvals_matrix);

figure(1)
scatter(real(eigvals), imag(eigvals))

xlabel('$Real$','Interpreter','Latex','FontSize',12);
ylabel('$Imaginary$','Interpreter','latex','FontSize',12);
title('Evolution of $v(t)$','Interpreter','latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',16)
set(gcf,'color','w', 'Position', [10 10 900 600])
