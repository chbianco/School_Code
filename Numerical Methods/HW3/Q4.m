%% Preamble
close all; clear variables; clc;
set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);

%% Grid parameters
Nx = 11; Ny = 11;
hx = 2 / (Nx - 1); 
hy = 1 / (Ny - 1);
x  = linspace(0, 2, Nx);
y  = linspace(0, 1, Ny);

%% Interior indices
xi = 2:Nx-1;
yi = 2:Ny-1;
Nxi = length(xi);
Nyi = length(yi);
N = Nxi * Nyi;

idx = @(p,q) (q-1)*Nxi + p;

%% Assemble sparse system Au = b
A = sparse(N, N);
b = zeros(N, 1);

rx = 1/hx^2;  % coefficient for x-terms
ry = 1/hy^2;  % coefficient for y-terms

for q = 1:Nyi % local y-index (interior)
    j = yi(q); % global j
    for p = 1:Nxi % local x-index (interior)
        i = xi(p); % global i
        k = idx(p,q);

        A(k,k) = A(k,k) - 2*rx;

        if p > 1
            A(k, idx(p-1,q)) = A(k, idx(p-1,q)) + rx;
        else
            b(k) = b(k) - rx * 0; % u_left = 0
        end

        if p < Nxi
            A(k, idx(p+1,q)) = A(k, idx(p+1,q)) + rx;
        else
            b(k) = b(k) - rx * y(j);
        end

        if q > 1 && q < Nyi
            A(k,k) = A(k,k) - 2*ry;
            A(k, idx(p,q-1)) = A(k, idx(p,q-1)) + ry;
            A(k, idx(p,q+1)) = A(k, idx(p,q+1)) + ry;
        elseif q == 1
            A(k,k) = A(k,k) - 2*ry;
            A(k, idx(p,q+1)) = A(k, idx(p,q+1)) + 2*ry;
        elseif q == Nyi
            A(k,k) = A(k,k) - 2*ry;
            A(k, idx(p,q-1)) = A(k, idx(p,q-1)) + 2*ry;
        end
    end
end

%% Solve
u_vec = A \ b;

%% Reshape into 2D grid (interior only), then fill full grid
U = zeros(Ny, Nx);

% Fill interior
for q = 1:Nyi
    for p = 1:Nxi
        U(yi(q), xi(p)) = u_vec(idx(p,q));
    end
end

% Apply Dirichlet boundary conditions
U(:, 1)  = 0;           % left:  u = 0
U(:, Nx) = y(:);        % right: u = y

%% Exact solution
u_exact = @(xv,yv) xv/4 - 4*sum( ...
    (1./((2*(1:50)-1)*pi).^2 ./ sinh(2*(2*(1:50)-1)*pi)) .* ...
    sinh((2*(1:50)-1)*pi*xv) .* cos((2*(1:50)-1)*pi*yv), 2);

[XX, YY] = meshgrid(x, y);
U_exact  = arrayfun(@(xv,yv) u_exact(xv,yv), XX, YY);

%% Plots
figure(1);
subplot(1,2,1)
surf(XX, YY, U, 'EdgeColor','none');
xlabel('x'); ylabel('y'); zlabel('u');
title('FD Numerical Solution'); colorbar;

subplot(1,2,2)
surf(XX, YY, U_exact, 'EdgeColor','none');
xlabel('x'); ylabel('y'); zlabel('u');
title('Exact Solution'); colorbar;

figure(2);
contourf(XX, YY, abs(U - U_exact), 20);
xlabel('x'); ylabel('y');
title('Pointwise Error |u_{FD} - u_{exact}|'); colorbar;

fprintf('Max error: %.6f\n', max(abs(U(:) - U_exact(:))));