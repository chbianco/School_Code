%% PREAMBLE
close all;
clear variables;
clc;

set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);

%% Solving 

Ns = [256, 512, 1024, 2048];
f = @(x) sin(2.*x).*cos(20.*x)+exp(sin(2.*x));
fp = @(x) -20*sin(2.*x).*sin(20.*x)+2.*cos(2.*x).*cos(20.*x)+2*cos(2.*x).*exp(sin(2.*x));

fo_fd = cell(4,1);
so_fd = cell(4,1);
fo_cd = cell(4,1);
fo_pa_to = cell(4,1);
fo_pa_pe = cell(4,1);

e_fo_fd = [0,0,0,0];
e_so_fd = [0,0,0,0];
e_fo_cd = [0,0,0,0];
e_fo_pa_to = [0,0,0,0];
e_fo_pa_pe = [0,0,0,0];


xis = cell(4,1);
hs = [0,0,0,0];

exact = cell(4,1);

for n = 1:4
    N = Ns(n);
    h = 2*pi/N;
    hs(n) = h;
    xis{n} = (0:N-1)*h;
    xi = xis{n};

    exact{n} = fp(xi);

    
    fo_fd{n} = 0.*xi;
    so_fd{n} = 0.*xi;
    fo_cd{n} = 0.*xi;
    b_to = 0.*xi;
    b_pe = 0.*xi;

    for j = 1:length(xi)-1 %First order finite difference
        fo_fd{n}(j) = (f(xi(j+1)) - f(xi(j)))/h; 
    end    
    %Impose periodicity
    fo_fd{n}(end) = (f(xi(1)) - f(xi(end)))/h;
    %error
    e_fo_fd(n) = max(abs(fo_fd{n} - exact{n}));

    for j = 2:length(xi)-1 %Second order finite difference
        so_fd{n}(j) = (f(xi(j+1)) - f(xi(j-1)))/(2*h); 
    end
    %Impose periodicity
    so_fd{n}(1) = (f(xi(2)) - f(xi(end)))/(2*h); 
    so_fd{n}(end) = (f(xi(1)) - f(xi(end-1)))/(2*h); 
    %error
    e_so_fd(n) = max(abs(so_fd{n} - exact{n}));

    for j = 3:length(xi)-2 %Fourth order central difference
        fo_cd{n}(j) = (f(xi(j-2)) - 8*f(xi(j-1)) + 8*f(xi(j+1))- f(xi(j+2)))/(12*h); 
    end
    %Impose periodicity
    fo_cd{n}(1) = (f(xi(end-1)) - 8*f(xi(end)) + 8*f(xi(2))- f(xi(3)))/(12*h); 
    fo_cd{n}(2) = (f(xi(end)) - 8*f(xi(1)) + 8*f(xi(3))- f(xi(4)))/(12*h);
    fo_cd{n}(end) = (f(xi(end-2)) - 8*f(xi(end-1)) + 8*f(xi(1))- f(xi(2)))/(12*h);
    fo_cd{n}(end-1) = (f(xi(end-3)) - 8*f(xi(end-2)) + 8*f(xi(end))- f(xi(1)))/(12*h);
    %error
    e_fo_cd(n) = max(abs(fo_cd{n} - exact{n}));

    for j = 2:length(xi)-1 %Fourth order pade
        b_to(j) = (3/h).*(f(xi(j+1))-f(xi(j-1)));
        b_pe(j) = (3/h).*(f(xi(j+1))-f(xi(j-1)));
    end
    %Periodicity BC
    b_pe(1) = (3/h).*(f(xi(2)) - f(xi(end)));
    b_pe(end) = (3/h).*(f(xi(1)) - f(xi(end-1)));
    %Third order BC
    b_to(1) = (1/h)*(-(5/2).*f(xi(1)) + 2.*f(xi(2)) + 0.5.*f(xi(3)));
    b_to(end) = (1/h)*((5/2).*f(xi(end)) - 2.*f(xi(end-1)) - 0.5.*f(xi(end-2)));

    e = ones(N,1);
    %Third order A
    A_to = diag(4*e) + diag(e(1:end-1),1) + diag(e(1:end-1),-1);
    % --- LEFT boundary row (j = 1)
    A_to(1,:) = 0;
    A_to(1,1) = 2;
    A_to(1,2) = 1;
    % --- RIGHT boundary row (j = N)
    A_to(end,:) = 0;
    A_to(end,end)   = 2;
    A_to(end,end-1) = 1;
    %Periodic A
    A_pe = diag(4*e) + diag(e(1:end-1),1) + diag(e(1:end-1),-1);
    A_pe(1,end) = 1;
    A_pe(end,1) = 1;
    fo_pa_to{n} = A_to \ b_to';
    fo_pa_pe{n} = A_pe \ b_pe';
    %error
    e_fo_pa_to(n) = max(abs(fo_pa_to{n}' - exact{n}));
    e_fo_pa_pe(n) = max(abs(fo_pa_pe{n}' - exact{n}));


end


%% Plots
figure(1); clf; hold on; grid on
set(gca,'XScale','log','YScale','log') 
xlabel('h')
ylabel('$L_\infty$ error')

loglog(hs, e_fo_fd, 'ob', 'MarkerSize', 20)
loglog(hs, e_fo_cd, '*c', 'MarkerSize', 20)
loglog(hs, e_so_fd, 'xr', 'MarkerSize', 20)
loglog(hs, e_fo_pa_to, '+k', 'MarkerSize', 20)
loglog(hs, e_fo_pa_pe, '^k', 'MarkerSize', 20)

h_ref = hs(end);

% Reference errors
C1 = e_fo_fd(end)    / h_ref^1;  % 1st order
C2 = e_so_fd(end)    / h_ref^2;  % 2nd order
C3 = e_fo_pa_to(end) / h_ref^1.75;  % 3rd order, pade third order
C4 = e_fo_cd(end)    / h_ref^4;  % 4th order, central difference
C5 = e_fo_pa_pe(end)    / h_ref^4;  % 4th order, pade periodic 

% Asymptote lines
loglog(hs, C1*hs.^1, '--b')
loglog(hs, C2*hs.^2, '--r')
loglog(hs, C3*hs.^1.75, '--k')
loglog(hs, C4*hs.^4, '--c')
loglog(hs, C5*hs.^4, '--k' )

legend({'FD1', 'CD4', 'FD2', 'PA4, third order BC', 'PA4, periodic BC'}, 'Location','southeast')

figure(2); clf
% ---- line styles 
ls_exact = 'k-';   lw_exact = 2;
ls_fd1   = 'b--';
ls_fd2   = 'r-.';
ls_cd4   = 'c-';
ls_pa3   = 'm:';
ls_pa4   = 'g-';

for n = 1:4
    subplot(2,2,n); hold on; grid on

    xi = xis{n};

    % Exact solution
    plot(xi, exact{n}, ls_exact, 'LineWidth', lw_exact)

    % Numerical schemes
    plot(xi, fo_fd{n},    ls_fd1)
    plot(xi, so_fd{n},    ls_fd2)
    plot(xi, fo_cd{n},    ls_cd4)
    plot(xi, fo_pa_to{n}, ls_pa3)
    plot(xi, fo_pa_pe{n}, ls_pa4)

    title(sprintf('N = %d,  h = %.2e', Ns(n), hs(n)))
    xlabel('x')
    ylabel('df/dx')

    if n == 1
        legend({ ...
            'Exact', ...
            'FD1', ...
            'FD2', ...
            'CD4', ...
            'Pade4 (3rd-order BC)', ...
            'Pade4 (periodic BC)' ...
        }, 'Location','best')
    end
end

