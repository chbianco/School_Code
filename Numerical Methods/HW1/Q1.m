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
f = @(x) sin(2.*x)*cos(20.*x)+exp(sin(2.*x));
fp = @(x) -20*sin(2.*x).*sin(20.*x)+2.*cos(2.*x).*cos(20.*x)+2.*exp(sin(2.*x).*cos(2.*x));

fo_fd = cell(4,1);
so_fd = cell(4,1);
fo_cd = cell(4,1);
fo_pa = cell(4,1);

e_fo_fd = [0,0,0,0];
e_so_fd = [0,0,0,0];
e_fo_cd = [0,0,0,0];
e_fo_pa = [0,0,0,0];

xis = cell(4,1);
hs = [0,0,0,0];

exact = cell(4,1);

for n = 1:4
    N = Ns(n);
    h = 2*pi/N;
    hs(n) = h;
    xis{n} = 0: h: 2*pi;
    xi = xis{n};

    exact = fp(xi);

    
    fo_fd{n} = 0.*xi;
    so_fd{n} = 0.*xi;
    fo_cd{n} = 0.*xi;
    b = 0.*xi;

    for j = 1:length(xi)-1 %First order finite difference
        fo_fd{n}(j) = (f(xi(j+1)) - f(xi(j)))/h; 
    end    
    %Impose periodicity
    fo_fd{n}(end) = (f(xi(2)) - f(xi(end)))/h;
    %error
    e_fo_fd(n) = max(abs(fo_fd{n} - exact));

    for j = 2:length(xi)-1 %Second order finite difference
        so_fd{n}(j) = (f(xi(j+1)) - f(xi(j-1)))/(2*h); 
    end
    %Impose periodicity
    so_fd{n}(1) = (f(xi(2)) - f(xi(end-1)))/(2*h); 
    so_fd{n}(end) = so_fd{n}(1);
    %error
    e_so_fd(n) = max(abs(so_fd{n} - exact));

    for j = 3:length(xi)-2 %Fourth order central difference
        fo_cd{n}(j) = (f(xi(j-2)) - 8*f(xi(j-1)) + 8*f(xi(j+1))- f(xi(j+2)))/(12*h); 
    end
    %Impose periodicity
    fo_cd{n}(1) = (f(xi(end-1)) - 8*f(xi(end)) + 8*f(xi(2))- f(xi(3)))/(12*h); 
    fo_cd{n}(2) = (f(xi(end-3)) - 8*f(xi(1)) + 8*f(xi(3))- f(xi(4)))/(12*h);
    fo_cd{n}(end) = (f(xi(end-2)) - 8*f(xi(end-1)) + 8*f(xi(1))- f(xi(2)))/(12*h);
    fo_cd{n}(end-1) = (f(xi(end-3)) - 8*f(xi(end-2)) + 8*f(xi(end))- f(xi(1)))/(12*h);
    %error
    e_fo_cd(n) = max(abs(fo_cd{n} - exact));

    for j = 2:length(xi)-1 %Fourth order pade
        b(j) = (3/h).*(f(xi(j+1))-f(xi(j-1)));
    end
    b(1) = (3/h).*(f(xi(2)) - f(xi(end)));
    b(end) = (3/h).*(f(xi(1)) - f(xi(end-1)));

    e = ones(N+1,1);
    A = diag(4*e) + diag(e(1:end-1),1) + diag(e(1:end-1),-1);
    fo_pa{n} = A \ b';
    %error
    e_fo_pa(n) = max(abs(fo_pa{n}' - exact));

end


%% Plots
figure(1); hold on;
grid on
xlabel('x')
ylabel('$f^\prime(x)$')
plot(xis{1}, fo_fd{1}, '-b')
plot(xis{2}, fo_fd{2}, '--b')
plot(xis{3}, fo_fd{3}, ':b')
plot(xis{4}, fo_fd{4}, '-.b')

plot(xis{1}, so_fd{1}, '-r')
plot(xis{2}, so_fd{2}, '--r')
plot(xis{3}, so_fd{3}, ':r')
plot(xis{4}, so_fd{4}, '-.r')

plot(xis{1}, so_fd{1}, '-c')
plot(xis{2}, so_fd{2}, '--c')
plot(xis{3}, so_fd{3}, ':c')
plot(xis{4}, so_fd{4}, '-.c')

plot(xis{1}, fo_pa{1}, '-k')
plot(xis{2}, fo_pa{2}, '--k')
plot(xis{3}, fo_pa{3}, ':k')
plot(xis{4}, fo_pa{4}, '-.k')

plot(xis{4}, exact, 'y')

figure(2); clf; hold on; grid on
set(gca,'XScale','log','YScale','log') 
loglog(hs, e_fo_fd)
loglog(hs, e_fo_cd)
loglog(hs, e_so_fd)
loglog(hs, e_fo_pa)
legend({'FD1', 'CD4', 'FD2', 'PA4'})

