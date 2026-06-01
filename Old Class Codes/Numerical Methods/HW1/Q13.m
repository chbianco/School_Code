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

%% Solve
N_vec = [8,16,32,64,128,256,512];
f = @(y) 100./(sqrt(y+0.01)) +1./((y-0.3).^2 + .001) - pi;
exact = 272.4531323160361;
err_trap = zeros(size(N_vec));
err_sim = err_trap;
err_trapcorr = err_trap;

for i = 1: length(N_vec)
    N = N_vec(i);
    h = 1/N;
    x = linspace(0,1,N+1);
    fx = f(x);
    
    %Trapezoidal
    trap = h*(0.5*fx(1) + sum(fx(2:end-1)) + 0.5*fx(end));
    err_trap(i) = abs(((trap - exact)/exact))*100; 
    
    %Simpson
    sim = h/3*(fx(1)+fx(end)+4*sum(fx(2:2:end-1))+2*sum(fx(3:2:end-2)));
    err_sim(i) = abs((sim-exact)/exact)*100;

    %Trapezoidal w/ end correction
    fp = @(y) -50*(y+0.01).^(-3/2) - 2*(y-0.3)./(( (y-0.3).^2 + 0.001 ).^2);
    trapcorr = trap - h^2/12*(fp(1) - fp(0));
    err_trapcorr(i) = abs((trapcorr-exact)/exact)*100;

    
end

%% Adaptive Simpson
a = 0; 
b = 1;

tols = [1e-2 1e-4 1e-6 1e-8];

for k = 1:length(tols)

    tol = tols(k);
    xeval = []; 

    [I, xeval] = adaptsimp(f,a,b,tol,xeval);

    fprintf('tol = %e,  I = %.12f,  Npoints = %d\n', ...
        tol, I, length(unique(xeval)));

    %Plot 
    figure(k)
    xx = linspace(0,1,2000);
    plot(xx,f(xx),'b','LineWidth',1.5)
    hold on
    plot(unique(xeval), f(unique(xeval)), 'ro','MarkerSize',6)
    xlabel('y')
    ylabel('f(y)')
    title(['Adaptive Simpson, tol = ' num2str(tol)])
    grid on

end



%% Plot
figure(5); hold on;
grid on
loglog(N_vec, err_trap)
loglog(N_vec, err_sim)
loglog(N_vec, err_trapcorr)
xlabel('N')
ylabel('Percent error')
legend({'Trapezoidal', 'Simpson', 'Trapezoidal w/ correction'})