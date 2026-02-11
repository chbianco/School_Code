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
exact = -25*sin(5*1.5);
f = @(x) sin(5*x);

h = logspace(-4, 0, 1000);

fpp_pop = zeros(1, length(h));
fpp_pop_e = zeros(1, length(h));

fpp_new = zeros(1, length(h));
fpp_new_e = zeros(1, length(h));

for j = 1:length(h)

fpp_pop(j) = (f(1.5 + h(j)) - 2*f(1.5) + f(1.5 - h(j)))/h(j)^2;
fpp_pop_e(j) = abs(fpp_pop(j) - exact);

fpp_new(j) = (f(1.5 + 2*h(j)) - 2*f(1.5) + f(1.5 - 2*h(j)))/(4*h(j)^2);
fpp_new_e(j) = abs(fpp_new(j) - exact);

end

%% Plot
figure(1); clf; hold on; grid on
set(gca,'XScale','log','YScale','log') 
xlabel('h')
ylabel('error')
loglog(h, fpp_pop_e)
loglog(h, fpp_new_e)

h_ref = h(end);
C2 = 10 / h_ref^2;  % 2nd order
loglog(h, C2*h.^2, '--k')


legend({'Popular (given) formula', 'Derived formula', 'Second order power law'}, 'Location','best')