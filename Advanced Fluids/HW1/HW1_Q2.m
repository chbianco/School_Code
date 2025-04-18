clear all
close all

ar = linspace(0.6392,2,1000);
cr1 = sqrt( (1-1./(2.*ar)).^2 - (1./(4.*ar.^2)).*exp(-4.*ar));
cr2 = -sqrt((1-1./(2.*ar)).^2 - (1./(4.*ar.^2)).*exp(-4.*ar));

xr = linspace(0, 0.6392, 1000);
zero = linspace(0,0,1000);

ai = linspace(.0001, 0.6392,1000);
ci1 = sqrt((1-1./(2.*ai)).^2 - (1./(4.*ai.^2)).*exp(-4.*ai));
ci2 = -sqrt((1-1./(2.*ai)).^2 - (1./(4.*ai.^2)).*exp(-4.*ai));

xi = linspace(.6392, 2, 1000);

figure(1)
plot(ar, cr1, 'LineWidth', 2, 'Color', 'b')
hold on
plot(ar, cr2, 'LineWidth', 2, 'Color', 'b')
hold on
plot(xr, zero,'LineWidth', 2, 'Color', 'b')
ylim([-1,1])
grid on;
xlabel('$\alpha$','Interpreter','Latex','FontSize',12);
ylabel('$c_r$','Interpreter','latex','FontSize',12);
%title('Lift Comparison, Bands with Wind','Interpreter','latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',16)
set(gcf,'color','w')

figure(2)
plot(ai, -1i.*ci1, 'LineWidth', 2, 'Color', 'b')
hold on
plot(ai, -1i.*ci2, 'LineWidth', 2, 'Color', 'b')
hold on
plot(xi, zero,'LineWidth', 2, 'Color', 'b')
ylim([-1,1])
grid on;
xlabel('$\alpha$','Interpreter','Latex','FontSize',12);
ylabel('$c_i$','Interpreter','latex','FontSize',12);
%title('Lift Comparison, Bands with Wind','Interpreter','latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',16)
set(gcf,'color','w')


theta = linspace(0, pi, 100); % Angles from 0 to pi (half-circle)
x = cos(theta); % X-coordinates
y = sin(theta); % Y-coordinates

figure(3)
plot(cr1, -1i.*ci1,'LineWidth', 2, 'Color', 'b')
hold on
plot(x, y, 'b', 'LineWidth', 2, 'Color', 'r');
ylim([0,1.5])
xlim([-1, 1.5])
grid on;
xlabel('$c_r$','Interpreter','Latex','FontSize',12);
ylabel('$c_i$','Interpreter','latex','FontSize',12);
%title('Lift Comparison, Bands with Wind','Interpreter','latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',16)
set(gcf,'color','w')
legend('Eigenvalues', 'Howards Semicircle')