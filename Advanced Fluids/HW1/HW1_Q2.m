ar = linspace(0.6392,2,1000);
cr1 = sqrt((1-1./(2.*a)).^2 - (1./(16.*a.^2)).*exp(-4.*a));
cr2 = -sqrt((1-1./(2.*a)).^2 - (1./(16.*a.^2)).*exp(-4.*a));

ai = linspace(0, 0.6392,1000);
ci1 = sqrt((1-1./(2.*a)).^2 - (1./(16.*a.^2)).*exp(-4.*a));
ci2 = -sqrt((1-1./(2.*a)).^2 - (1./(16.*a.^2)).*exp(-4.*a));

figure(1)
plot(ar, cr1, 'LineWidth', 2, 'Color', 'b')
hold on
plot(ar, cr2, 'LineWidth', 2, 'Color', 'b')
ylim([-1,1])
grid on;
xlabel('$\alpha$','Interpreter','Latex','FontSize',12);
ylabel('$c_r$','Interpreter','latex','FontSize',12);
%title('Lift Comparison, Bands with Wind','Interpreter','latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',16)
set(gcf,'color','w')

figure(2)
plot(ai, ci1, 'LineWidth', 2, 'Color', 'b')
hold on
plot(ai, ci2, 'LineWidth', 2, 'Color', 'b')
ylim([-1,1])
grid on;
xlabel('$\alpha$','Interpreter','Latex','FontSize',12);
ylabel('$c_i$','Interpreter','latex','FontSize',12);
%title('Lift Comparison, Bands with Wind','Interpreter','latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',16)
set(gcf,'color','w')