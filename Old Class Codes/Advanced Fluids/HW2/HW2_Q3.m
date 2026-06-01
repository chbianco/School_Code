close all
clear all
t = linspace(0, 2000, 100000);

Re = [100, 500, 1000];



for i = 1:3
v = exp(-t./Re(i));
n = Re(i).*exp(-t./Re(i)) - Re(i).*exp(-2.*t./Re(i));

figure(1)
plot(t, v, 'LineWidth', 2, 'DisplayName',['Re = ', num2str(Re(i))])
hold on

figure(2)
plot(t, n, 'LineWidth', 2,  'DisplayName',['Re = ', num2str(Re(i))])
hold on
end

figure(1)
xlabel('$t$','Interpreter','Latex','FontSize',12);
ylabel('$v$','Interpreter','latex','FontSize',12);
title('Evolution of $v(t)$','Interpreter','latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',16)
set(gcf,'color','w', 'Position', [10 10 900 600])
legend('show')

figure(2)
xlabel('$t$','Interpreter','Latex','FontSize',12);
ylabel('$\eta$','Interpreter','latex','FontSize',12);
title('Evolution of $\eta(t)$','Interpreter','latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',16)
set(gcf,'color','w', 'Position', [10 10 900 600])
legend('show')