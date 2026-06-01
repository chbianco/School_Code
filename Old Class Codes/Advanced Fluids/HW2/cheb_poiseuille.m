clear all
close all
%% Doing the math
% CHANGE THIS
n = 150;
Re = 6000;
alp = 1; %Least stable alpha for Re = 2000 is 1.5102
beta = 0;

zi=sqrt(-1);
% mean velocity
ak2=alp^2+beta^2;
Nos=n+1;
Nsq=n+1 ;
vec=(0:1:n)'; 
u = (ones(length(vec),1)-cos(pi*vec/n).^2);
du=-2*cos(pi*vec/n); 

% Set up Chebyshev collocation points and differentiation matrices
[y, D0, D1, D2, D3, D4] = setup_cheb(n);

% set up Orr-Sommerfeld matrix
B11=D2-ak2*D0;
A11=-(D4-2*ak2*D2+(ak2^2)*D0)/(zi*Re);
A11=A11+alp*(u*ones(1,length(u))).*B11+alp*2*D0;
er=-200*zi;
A11=[er*D0(1,:); er*D1(1,:); A11(3:Nos-2,:); er*D1(Nos,:); er*D0(Nos,:)];
B11=[D0(1,:); D1(1,:); B11(3:Nos-2,:);D1(Nos,:); D0(Nos,:)]; 

% set up Squire matrix and cross-term matrix
A21=beta*(du*ones(1,length(u))).*D0(1:Nos,:);
A22=alp*(u*ones(1,length(u))).*D0-(D2-ak2*D0)/(zi*Re);
B22=D0;
A22=[er*D0(1,:); A22(2:Nsq-1,:); er*D0(Nsq,:)];
A21=[zeros(1,Nos); A21(2:Nsq-1,:); zeros(1,Nos)];

% combine all the blocks
A=[A11 zeros(Nos,Nsq); A21 A22];
B=[B11 zeros(Nos,Nsq); zeros(Nsq,Nos) B22];

% Solve the eigenvalue problem
d=inv(B)*A;
[eigvecs, eigvals_matrix] = eig(d);
eigvals = diag(eigvals_matrix);

%Sorting our OS vs Squire
% ci_comp = sort(imag(eigvals), 'desc');
% cr_comp = sort(real(eigvals), 'desc');
% 
% ci_os = ci_comp;
% ci_os(1:2:end) = [];

%% Calculate Velocities
modes = 5; %number of modes you want
[~, idx] = maxk(imag(eigvals), modes);

for j = 1:modes

    tk = eigvecs(:, idx(j));
    tk = tk(tk ~= 0);
    vhat = D0*tk;
    vprime = D1*tk;
    uhat = vprime.*(1i*alp);

    figure(2)
    plot(abs(uhat), y, 'LineWidth', 2, 'DisplayName',['Mode ', num2str(j)])
    hold on
    
    figure(3)
    plot(abs(vhat), y, 'LineWidth', 2, 'DisplayName',['Mode ', num2str(j)])
    hold on
   

end

figure(2)
xlabel('$|\hat{u}|$','Interpreter','Latex','FontSize',12);
ylabel('$y$','Interpreter','latex','FontSize',12);
% xlim([0,2])
% ylim([-1, 0])
title('x-velocity','Interpreter','latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',16)
legend('show')

figure(3)
xlabel('$|\hat{v}|$','Interpreter','Latex','FontSize',12);
ylabel('$y$','Interpreter','latex','FontSize',12);
% xlim([0,2])
% ylim([-1, 0])
title('y-velocity','Interpreter','latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',16)
legend('show')

%% Plotting
figure(1)
scatter(real(eigvals), imag(eigvals))
xlabel('$c_r$','Interpreter','Latex','FontSize',12);
ylabel('$c_i$','Interpreter','latex','FontSize',12);
xlim([0,1])
ylim([-1, 0])
title(['Eigenvalue Spectrum for $\alpha$ = ',num2str(alp), ', Re = ',num2str(Re)],'Interpreter','latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',16)
% %set(gcf,'color','w', 'Position', [10 10 900 600])