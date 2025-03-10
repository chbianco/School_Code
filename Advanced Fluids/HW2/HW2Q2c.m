n = 150;

Re_vec = [500, 2000, 6000];
alph_vec = linspace(0.6, 1.4, 30);
beta = 0;

for j = 1:3
    Re = Re_vec(j);
    ci_vec = zeros(length(alph_vec),1);
    cr_vec = zeros(length(alph_vec),1);

    for m = 1:length(alph_vec)
        alp = alph_vec(m);

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
        [eigvecs, eigvals_matrix] = eig(A, B);
        eigvals = diag(eigvals_matrix);

        eigvals = eigvals(real(eigvals) >= 0 & real(eigvals) <= 1.5);


        ci_vec(m) = max(imag(eigvals));
        cr_vec(m) = max(real(eigvals));

    end

    figure(1)
    plot(alph_vec, ci_vec, 'LineWidth', 2, 'DisplayName',['Re = ', num2str(Re_vec(j))])
    hold on

    figure(2)
    plot(alph_vec, cr_vec, 'LineWidth', 2,  'DisplayName',['Re = ', num2str(Re_vec(j))])
    hold on
end

figure(1)
xlabel('$\alpha$','Interpreter','Latex','FontSize',12);
ylabel('$c_i$','Interpreter','latex','FontSize',12);
title('Imaginary Dispersion Relation','Interpreter','latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',16)
set(gcf,'color','w', 'Position', [10 10 900 600])
legend('show')

figure(2)
xlabel('$\alpha$','Interpreter','Latex','FontSize',12);
ylabel('$c_r$','Interpreter','latex','FontSize',12);
title('Real Dispersion Relation','Interpreter','latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',16)
set(gcf,'color','w', 'Position', [10 10 900 600])
legend('show')
