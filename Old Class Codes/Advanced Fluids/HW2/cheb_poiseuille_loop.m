clear all
close all

% CHANGE THIS
n = 150;
Re = 6000;
beta = 0;

alph_vec = linspace(0.6,1.4, 50);

ci_vec = zeros(length(alph_vec),1);
cr_vec = zeros(length(alph_vec),1);

for i = 1:length(alph_vec)
alp = alph_vec(i);

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

ci_vec(i) = max(imag(eigvals));
cr_vec(i) = max(real(eigvals));

end

scatter(alph_vec, ci_vec)
xlabel('$\alpha$','Interpreter','Latex','FontSize',12);
ylabel('$c_i$','Interpreter','latex','FontSize',12);


