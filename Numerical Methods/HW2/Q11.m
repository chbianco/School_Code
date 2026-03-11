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

%% Setup

beta_I = 0.160;       
beta_H = 0.062;       
beta_F = 0.489;       
alpha = 1/12;       
gamma_H = 1/3.24;     
gamma_DH = 1/10.07;
gamma_F = 1/2.01;  
gamma_I = 1/15.00;  
gamma_D = 1/13.31; 
gamma_IH = 1/15.88;  
theta1 = 0.197;  
delta1 = 0.500; 
delta2 = 0.500;  
N = 4.3e6;  

E0 = 2;
S0 = N - E0;
I0 = 0; H0 = 0; F0 = 0; R0 = 0;
u0 = [S0; E0; I0; H0; F0; R0];

%% Jacobian
J = zeros(6,6);
J(1,3) = -beta_I*S0/N;
J(1,4) = -beta_H*S0/N;
J(1,5) = -beta_F*S0/N;
J(2,2) = -alpha;
J(2,3) =  beta_I*S0/N;
J(2,4) =  beta_H*S0/N;
J(2,5) =  beta_F*S0/N;
J(3,2) =  alpha;
J(3,3) = -gamma_H*theta1 + gamma_I*(1-theta1)*(1-delta1) ...
          + gamma_D*(1-theta1)*delta1;
J(4,3) =  gamma_H*theta1;
J(4,4) = -gamma_DH*delta2 + gamma_IH*(1-delta2);
J(5,3) =  gamma_D*(1-theta1)*delta1;
J(5,4) =  gamma_DH*delta2;
J(5,5) = -gamma_F;
J(6,3) =  gamma_I*(1-theta1)*(1-delta1);
J(6,4) =  gamma_IH*(1-delta2);
J(6,5) =  gamma_F;

eigs_J = eig(J);
min_eig = min(eigs_J)

%% Solving
%RHS
function dudt = ebola_rhs(u, params)
    S = u(1); E = u(2); I = u(3);
    H = u(4); F = u(5); R = u(6);

    p = params;
    lambda_force = (p.beta_I*I + p.beta_H*H + p.beta_F*F) / p.N;

    dS = -lambda_force * S;
    dE =  lambda_force * S - p.alpha * E;
    dI =  p.alpha * E ...
         - (p.gamma_H*p.theta1 ...
         +  p.gamma_I*(1-p.theta1)*(1-p.delta1) ...
         +  p.gamma_D*(1-p.theta1)*p.delta1) * I;
    dH =  p.gamma_H*p.theta1*I ...
         - (p.gamma_DH*p.delta2 + p.gamma_IH*(1-p.delta2)) * H;
    dF =  p.gamma_D*(1-p.theta1)*p.delta1*I ...
         + p.gamma_DH*p.delta2*H - p.gamma_F*F;
    dR =  p.gamma_I*(1-p.theta1)*(1-p.delta1)*I ...
         + p.gamma_IH*(1-p.delta2)*H + p.gamma_F*F;

    dudt = [dS; dE; dI; dH; dF; dR];
end

%Consolidate parameters
params.beta_I = beta_I;  params.beta_H = beta_H;
params.beta_F = beta_F;  params.alpha = alpha;
params.gamma_H = gamma_H; params.gamma_DH = gamma_DH;
params.gamma_F = gamma_F; params.gamma_I = gamma_I;
params.gamma_D = gamma_D; params.gamma_IH = gamma_IH;
params.theta1 = theta1;  params.delta1  = delta1;
params.delta2 = delta2;  params.N = N;

%% Part b
h = 1; t0 = 0; tf = 200;
tspan = 0:h:tf;
[t, u_sol] = ode45(@(t,u) ebola_rhs(u, params), tspan, u0);

S = u_sol(:,1);
E = u_sol(:,2);
I = u_sol(:,3);
H = u_sol(:,4);
F = u_sol(:,5);
R = u_sol(:,6);

figure
plot(t, I)
xlabel('Time (days)')
ylabel('Infected Population')

figure
plot(t, cumsum(alpha.*E))
xlabel('Time (days)')
ylabel('Cumulative Infections')
ylim([0, 9000])

%% Part d
params.delta2 = 0.1;

[t, u_sol] = ode45(@(t,u) ebola_rhs(u, params), tspan, u0);

S = u_sol(:,1);
E = u_sol(:,2);
I = u_sol(:,3);
H = u_sol(:,4);
F = u_sol(:,5);
R = u_sol(:,6);

figure
plot(t, I)
xlabel('Time (days)')
ylabel('Infected Population')

figure
plot(t, cumsum(alpha.*E))
xlabel('Time (days)')
ylabel('Cumulative Infections')
ylim([0, 9000])

%% Part e
params.delta2 = 0.5;
params.gamma_H = 1/3.0;
params.theta1 = 0.25;

[t, u_sol] = ode45(@(t,u) ebola_rhs(u, params), tspan, u0);

S = u_sol(:,1);
E = u_sol(:,2);
I = u_sol(:,3);
H = u_sol(:,4);
F = u_sol(:,5);
R = u_sol(:,6);

figure
plot(t, I)
xlabel('Time (days)')
ylabel('Infected Population')

figure
plot(t, cumsum(alpha.*E))
xlabel('Time (days)')
ylabel('Cumulative Infections')
ylim([0, 9000])

%% Part f
params.delta2 = 0.1;
params.gamma_H = 1/3.0;
params.theta1 = 0.25;

[t, u_sol] = ode45(@(t,u) ebola_rhs(u, params), tspan, u0);

S = u_sol(:,1);
E = u_sol(:,2);
I = u_sol(:,3);
H = u_sol(:,4);
F = u_sol(:,5);
R = u_sol(:,6);

figure
plot(t, I)
xlabel('Time (days)')
ylabel('Infected Population')

figure
plot(t, cumsum(alpha.*E))
xlabel('Time (days)')
ylabel('Cumulative Infections')
ylim([0, 9000])