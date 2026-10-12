%% Preamble
close all;
clear variables;
clc;

set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextFontSize', 20);
set(groot, 'defaultAxesFontSize', 24);
set(groot, 'defaultLineLineWidth', 2);

%% Load Data
data = readmatrix('Turbulence_HW3_data.txt');

%% Zero out uw and vw

data(:, 13) = zeros(1, 541);
data(:, 14) = zeros(1, 541);

%% Deviatoric stress
y  = data(:, 1);
uu = data(:, 4).^2;
vv = data(:, 5).^2;
ww = data(:, 6).^2;
uv = data(:, 11);
uw = data(:, 12);
vw = data(:, 13);

N     = size(data, 1);
theta = zeros(N, 1);
for i = 1:N
    R = [uu(i) uv(i) uw(i);
         uv(i) vv(i) vw(i);
         uw(i) vw(i) ww(i)];
    k = 0.5*trace(R); % TKE
    a = R - (2/3)*k*eye(3); % anisotropic (deviatoric) part
    [V, D] = eig(a);
    [~, idx] = max(diag(D));  % major principal axis
    e1 = V(:, idx);
    theta(i) = atand(e1(2)/e1(1)); % angle w.r.t. x-axis
end

%Plot principal direction
figure;
plot(y, theta, 'k');
xlabel('$y/h$');
ylabel('$\theta$ ($^\circ$)');
grid on;

%% Eddy viscosity
yp  = data(:, 2);             % y+
Up  = data(:, 3);             % <u>+
uvp = data(:, 11);            % <u'v'>+

% dU+/dy+ : 2nd-order FD on non-uniform grid
dUdy = zeros(N, 1);
for i = 2:N-1
    h1 = yp(i) - yp(i-1);
    h2 = yp(i+1) - yp(i);
    dUdy(i) = -h2/(h1*(h1+h2))*Up(i-1) + (h2-h1)/(h1*h2)*Up(i) ...
              + h1/(h2*(h1+h2))*Up(i+1);
end
h1 = yp(2) - yp(1);           % one-sided forward (wall)
h2 = yp(3) - yp(2);
dUdy(1) = -(2*h1+h2)/(h1*(h1+h2))*Up(1) + (h1+h2)/(h1*h2)*Up(2) ...
          - h1/(h2*(h1+h2))*Up(3);
h1 = yp(N-1) - yp(N-2);       % one-sided backward (centerline)
h2 = yp(N) - yp(N-1);
dUdy(N) = h2/(h1*(h1+h2))*Up(N-2) - (h1+h2)/(h1*h2)*Up(N-1) ...
          + (h1+2*h2)/(h2*(h1+h2))*Up(N);



%Plot eddy viscosity
figure;
plot(y, nut/Re_tau, 'k');
xlabel('$y/h$');
ylabel('$\nu_t/(u_\tau h)$');
grid on;