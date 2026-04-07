%% Preamble
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

%% Load Data
tracks = load('jhtdb_long_LPT.mat').LPT;

t = tracks.t;
x = tracks.x;
y = tracks.y;
z = tracks.z; 

%% Test with one track
%Positions
x1 = rmmissing(x(:,1));
y1 = rmmissing(y(:,1));
z1 = rmmissing(z(:,1));

trk_lng = max([length(x1),length(y1), length(z1)]);
ux1 = zeros(1, trk_lng -1);
uy1 = zeros(1, trk_lng -1);
uz1 = zeros(1, trk_lng -1);
t_trk = t(1:trk_lng)';

%Velocities
for j = 1: trk_lng - 1
    dt = t_trk(j+1) - t_trk(j);
    ux1(j) = (x1(j+1) - x1(j))/dt;
    uy1(j) = (y1(j+1) - y1(j))/dt;
    uz1(j) = (z1(j+1) - z1(j))/dt;
end

speed = sqrt(ux1.^2 + uy1.^2 + uz1.^2)';


figure
plot3(x1, y1, z1)
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')
xlim([0 8*pi])
ylim([-1 1])
zlim([0 3*pi])

figure
hold on
plot(t_trk(1:end-1), ux1)
plot(t_trk(1:end-1), uy1)
plot(t_trk(1:end-1), uz1)
xlabel('t')
ylabel('U')
legend({'u', 'v','w'})
hold off

figure
scatter3(x1(1:end-1),y1(1:end-1),z1(1:end-1), 20, speed, 'filled')
c = colorbar;
c.Label.String = 'Speed';
c.Label.FontName = 'Latex';
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')
xlim([0 8*pi])
ylim([-1 1])
zlim([0 3*pi])

%% For n tracks
n = 100;
speed = NaN(4000, n);

u_vals = NaN (4000, n); 
v_vals = NaN (4000, n);
w_vals = NaN (4000, n);

figure(1)
hold on


for k = 1:n
    xn = rmmissing(x(:,k));
    yn = rmmissing(y(:,k));
    zn = rmmissing(z(:,k));

    trk_lng = min([length(xn),length(yn), length(zn)]);
    un = zeros(1, trk_lng -1);
    vn = zeros(1, trk_lng -1);
    wn = zeros(1, trk_lng -1);
    t_trk = t(1:trk_lng)';
    
    %Velocities
    for j = 1: trk_lng - 1
        dt = t_trk(j+1) - t_trk(j);
        un(j) = (xn(j+1) - xn(j))/dt;
        vn(j) = (yn(j+1) - yn(j))/dt;
        wn(j) = (zn(j+1) - zn(j))/dt;
    end
    u_vals(1:trk_lng-1, k) = un';
    v_vals(1:trk_lng-1, k) = vn';
    w_vals(1:trk_lng-1, k) = wn';
    
    speed(1:trk_lng-1, k) = sqrt(un.^2 + vn.^2 + wn.^2)';

    scatter3(xn(1:end-1),yn(1:end-1),zn(1:end-1), 20, speed(1:trk_lng-1,k), 'filled')


end 

c = colorbar;
c.Label.String = 'Speed';
c.Label.FontName = 'Latex';
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')
xlim([0 8*pi])
ylim([-1 1])
zlim([0 3*pi])

view(3)

hold off