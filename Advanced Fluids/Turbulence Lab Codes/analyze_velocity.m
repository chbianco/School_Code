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
set(groot, 'defaultFigureColor', 'white');

%% Load Velocity Data
fprintf('Select file for analysis');
[dataFile, dataPath] = uigetfile;
allData = load(strcat(dataPath, dataFile)); %%Load the file generated from xwire_processdata

Vels = allData.Vels; %Assign the velocity matrix to Vels
Speeds = allData.Speed; %Assign wind speed to Speeds
%% Mean Velocities
%Initialize velocity average vectors
u_ave = zeros(1,3);
v_ave = zeros(1,3);

%Compute average u and v for all three speeds
for i = 1:3
u_ave(i) = mean(Vels(:, 1, i));
v_ave(i) = mean(Vels(:, 2, i));
end
Speeds
u_ave
v_ave

%% Reynolds Stresses

%Initialize vectors
up = zeros(500000, 3);
vp = zeros(500000, 3);
uu = zeros(1, 3);
vv = zeros(1,3);
uv = zeros(1,3);

%Compute Reynolds Stresses all three speeds
for i = 1:3
    %Calculate fluctuating velocity  
    up(:, i) = Vels(:,1,i) - u_ave(i);
    vp(:, i) = Vels(:,2,i) - v_ave(i);

    %Calculate Reynolds Stresses
    uu(i) = mean(up(i).*up(i))./(Speeds(i)^2);
    vv(i) = mean(vp(i).*vp(i))./(Speeds(i)^2);
    uv(i) = mean(up(i).*vp(i))./(Speeds(i)^2);
end

uu
vv
uv

%% Skewness and Flatness

%Initialize vectors
S_up = zeros(1, 3);
F_up = zeros(1,3);
S_vp = zeros(1, 3);
F_vp = zeros(1, 3);

%Compute Skewness and Flatness for all three speeds
for i = 1:3
    %Skewness and Flatness for uprime
    S_up(i) = mean(up(:,i).^3)/(  mean(up(:, i).^2)^(3/2) );
    F_up(i) = mean(up(:,i).^4)/(  mean(up(:, i).^2)^(2) );
    
    %Skewness and Flatness for vprime
    S_vp(i) = mean(vp(:,i).^3)/(  mean(vp(:, i).^2)^(3/2) );
    F_vp(i) = mean(vp(:,i).^4)/(  mean(vp(:, i).^2)^(2) );
end

S_up 
F_up 
S_vp 
F_vp 

%% Temporal Correlation Functions

u_t_corrs = cell(1,3);
v_t_corrs = cell(1,3);
uv_t_corrs = cell(1,3);

u_t_lags = cell(1,3);
v_t_lags = cell(1,3);
uv_t_lags = cell(1,3);


for i = 1:3
    %Get autocorrelation for u
    [u_autocorr, u_lags] = xcorr(up(:,i), 'coeff');
    [v_autocorr, v_lags] = xcorr(vp(:,i), 'coeff');
    [uv_autocorr, uv_lags] = xcorr(up(:,i), vp(:,i), 'coeff');

    % Convert lags into physical time using digitizing rate
    dt = 1/51200; %Sampling interval is 1/digitizing_rate

    u_time_lags = u_lags * dt;
    v_time_lags = v_lags * dt;
    uv_time_lags = v_lags * dt;

    % Plot temporal autocorrelation
    figure(i)

    subplot(1,3,1)
    plot(u_time_lags, u_autocorr, 'LineWidth', 2)
    xlabel('Time lag, $\tau$ [s]')
    ylabel('Temporal Autocorrelation, R($\tau$)')
    grid on
    xlim([-max(u_time_lags)/2, max(u_time_lags)/2])
    title(strcat('u Fluctuation, Speed = ', num2str(Speeds(i)), ' m/s'))

    subplot(1,3,2)
    plot(v_time_lags, v_autocorr, 'LineWidth', 2)
    xlabel('Time lag, $\tau$ [s]')
    ylabel('Temporal Autocorrelation, R($\tau$)')
    grid on
    xlim([-max(v_time_lags)/2, max(v_time_lags)/2])
    title(strcat('v, Speed = ', num2str(Speeds(i)), ' m/s'))

    subplot(1,3,3)
    plot(uv_time_lags, uv_autocorr, 'LineWidth', 2)
    xlabel('Time lag, $\tau$ [s]')
    ylabel('Temporal Autocorrelation, R($\tau$)')
    grid on
    xlim([-max(uv_time_lags)/2, max(uv_time_lags)/2])
    title(strcat('uv, Speed = ', num2str(Speeds(i)), ' m/s'))

    %Store corrs and lags
    u_t_corrs{i} = u_autocorr;
    v_t_corrs{i} = v_autocorr;
    uv_t_corrs{i} = uv_autocorr;

    u_t_lags{i} = u_time_lags;
    v_t_lags{i} = v_time_lags;
    uv_t_lags{i} = uv_time_lags;

end

