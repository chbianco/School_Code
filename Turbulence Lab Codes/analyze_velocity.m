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
    xlim([0, 0.01])
    title(strcat('u, Speed = ', num2str(Speeds(i)), ' m/s'))

    subplot(1,3,2)
    plot(v_time_lags, v_autocorr, 'LineWidth', 2)
    xlabel('Time lag, $\tau$ [s]')
    ylabel('Temporal Autocorrelation, R($\tau$)')
    grid on
    xlim([0, 0.01])
    title(strcat('v, Speed = ', num2str(Speeds(i)), ' m/s'))

    subplot(1,3,3)
    plot(uv_time_lags, uv_autocorr, 'LineWidth', 2)
    xlabel('Time lag, $\tau$ [s]')
    ylabel('Temporal Autocorrelation, R($\tau$)')
    grid on
    xlim([0, 0.01])
    title(strcat('uv, Speed = ', num2str(Speeds(i)), ' m/s'))

    %Store corrs and lags
    u_t_corrs{i} = u_autocorr;
    v_t_corrs{i} = v_autocorr;
    uv_t_corrs{i} = uv_autocorr;

    u_t_lags{i} = u_time_lags;
    v_t_lags{i} = v_time_lags;
    uv_t_lags{i} = uv_time_lags;

end

%% Spatial Correlation Functions

u_x_lags = cell(1,3);
v_x_lags = cell(1,3);
uv_x_lags = cell(1,3);

for i = 1:3
    U = Speeds(i);
    
    %Change lags to spatial by multiplying by veloccity
    u_pos_lags = u_t_lags{i} .* U;
    v_pos_lags = v_t_lags{i} .* U;
    uv_pos_lags = uv_t_lags{i} .* U;
    
    % Plot spatial autocorrelation
    figure(i)

    subplot(1,3,1)
    plot(u_pos_lags, u_autocorr, 'LineWidth', 2)
    xlabel('Spatial lag [m]')
    ylabel('Spatial Autocorrelation')
    grid on
    xlim([0, 0.2])
    title(strcat('u, Speed = ', num2str(Speeds(i)), ' m/s'))

    subplot(1,3,2)
    plot(v_pos_lags, v_autocorr, 'LineWidth', 2)
    xlabel('Spatial lag [m]')
    ylabel('Spatial Autocorrelation')
    grid on
    xlim([0, 0.2])
    title(strcat('v, Speed = ', num2str(Speeds(i)), ' m/s'))

    subplot(1,3,3)
    plot(uv_pos_lags, uv_autocorr, 'LineWidth', 2)
    xlabel('Spatial lag [m]')
    ylabel('Spatial Autocorrelation')
    grid on
    xlim([0, 0.2])
    title(strcat('uv, Speed = ', num2str(Speeds(i)), ' m/s'))

    %Store spatial lags
    u_x_lags{i} = u_pos_lags;
    v_x_lags{i} = v_pos_lags;
    uv_x_lags{i} = uv_pos_lags;

end


%% Taylor and Integral Length Scales, Kolmogorov timescales, Reynolds Numbers
%Initialize vectors
Taylor_t = zeros(1,3); % Taylor time scale
Taylor_x = zeros(1,3);
integral_t = zeros(1,3);
integral_x = zeros(1,3);
eps = zeros(1,3);
nK = zeros(1,3);
tK = zeros(1,3);
Re_macro = zeros(1,3);
Re_turb = zeros(1,3);

nu = 1.48*10^(-5); %Kinematic viscosity of air
D = 0.01; %Diameter of the cylinder, meters

for i = 1:3
    %Find parabola describing the curvature of the correlation at zero 
    zero_lag_idx = find(u_t_lags{i} == 0);
    fit_range = (zero_lag_idx-10):(zero_lag_idx+10); %This is kind of arbitrary 
    p = polyfit(u_t_lags{i}(fit_range), u_t_corrs{i}(fit_range), 2); %Fit the points to a parabola

    %Compute Taylor timescale and length scale 
    Taylor_t(i) = sqrt(-1/(2*p(1)));% Taylor time scale MIGHT BE SUS BUT ORDER OF MAGNITUDE IS SOLID
    Taylor_x(i) = Taylor_t(i) * Speeds(i);
    
    %Compute integral length scale
    zero_cross =  find(u_t_corrs{i}(zero_lag_idx:end) <= 0, 1, 'first') + zero_lag_idx - 1; %Get first zero crossing 
    integral_t(i) = trapz(u_t_lags{i}(zero_lag_idx:zero_cross), u_t_corrs{i}(zero_lag_idx:zero_cross)); %Integrate Ruu to find timescale
    integral_x(i) = integral_t(i) * Speeds(i);

    %Estimate dissipation rate, eps
    eps(i) = 15 * nu * mean(up(:,i).^2) ./ Taylor_x(i); %MIGHT NEED TO DIVIDE BY TAYLOR_t INSTEAD

    %Find Kolmogorov length and time scales
    nK(i) = (nu^3 / eps(i))^(1/4);
    tK(i) = sqrt(nu/eps(i));
    
    %Reynolds Numbers
    Re_macro(i) = Speeds(i)*D/nu;
    Re_turb(i) = rms(up(:,i))*Taylor_x(i)/nu;
end

%% Longitudinal and transverse functions
f_vec = cell(1,3);
g_vec = cell(1,3);

for i = 1:3
%Compute f and g
f = u_t_corrs{i};
g = v_t_corrs{i};
r = u_x_lags{i}';

%Compute predicted g from the isotropic relations
dfdr = gradient(f, r);
g_pred = f + 0.5.*r.*dfdr;

%Plot g and g predicted
figure;
plot(r, g_pred, 'r--', r, g, 'b-', 'LineWidth',1.5)
xlabel('Spatial lag r [m]')
ylabel('Transverse correlation g(r)')
legend('measured g','predicted g','Location','Best')
grid on
title('Isotropy test: measured vs.\ predicted transverse correlation')

end

%% Spatial Power Spectra

%initialize values for the spectrum
nfft = 2^14;
window = hamming(nfft);
nover = nfft/2;

digiRate = 51200;
dt = 1/digiRate;
colors = lines(3);

figure(1); hold on
figure(2); hold on
figure(3); hold on

for i = 1:3
    U = Speeds(i);

    %Temporal PSD
    [Suu, freq] = pwelch(up(:, i), window, nover, nfft, digiRate, 'one-sided');

    % Map to spatial k with Taylor hypothesis
    k = 2*pi * freq / U; % [rad/m]
    E  = U * Suu; % spatial spectrum [ (m^3/s^2) * m ]

    % Plot raw spatial spectra
    figure(1)
    subplot(1, 3, i)
    loglog(k, E, 'Color', colors(i,:));
    title(sprintf('U=%.1f m/s',U))
    xlim([0, 10^5])
    ylim([10^(-13), 10^(-2)])
    xlabel('k')
    ylabel('Spatial Spectrum')

    %Scale spectra by inertia/Kolmogorov
    E_nk = E /((eps(i)*nu^5)^1/4);
    k_nk = k*nK(i);

    %Plot inertial scaled PSD
    figure(2)
    subplot(1, 3, i)
    loglog(k_nk, E_nk, 'Color', colors(i,:));
    title(sprintf('U=%.1f m/s',U))
    %xlim([0, 10^5])
    %ylim([10^(-13), 10^(-2)])
    xlabel('$\eta_k$k')
    ylabel('$E(k)/(\epsilon\nu^5)^{1/4}$')

    %Scale spectra by outer
    E_D = E / (Speeds(i)*D);
    k_D = k*D;

    %Plot outer scaled PSD
    figure(3)
    subplot(1, 3, i)
    loglog(k_D, E_D, 'Color', colors(i,:));
    title(sprintf('U=%.1f m/s',U))
    %xlim([0, 10^5])
    %ylim([10^(-13), 10^(-2)])
    xlabel('Dk')
    ylabel('$E(k)/(UD)$')
    

end
