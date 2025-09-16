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

%% Loading Data 
%REMEMBER TO CHANGE THIS IF ON WINDOWS cause apples hates me 
fileDir = '/Users/christopherbianco/Desktop/School_Code/Wind Physics/HW1';

data = load(fullfile(fileDir, '08_28_2019_22_00_00_000.mat'));

%Convert UTC datenum to datetime. Don't think we need but keeping just in case
% thyme = time_UTC;
% time = NaT(size(thyme));
% for ti = 1 : length(thyme)
%     time(ti) = datetime(floor(thyme(ti)), 'ConvertFrom', 'datenum') ...
%         + days(rem(thyme(ti),1)); % convert datenum to datetime object
% end
% time.TimeZone = 'Etc/UTC'; % add TimeZone field (UTC time)
% time.TimeZone = 'America/Denver'; % shift to NREL time zone



%% Part a: Horizontal wind speed profile
%Initialize measurement height
sonic_heights = [15,41,61, 74, 100, 119];
U_av_sonic = NaN(length(sonic_heights));
U_std_sonic = NaN(length(sonic_heights));

cup_heights = [3, 10, 30, 38, 55, 80, 87, 105, 122, 130];
U_av_cup = NaN(length(cup_heights));
U_std_cup = NaN(length(cup_heights)); 

%Calculate horizontal wind speed from sonic anemometers
for i = 1:length(sonic_heights)
    sonic_height = sonic_heights(i);
    U = sqrt(data.(strcat('Sonic_u_',num2str(sonic_height),'m')).val.^2 + data.(strcat('Sonic_v_',num2str(sonic_height),'m')).val.^2);
    
    U_av_sonic(i) = mean(U);
    U_std_sonic(i) = std(U);
end

%Extract horizontal wind speed from cup anemometers
for i = 1:length(cup_heights)
    cup_height = cup_heights(i);
    
    %Now, we need to take into account the naming convention
    fn = fieldnames(data);
    matchIdx = contains(fn, 'Cup_WS_') & contains(fn, strcat(num2str(cup_height),'m'));
    match = fn(matchIdx);
    
    %Extract U
    U = data.(match{1}).val;
    
    %Do mean and standard deviation
    U_av_cup(i) = mean(U);
    U_std_cup(i) = std(U);
end

%Make the figure
figure(1); hold on; 
xlabel('Horizontal Wind Speed (m/s)');
ylabel('z (m)');
grid on 

%Plot sonic
e1 = errorbar(U_av_sonic, sonic_heights, U_std_sonic, 'horizontal', '*c', 'LineWidth',2);
%Plot cup
e2 = errorbar(U_av_cup, cup_heights, U_std_cup, 'horizontal', '*k', 'LineWidth',2);

% Dummy plots for legend only
h1 = plot(nan, nan, '*c', 'LineWidth', 2);
h2 = plot(nan, nan, '*k', 'LineWidth', 2);

% Legend
legend([h1 h2], {'Sonic','Cup'}, 'Location', 'best')

hold off

%% Part b: Temperature readings
%Initialize measurement height
sonic_heights = [15,41,61,74,100,119];
temp_av_sonic = NaN(length(sonic_heights));
temp_std_sonic = NaN(length(sonic_heights));

solo_heights = [3, 38, 87];
temp_av_solo = NaN(length(solo_heights));
temp_std_solo = NaN(length(solo_heights));

%Extract temp from sonic anemometers
for i = 1:length(sonic_heights)
    sonic_height = sonic_heights(i);
    temp = data.(strcat('Sonic_Temp_clean_', num2str(sonic_height), 'm')).val;
    
    temp_av_sonic(i) = mean(temp);
    temp_std_sonic(i) = std(temp);
end

%Extract temp from stand alone measurements
for i = 1:length(solo_heights)
    solo_height = solo_heights(i);
    temp = data.(strcat('Air_Temp_', num2str(solo_height), 'm')).val;
    
    temp_av_solo(i) = mean(temp);
    temp_std_solo(i) = std(temp);
end

%Make figure
figure(2); hold on;
xlabel('Temperature (C)');
ylabel('z (m)');
grid on

%Plot sonic
e1 = errorbar(temp_av_sonic, sonic_heights, temp_std_sonic, 'horizontal', '*c', 'LineWidth',2);
%Plot cup
e2 = errorbar(temp_av_solo, solo_heights, temp_std_solo, 'horizontal', '*b', 'LineWidth',2);

% Dummy plots for legend only
h1 = plot(nan, nan, '*c', 'LineWidth', 2);
h2 = plot(nan, nan, '*k', 'LineWidth', 2);

% Legend
legend([h1 h2], {'Sonic','Stand Alone Measurement'}, 'Location', 'southwest')

hold off

%Based on these graphs, the temperature readings from the sonic anemometers
%at z = 41m and 74m are quite low, and deviate wildly from all other
%measurements. For that reason, these points look unreliable. 

%% Part c: Adiabatic lapse rate 



