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

%% Load Data
% Ask user to select a .mat file for analysis
[file, path] = uigetfile('*.mat', 'Select a MAT file for analysis');
if isequal(file, 0)
    disp('No file selected. Exiting script.');
    return;
end
fullFileName = fullfile(path, file);

% Load the selected file
data = load(fullFileName);
if isfield(data, 'trackingData')
    trackingData = data.trackingData;
else
    error('The variable ''trackingData'' was not found in the selected file.');
end

%% Extract data and analyze

% Extract time, x, and y from trackingData
time = trackingData(:, 1).*0.35./0.1; %Non-dim by speed in m/s and chord in m
x = trackingData(:, 2);
y = trackingData(:, 3);

% Normalize positions by subtracting the starting values
x0 = x(1);
y0 = y(1);
normX = (x - x0)/10; %Non-dim by chord in cm
normY = (y - y0)/10;

% Calculate the distance from the starting position for each time point
normDistance = sqrt(normX.^2 + normY.^2);

%% Plot
% Create subplots for the three graphs
figure;

% Subplot 1: Normalized X-position vs. Time
subplot(1, 3, 1);
plot(time, normX, '-b', 'LineWidth', 1.5);
xlabel('$t_c = \frac{tU}{L}$');
ylabel('$x^* = x/L$');
ylim([-0.4 .4])
title('x Position');
grid on;

% Subplot 2: Normalized Y-position vs. Time
subplot(1, 3, 2);
plot(time, normY, '-r', 'LineWidth', 1.5);
xlabel('$t_c = \frac{tU}{L}$');
ylim([-0.4 .4])
ylabel('$y^* = y/L$');
title('y Position');
grid on;

% Subplot 3: Distance from Starting Position vs. Time
subplot(1, 3, 3);
plot(time, normDistance, '-k', 'LineWidth', 1.5);
xlabel('$t_c = \frac{tU}{L}$');
ylabel('$r^* = \sqrt{x^2 + y^2}/L$');
ylim([0 0.6])
title('Position from Start');
grid on;
