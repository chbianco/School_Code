%% PREAMBLE
close all;
clear variables;
clc;

% Set default formatting properties for figures using LaTeX formatting
set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultFigureColor', 'white');

%% Number of runs, color, legend
n = 2;
color_vec = ["#0072BD", "#D95319", "#EDB120"];  % String array
legend_vec = ["High Flex", "Low Flex"]; % String array

%% Create Figure and Configure Subplots
% Make the figure much larger by setting its units to normalized and outerposition to full screen.
figure('units','normalized','outerposition',[0 0 1 1]);

% Subplot for normalized x-position vs. time
ax1 = subplot(1, 3, 1);
hold(ax1, 'on');
xlabel(ax1, '$t_c = \frac{tU}{L}$');
ylabel(ax1, '$x^* = x/L$');
title(ax1, 'x Position');
xlim(ax1, [10 40]);
ylim(ax1, [-0.4 0.4]);
grid(ax1, 'on');

% Subplot for normalized y-position vs. time
ax2 = subplot(1, 3, 2);
hold(ax2, 'on');
xlabel(ax2, '$t_c = \frac{tU}{L}$');
ylabel(ax2, '$y^* = y/L$');
title(ax2, 'y Position');
xlim(ax2, [10 40]);
ylim(ax2, [-0.4 0.4]);
grid(ax2, 'on');

% Subplot for distance from the starting position vs. time
ax3 = subplot(1, 3, 3);
hold(ax3, 'on');
xlabel(ax3, '$t_c = \frac{tU}{L}$');
ylabel(ax3, '$r^* = \sqrt{x^2+y^2}/L$');
title(ax3, 'Position from Start');
xlim(ax3, [10 40]);
ylim(ax3, [0 0.6]);
grid(ax3, 'on');

%% Loop Over Runs and Plot Data from Each File
for i = 1:n
    % Select a MAT file for the current run
    [file, path] = uigetfile('*.mat', sprintf('Select MAT file for run %d', i));
    if isequal(file, 0)
        disp('No file selected. Skipping this run.');
        continue;
    end
    fullFileName = fullfile(path, file);

    % Extract the video delay from the filename
    timeStart = split(file, '=');
    timeStart = split(timeStart(2), '.');
    timeStart = strcat('0.', timeStart(1));
    timeStart = str2double(timeStart);
    
    % Load the selected file and check if it contains 'trackingData'
    data = load(fullFileName);
    if isfield(data, 'trackingData')
        trackingData = data.trackingData;
    else
        warning('The variable ''trackingData'' was not found in file: %s. Skipping this run.', file);
        continue;
    end

    %% Extract Data and Compute Non-dimensional Quantities
    % Non-dimensionalize time: multiply by (speed)/(chord) and add offset from file delay and starting time scaling
    time = (trackingData(:, 1) + timeStart) .* 0.35 ./ 0.1 + 10; 

    % Extract x and y positions
    x = trackingData(:, 2);
    y = trackingData(:, 3);

    % Normalize positions by subtracting the starting values and non-dimensionalize by chord (in cm)
    x0 = x(1);
    y0 = y(1);
    normX = (x - x0) / 10;
    normY = (y - y0) / 10;

    % Compute current non-dimensional distance from the starting position
    normDistance = sqrt(normX.^2 + normY.^2);

    %% Plot Data for Current Run on the Same Axes
    % Add a run-specific label using DisplayName so legends can distinguish runs.
    plot(ax1, time, normX, 'DisplayName', legend_vec(i), 'Color', color_vec(i));
    plot(ax2, time, normY, 'DisplayName', legend_vec(i), 'Color', color_vec(i));
    plot(ax3, time, normDistance, 'DisplayName', legend_vec(i), 'Color', color_vec(i));
end

% Add legends to each subplot to differentiate between runs
legend(ax1, 'show');
legend(ax2, 'show');
legend(ax3, 'show');

%% Save Graph
filename = input('Enter filename (without extension): ', 's');
saveas(gcf, [filename, '.png']); % Save the figure as a PNG file
saveas(gcf, [filename, '.fig']); % Save the figure as a FIG file