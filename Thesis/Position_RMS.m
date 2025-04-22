clear 
close all
%% Select a MAT file for the current run
    [file, path] = uigetfile;
    fullFileName = fullfile(path, file);

    % Extract the video delay from the filename
    timeStart = split(file, '=');
    timeStart = split(timeStart(2), '.');
    timeStart = strcat('0.', timeStart(1));
    timeStart = str2double(timeStart);
    
    % Load the selected file and check if it contains 'trackingData'
    data = load(fullFileName);
    
    trackingData = data.trackingData;

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
%% RMS
% Example input (replace this with your actual data)
signal = normDistance(1:320);  % Replace with your variable
t = time(1:320);         % Optional, if you have a time vector

% signal(isnan(signal)) = 0.3;

% Step 2: Find the peak region
[~, peak_idx] = max(signal);  % Index of the global maximum

% Step 3: Define a threshold (e.g., % of max) to capture peak region
threshold = 0.025 * signal(peak_idx);  % 70% of peak value
peak_mask = signal > threshold;

% Step 4: Group the continuous regions above the threshold
cc = bwconncomp(peak_mask);  % Find connected regions

% Step 5: Identify the connected region that contains the peak
for k = 1:cc.NumObjects
    if any(cc.PixelIdxList{k} == peak_idx)
        peak_region = cc.PixelIdxList{k};
        break;
    end
end

% Step 6: Separate the peak and non-peak data
signal_peak = signal(peak_region);
nonpeak_region = setdiff(1:length(signal), peak_region);
signal_nonpeak = signal(nonpeak_region);

% Step 7: Compute RMS values
rms_peak = sqrt(mean(signal_peak.^2));
rms_nonpeak = sqrt(mean(signal_nonpeak.^2));

% Output
fprintf('RMS of peak region: %.4f\n', rms_peak);
fprintf('RMS of non-peak region: %.4f\n', rms_nonpeak);


%% Plotting
figure;
hold on;
plot(t, signal, 'k--', 'DisplayName', 'Original Signal');
plot(t(peak_region), signal(peak_region), 'r', 'LineWidth', 2, 'DisplayName', 'Peak Region');
plot(t(nonpeak_region), signal(nonpeak_region), 'b', 'DisplayName', 'Non-Peak Region');
xlabel('Time');
ylabel('Signal');
legend('Location', 'best');
title(sprintf('RMS Peak = %.4f | RMS Non-Peak = %.4f', rms_peak, rms_nonpeak));
grid on;