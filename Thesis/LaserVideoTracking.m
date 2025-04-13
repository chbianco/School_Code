%% Laser Dot Tracker Script with Calibration, Plotting, and Data Saving
% This script reads a video file, brightens each frame, and then performs:
% 1) Calibration: The user selects two points known to be 3 cm apart.
%    The pixel distance between these points is used to compute a scale factor.
% 2) Laser Dot Selection: The user selects the laser dot on the first frame.
% 3) Tracking: The script tracks the laser dot in subsequent frames using an
%    ROI around the last known dot location.
% 4) Plotting: It plots the x and y positions (in cm) versus time.
% 5) Saving: The final results (time, x, and y in cm) are saved to a MAT file.

clear; clc; close all;

%% PARAMETERS
[filename, pathname] = uigetfile({ ...
    '*.mp4;*.avi;*.mov;*.mkv','Video Files (*.mp4, *.avi, *.mov, *.mkv)'; ...
    '*.*','All Files (*.*)'}, ...
    'Select a video file');
if isequal(filename,0)
    error('No video file selected.');
end
videoFile = fullfile(pathname, filename);

thresholdLevel = 0.8;          % Normalized threshold for binarizing (0 to 1)
ROI_radius = 50;               % ROI radius (in pixels)

%% READ VIDEO & INITIAL FRAME
vidObj = VideoReader(videoFile);

% Read and brighten the first frame for calibration and laser dot selection.
firstFrame = readFrame(vidObj);
brightFirstFrame = imadjust(firstFrame, stretchlim(firstFrame), []);

%% CALIBRATION STEP: Select two points 10 cm apart
figure;
imshow(brightFirstFrame);
title('Calibration: Click TWO points that are 10 cm apart, then press Enter');
% Get separate x and y outputs.
[xCal, yCal] = getpts;
if length(xCal) < 2 || length(yCal) < 2
    error('Calibration: Please select TWO points for calibration.');
end
% Combine into two calibration points.
calibPt1 = [xCal(1), yCal(1)];
calibPt2 = [xCal(2), yCal(2)];
pixelDistance = sqrt((calibPt1(1) - calibPt2(1))^2 + (calibPt1(2) - calibPt2(2))^2);
scaleFactor = 10 / pixelDistance;  % cm per pixel conversion factor
close;

%% LASER DOT SELECTION: Select the laser dot manually on the first frame
figure;
imshow(brightFirstFrame);
title('Select the Laser Dot (click once then press Enter)');
[xInit, yInit] = ginput(1);
laserDot = [xInit, yInit];  % Laser dot initial position in pixels
close;

% Set the last known position using the manual selection.
lastDotPos = laserDot;

%% PREPARE FOR TRACKING
% Initialize tracking data storage: each row is [time (s), x (cm), y (cm)]
trackingData = [];

%% PROCESS EACH FRAME
figure;  % For real-time visualization

while hasFrame(vidObj)
    % Read the next frame and enhance its brightness.
    frame = readFrame(vidObj);
    brightFrame = imadjust(frame, stretchlim(frame), []);
    
    % Convert the frame to grayscale if needed.
    if size(brightFrame, 3) == 3
        grayFrame = rgb2gray(brightFrame);
    else
        grayFrame = brightFrame;
    end
    
    % Apply a median filter to reduce noise.
    grayFrameFiltered = medfilt2(grayFrame, [3 3]);
    
    % Define the ROI based on the last known laser dot position.
    [rows, cols] = size(grayFrameFiltered);
    xCenter = round(lastDotPos(1));
    yCenter = round(lastDotPos(2));
    xMin = max(xCenter - ROI_radius, 1);
    xMax = min(xCenter + ROI_radius, cols);
    yMin = max(yCenter - ROI_radius, 1);
    yMax = min(yCenter + ROI_radius, rows);
    
    % Extract the ROI from the filtered frame.
    roiFrame = grayFrameFiltered(yMin:yMax, xMin:xMax);
    
    % Binarize the ROI using the specified threshold.
    bwROI = imbinarize(roiFrame, thresholdLevel);
    
    % Label connected components (blobs) in the ROI.
    cc = bwconncomp(bwROI);
    stats = regionprops(cc, roiFrame, 'Centroid', 'MaxIntensity');
    
    % Initialize current laser dot position.
    currentDotPos = [];
    
    if ~isempty(stats)
        % Choose the blob with maximum intensity as the laser dot.
        maxIntensities = [stats.MaxIntensity];
        [~, idx] = max(maxIntensities);
        roiCentroid = stats(idx).Centroid;
        % Convert ROI centroid to full image coordinates.
        currentDotPos = [roiCentroid(1) + xMin - 1, roiCentroid(2) + yMin - 1];
        lastDotPos = currentDotPos;
    end
    
    % Get the current time in seconds.
    currentTime = vidObj.CurrentTime;
    
    % Save the current tracking result, converting pixel coordinates to cm.
    if isempty(currentDotPos)
        trackingData = [trackingData; currentTime, NaN, NaN];
    else
        trackingData = [trackingData; currentTime, currentDotPos(1) * scaleFactor, currentDotPos(2) * scaleFactor];
    end
    
    % Real-time visualization.
    imshow(brightFrame);
    hold on;
    rectangle('Position', [xMin, yMin, (xMax - xMin), (yMax - yMin)], 'EdgeColor', 'g');
    if ~isempty(currentDotPos)
        plot(currentDotPos(1), currentDotPos(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    end
    hold off;
    title(sprintf('Time: %.2f s', currentTime));
    drawnow;
end

%% FINAL PLOTTING: X and Y Positions vs Time
timeArray = trackingData(:,1);
xPositions = trackingData(:,2);
yPositions = trackingData(:,3);

figure;
subplot(2,1,1);
plot(timeArray, xPositions, '-o');
xlabel('Time (s)');
ylabel('X Position (cm)');
title('X Position vs Time');

subplot(2,1,2);
plot(timeArray, yPositions, '-o');
xlabel('Time (s)');
ylabel('Y Position (cm)');
title('Y Position vs Time');

%% SAVE FINAL TRACKING RESULTS
dataName = split(filename, '.');
save(strcat(dataName{1},'_time=', dataName{2},'.mat'), 'trackingData');
disp('Laser dot tracking completed. Tracking data saved in finalTrackingResults.mat');