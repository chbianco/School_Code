% =========================================================================
% ThresholdTester.m
% =========================================================================
% PURPOSE: Run a diagnostic sweep over different EdgeThreshold values 
%          to balance imfindcircles processing speed and detection accuracy.
% =========================================================================
clear; clc; close all;

% --- Configuration ---
testImagePath = 'wand_images/cam1/frame_00172.tif'; % Path to your uploaded test frame
radiusRange   = [20, 90];
methodToTest  = 'TwoStage';  %'TwoStage' or 'PhaseCode' 

% Test a sweep from very permissive (0.03) to aggressive (0.25)
thresholdsToTest = [ 0.08, 0.1, 0.125 0.15];

% --- Image Preparation ---
if ~exist(testImagePath, 'file')
    error('Could not find test image: %s. Please update the path.', testImagePath);
end

im = imread(testImagePath);
if size(im, 3) == 3
    im_gray = rgb2gray(im); % Convert 1920x1088 RGB to grayscale
else
    im_gray = im;
end

% Boost visualization contrast slightly since your raw frame is naturally dark
im_disp = imadjust(im_gray);

figure('Name', 'EdgeThreshold Diagnostic Sweep', 'Position', [100, 100, 1200, 800]);

% --- Run Sweep ---
for i = 1:numel(thresholdsToTest)
    t = thresholdsToTest(i);
    
    % Time the execution of each threshold pass
    tic;
    [centers, radii, metrics] = imfindcircles(im_gray, radiusRange, ...
        'ObjectPolarity', 'dark', ...
        'Method', methodToTest, ...
        'EdgeThreshold', t);
    elapsedTime = toc;
    
    % Plot the results
    ax = subplot(2, 2, i);
    imshow(im_disp, 'Parent', ax);
    hold(ax, 'on');
    
    numFound = size(centers, 1);
    
    if numFound >= 2
        % Draw detected circles (Green = Top 2 strongest, Red = Extra noise circles)
        for c = 1:numFound
            if c <= 2
                viscircles(ax, centers(c,:), radii(c), 'EdgeColor', 'g', 'LineWidth', 2);
                plot(ax, centers(c,1), centers(c,2), 'gx', 'MarkerSize', 10, 'LineWidth', 2);
            else
                viscircles(ax, centers(c,:), radii(c), 'EdgeColor', 'r', 'LineWidth', 1, 'LineStyle', '--');
            end
        end
        titleStr = sprintf('Threshold: %.2f | Time: %.4fs\nFound %d circles (Top 2 in Green)', ...
            t, elapsedTime, numFound);
    else
        titleStr = sprintf('Threshold: %.2f | Time: %.4fs\nFAILED: Found only %d circle(s)', ...
            t, elapsedTime, numFound);
    end
    
    title(ax, titleStr, 'FontSize', 10, 'FontWeight', 'bold');
    hold(ax, 'off');
end

sgtitle(sprintf('imfindcircles Method: %s | Range: [%d, %d]', ...
    methodToTest, radiusRange(1), radiusRange(2)), 'FontSize', 14);