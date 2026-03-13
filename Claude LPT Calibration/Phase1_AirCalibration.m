% =========================================================================
% Phase1_AirCalibration.m
% =========================================================================
% PURPOSE:
%   Estimate intrinsic parameters (focal length, principal point, radial
%   and tangential distortion) for each camera using checkerboard images
%   taken IN AIR before the flume is filled.
%
%   Uses MATLAB's Camera Calibration Toolbox functions (detectCheckerboardPoints,
%   estimateCameraParameters) — part of the Computer Vision Toolbox.
%
% OUTPUTS (saved to results/):
%   intrinsics_camN.mat   — cameraIntrinsics object per camera
%   phase1_report.mat     — full calibration structs + reprojection errors
%
% RECOMMENDATION:
%   Collect 40-80 checkerboard images per camera, covering:
%     • All regions of the image (tilt left/right/up/down)
%     • Multiple distances (not just one working distance)
%     • ±30° tilt in pitch and yaw
%   This ensures the distortion model is well-conditioned.
% =========================================================================

clear; clc;
LPT_Config;   % loads cfg

fprintf('=== Phase 1: Air Intrinsic Calibration ===\n\n');

nCams       = cfg.nCams;
resultsDir  = cfg.resultsDir;
squareSize  = cfg.cbSquareSize * 1000;  % detectCheckerboardPoints expects mm
boardSize   = cfg.cbBoardSize;
airDir      = cfg.airCalibDir;

intrinsics  = cell(nCams, 1);
calib       = cell(nCams, 1);

for c = 1:nCams
    fprintf('--- Camera %d ---\n', c);
    camDir = fullfile(airDir, sprintf('cam%d', c));

    if ~exist(camDir, 'dir')
        error('Air calibration image folder not found: %s', camDir);
    end

    % Collect images
    exts   = {'*.png','*.jpg','*.tif','*.bmp'};
    imgFiles = [];
    for e = 1:numel(exts)
        found = dir(fullfile(camDir, exts{e}));
        imgFiles = [imgFiles; found]; %#ok<AGROW>
    end
    if isempty(imgFiles)
        error('No images found in %s', camDir);
    end
    fprintf('  Found %d images\n', numel(imgFiles));

    % Detect checkerboard corners in all images
    imagePoints = [];
    validIdx    = [];
    for i = 1:numel(imgFiles)
        fpath = fullfile(imgFiles(i).folder, imgFiles(i).name);
        img   = imread(fpath);
        if size(img,3) == 3, img = rgb2gray(img); end

        [pts, boardSizeDetected] = detectCheckerboardPoints(img);

        % Validate detected board matches expected size
        if isequal(boardSizeDetected, boardSize)
            imagePoints(:,:,end+1) = pts; %#ok<AGROW>
            validIdx(end+1)        = i;   %#ok<AGROW>
        else
            fprintf('  [WARN] Image %d: board not detected or wrong size, skipping\n', i);
        end
    end

    nValid = size(imagePoints, 3);
    fprintf('  Valid images: %d / %d\n', nValid, numel(imgFiles));
    if nValid < 10
        warning('Camera %d: fewer than 10 valid images — calibration may be unreliable', c);
    end

    % Generate world points for the checkerboard
    worldPoints = generateCheckerboardPoints(boardSize, squareSize);

    % Calibrate — estimate K, distortion (radial k1,k2,k3 + tangential p1,p2)
    imgSize = fliplr(cfg.imageSize);   % [height width] for estimateCameraParameters
    [camParams, imagesUsed, estErrors] = estimateCameraParameters( ...
        imagePoints, worldPoints, ...
        'ImageSize',            imgSize, ...
        'EstimateSkew',         false, ...
        'NumRadialDistortionCoefficients', 3, ...
        'EstimateTangentialDistortion',    true, ...
        'WorldUnits',           'millimeters');

    % Report
    meanErr = mean(camParams.ReprojectionErrors(:));
    fprintf('  Mean reprojection error: %.4f px\n', meanErr);
    if meanErr > 0.5
        warning('Camera %d: mean reprojection error > 0.5 px — check images', c);
    end

    % Store
    intrinsics{c} = camParams.Intrinsics;
    calib{c}      = struct( ...
        'params',       camParams, ...
        'imagesUsed',   imagesUsed, ...
        'errors',       estErrors, ...
        'meanReprojErr',meanErr);

    % Save individual intrinsics
    outFile = fullfile(resultsDir, sprintf('intrinsics_cam%d.mat', c));
    K            = camParams.Intrinsics.K;
    distCoeffs   = [camParams.RadialDistortion, camParams.TangentialDistortion];
    imageSize    = cfg.imageSize;
    save(outFile, 'K', 'distCoeffs', 'imageSize', 'camParams');
    fprintf('  Saved: %s\n', outFile);

    % Visualise reprojection errors
    fig = figure('Visible','off','Name', sprintf('Cam %d Reprojection', c));
    showReprojectionErrors(camParams);
    title(sprintf('Camera %d — Mean error = %.3f px', c, meanErr));
    saveas(fig, fullfile(resultsDir, sprintf('cam%d_reprojection.png', c)));
    close(fig);

    % Visualise extrinsic poses used during air calibration
    fig2 = figure('Visible','off','Name', sprintf('Cam %d Poses', c));
    showExtrinsics(camParams, 'PatternCentric');
    title(sprintf('Camera %d — Calibration Board Poses', c));
    saveas(fig2, fullfile(resultsDir, sprintf('cam%d_airPoses.png', c)));
    close(fig2);
end

% Save combined report
save(fullfile(resultsDir, 'phase1_report.mat'), 'intrinsics', 'calib', 'cfg');
fprintf('\n[Phase 1] Complete. Results saved to %s\n', resultsDir);

% Print summary table
fprintf('\n%s\n', repmat('=',1,55));
fprintf('%-10s %-20s %-20s\n','Camera','Focal Length (px)','Mean Reproj Err (px)');
fprintf('%s\n', repmat('-',1,55));
for c = 1:nCams
    fl = intrinsics{c}.FocalLength;
    fprintf('  Cam %d     [%.1f  %.1f]           %.4f\n', ...
        c, fl(1), fl(2), calib{c}.meanReprojErr);
end
fprintf('%s\n', repmat('=',1,55));
