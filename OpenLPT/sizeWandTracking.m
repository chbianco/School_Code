% =========================================================================
% WandTracking.m
% =========================================================================
% PURPOSE:
%   Automatically detect the two LED tips of the calibration wand in every
%   synchronised image frame, then export detections to a CSV usable by
%   OpenLPT wand calibration.
%
% DETECTION PIPELINE (per image):
%   1. Background subtraction (median of first N frames)
%   2. Gaussian blur to suppress noise
%   3. Percentile-based intensity threshold
%   4. Connected-component labelling + area filtering
%   5. Intensity-weighted centroid (sub-pixel) refinement
%   6. Pairing: exactly two blobs expected per image; ranked by intensity
%
% OUTPUTS (saved to resultsDir):
%   wand_detections.mat    — struct array of detections per frame per camera
%   wand_points.csv        — OpenLPT-compatible wand point CSV
%   detection_overview.png — montage of sample detections
% =========================================================================
clear; clc; close all; 

%% =========================================================================
%  CONFIGURATION — edit this section to match your setup
% =========================================================================
% --- Paths ---
mainDir = uigetdir(pwd, 'Select calibration folder (containing cam0/, cam1/, etc.)');

% Check valid dir
if ischar(mainDir) 
    cd(mainDir);
    fprintf('Working directory successfully changed to:\n%s\n', mainDir);
else
    disp('Folder selection cancelled. Script aborted.');
    return;
end

wandDir    = mainDir; % Expects cam0/, cam1/, etc. to live directly inside this folder
resultsDir = mainDir; % Outputs (.mat, .csv) will be saved directly here

% --- Camera / frame settings ---
nCams = 4;       % number of cameras

% --- LED detection parameters ---
minRadius = 20;    % min ball radius (px)
maxRadius = 75;  % max ball radius (px)
detection_method = 'twostage'; %Set to 'phase' or 'twostage'. Phase is slightly more accurate, but much longer
% Edge threshold for imfindcircles, ONE VALUE PER CAMERA (cam0, cam1, ...).
% Lower = more permissive (finds weaker / lower-contrast circles). 0.125 works
% well as a default; thresholdTester.m can help determine per-camera values.
% cam1 has a strong illumination gradient + dark background structure, so it
% needs a lower threshold than the others. A single scalar is also accepted
% and will be applied to every camera.
edgeThresh = [0.15, 0.05, 0.10, 0.15];

% --- Post-detection filters ---
minWandPx = 20;   % min wand pixel length (px)
maxWandPx = 700;   % max wand pixel length (px)
maxJump   = 800;    % max per-frame displacement for temporal filter (px)

% --- CSV export ---
OUTPUT_CSV = fullfile(resultsDir, 'wand_points.csv');
%% =========================================================================
% Wand Tracking  
% =========================================================================
fprintf('=== Wand Tracking ===\n\n');
tic

% -------------------------------------------------------------------------
% 1.  Discover synchronised frames
% -------------------------------------------------------------------------
frameSets = cell(nCams, 1);
camDirs   = cell(nCams, 1);
for c = 1:nCams
    % MODIFIED: Changed 'cam%d' with 'c' to 'cam%d' with 'c - 1' to start from cam0
    camDirs{c} = fullfile(wandDir, sprintf('cam%d', c - 1));
    if ~exist(camDirs{c}, 'dir')
        error('Wand image folder not found: %s', camDirs{c});
    end
    files = dir(fullfile(camDirs{c}, '*.png'));
    if isempty(files), files = dir(fullfile(camDirs{c}, '*.tif')); end
    if isempty(files), files = dir(fullfile(camDirs{c}, '*.jpg')); end
    frameSets{c} = {files.name};
end

% Warn if frame counts differ across cameras
nFrames = numel(frameSets{1});
for c = 2:nCams
    if numel(frameSets{c}) ~= nFrames
        warning('Camera %d has %d frames; Camera 0 has %d. Using min.', ...
            c - 1, numel(frameSets{c}), nFrames);
    end
end
nFrames = min(cellfun(@numel, frameSets));
%% -------------------------------------------------------------------------
% 3.  Detection
% -------------------------------------------------------------------------
%Set detection method 
if strcmp(detection_method, 'phase')
    method = 'PhaseCode';
elseif strcmp(detection_method, 'twostage')
    method = 'TwoStage';
else
    fprintf('Detection method incorrect')
end

% Normalise edgeThresh to one value per camera (scalar -> applied to all)
if isscalar(edgeThresh)
    edgeThresh = repmat(edgeThresh, 1, nCams);
elseif numel(edgeThresh) ~= nCams
    error('edgeThresh must be a scalar or have one value per camera (nCams = %d), but has %d.', ...
        nCams, numel(edgeThresh));
end

detections = cell(nCams, 1);

for c = 1:nCams
    detections{c} = repmat(struct('pts',[],'radii',[],'areas',[],'valid',false,'frame',0, 'metrics', []), nFrames, 1);
end
fprintf('\nDetecting wand points via Circle Hough Transform (Parallelized)...\n');

% Outer loop over cameras (sequential)
for c = 1:nCams
    fprintf('  Processing Camera %d / %d (EdgeThreshold = %.3f)...\n', c - 1, nCams - 1, edgeThresh(c));

    % FOR PARFOR: Pull variables into local variables for clean broadcasting/slicing
    camDir      = camDirs{c};
    frameSet    = frameSets{c};
    edgeThreshC = edgeThresh(c);   % this camera's edge threshold
    
    % Create a local, flat struct array that parfor can slice easily
    camDetections = repmat(struct('pts',[],'radii',[],'areas',[],'valid',false,'frame',0, 'metrics', []), nFrames, 1);
    
    % Inner loop over frames (parallelized)
    parfor f = 1:nFrames
        fpath = fullfile(camDir, frameSet{f});
        im    = double(imread(fpath));
        if size(im,3)==3, im = mean(im,3); end
        
        % Find dark circles directly using the Image Processing Toolbox
        [centers, radii, metric] = imfindcircles(im, [minRadius, maxRadius], 'ObjectPolarity', 'dark', ...
           'Method', method, 'EdgeThreshold', edgeThreshC);
        
        % CHECK: If fewer than 2 circles are found, flag as invalid and skip
        if size(centers, 1) < 2
            camDetections(f).frame = f;
            continue;
        end
        
        % Extract the top 2 strongest detections based on Hough accumulator metric
        centers2 = centers(1:2, :);
        radii2   = radii(1:2);
        metrics2 = metric(1:2);
        
        % Sort these 2 circles by RADIUS (descending) 
        % This ensures row 1 = Large Ball, row 2 = Small Ball
        [~, sortIdx] = sort(radii2, 'descend');
        centers2     = centers2(sortIdx, :);
        radii2       = radii2(sortIdx);
        metrics2     = metrics2(sortIdx);
        
        % Assign directly to local parallel data structure
        camDetections(f).pts     = centers2;         
        camDetections(f).radii   = radii2;           
        camDetections(f).areas   = pi * (radii2.^2); 
        camDetections(f).valid   = true;
        camDetections(f).frame   = f;
        camDetections(f).metrics = metrics2;
        
        % NOTE: In parfor, printing triggers as workers complete chunks out-of-order,
        % but it still provides a useful gauge of background progress.
        if mod(f, 50) == 0
            fprintf('    Worker processed frame %d / %d (Cam %d)\n', f, nFrames, c - 1);
        end
    end
    
    % Assign the fully populated camera array back into the main cell array
    detections{c} = camDetections;
end
toc
fprintf('\n--- Raw detection counts ---\n');
for c = 1:nCams
    nValid = sum([detections{c}.valid]);
    fprintf('  Camera %d: %d / %d frames with valid detection\n', c - 1, nValid, nFrames);
end

%% -------------------------------------------------------------------------
% POST-DETECTION FILTERS
% -------------------------------------------------------------------------
% ---- Wand pixel length filter ----
fprintf('\n--- Wand pixel length filter [%d, %d] px ---\n', minWandPx, maxWandPx);
for c = 1:nCams
    nDropped = 0;
    for f = 1:nFrames
        if ~detections{c}(f).valid, continue; end
        d = norm(detections{c}(f).pts(1,:) - detections{c}(f).pts(2,:));
        if d < minWandPx || d > maxWandPx
            detections{c}(f).valid = false;
            nDropped = nDropped + 1;
        end
    end
    fprintf('  Cam %d: dropped %d frames (bad wand length)\n', c - 1, nDropped);
end

% ---- Temporal consistency filter ----
fprintf('\n--- Temporal consistency filter [max jump = %d px] ---\n', maxJump);
for c = 1:nCams
    nDropped = 0;
    prevValid = 0;
    for f = 1:nFrames
        if ~detections{c}(f).valid, continue; end
        if prevValid == 0
            prevValid = f;
            continue;
        end
        pts_cur  = detections{c}(f).pts;
        pts_prev = detections{c}(prevValid).pts;
        frameGap        = f - prevValid;
        adjustedMaxJump = maxJump * min(frameGap, 5);
        jumpA = norm(pts_cur(1,:) - pts_prev(1,:));
        jumpB = norm(pts_cur(2,:) - pts_prev(2,:));
        if jumpA > adjustedMaxJump || jumpB > adjustedMaxJump
            detections{c}(f).valid = false;
            nDropped = nDropped + 1;
        else
            prevValid = f;
        end
    end
    fprintf('  Cam %d: dropped %d frames (temporal jump)\n', c - 1, nDropped);
end

% ---- Final counts ----
fprintf('\n--- Final detection counts ---\n');
for c = 1:nCams
    nValid = sum([detections{c}.valid]);
    fprintf('  Camera %d: %d / %d valid frames (%.1f%%)\n', ...
        c - 1, nValid, nFrames, 100*nValid/nFrames);
end

% ---- Wand pixel length stats ----
fprintf('\n--- Wand pixel length stats (post-filter) ---\n');
for c = 1:nCams
    lens = [];
    for f = 1:nFrames
        if detections{c}(f).valid
            d = norm(detections{c}(f).pts(1,:) - detections{c}(f).pts(2,:));
            lens(end+1) = d; 
        end
    end
    if ~isempty(lens)
        fprintf('  Cam %d: %.0f +/- %.0f px (range %.0f-%.0f)\n', ...
            c - 1, mean(lens), std(lens), min(lens), max(lens));
    end
end

%% -------------------------------------------------------------------------
% 4.  Save .mat
% -------------------------------------------------------------------------
outFile = fullfile(resultsDir, 'wand_detections.mat');
save(outFile, 'detections', 'nFrames', 'nCams', 'camDirs', 'frameSets');
fprintf('\n[WandTracking] Detections saved to %s\n', outFile);

%% -------------------------------------------------------------------------
% 5.  Export OpenLPT CSV
% -------------------------------------------------------------------------
fprintf('\n--- Exporting OpenLPT CSV ---\n');

% Find frames valid across ALL cameras
common_frames = [];
for c = 1:nCams
    cam_valid_frames = find([detections{c}.valid]);
    if c == 1
        common_frames = cam_valid_frames;
    else
        common_frames = intersect(common_frames, cam_valid_frames);
    end
end

common_frames = sort(common_frames);
fprintf('  Frames valid across all %d cameras: %d\n', nCams, numel(common_frames));
fid = fopen(OUTPUT_CSV, 'w');
fprintf(fid, 'Frame,Camera,Status,PointIdx,X,Y,Radius,Metric\n');

for fi = 1:numel(common_frames)
    f = common_frames(fi);
    for c = 1:nCams
        % Data structure coordinates derived from the imfindcircles pass:
        % row 1 = Large Dark Ball [x, y], row 2 = Small Dark Ball [x, y]
        pts     = detections{c}(f).pts;   
        x_large = pts(1,1);  y_large = pts(1,2);
        x_small = pts(2,1);  y_small = pts(2,2);
        
        cam_id  = c - 1;   % 0-indexed camera configuration mapping for OpenLPT
        
        r_large = detections{c}(f).radii(1);
        r_small = detections{c}(f).radii(2);
        m_large = detections{c}(f).metrics(1);
        m_small = detections{c}(f).metrics(2);
        
        % Small written first (insertion sequence sorting is critical to the OpenLPT tracking engine)
        fprintf(fid, '%d,%d,Filtered_Small,0,%.6f,%.6f,%.4f,%.1f\n', ...
            f, cam_id, x_small, y_small, r_small, m_small);
        fprintf(fid, '%d,%d,Filtered_Large,1,%.6f,%.6f,%.4f,%.1f\n', ...
            f, cam_id, x_large, y_large, r_large, m_large);
    end
end
fclose(fid);
fprintf('  Wrote %d frames x %d cameras to %s\n', numel(common_frames), nCams, OUTPUT_CSV);