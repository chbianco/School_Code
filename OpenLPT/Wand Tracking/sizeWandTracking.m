% =========================================================================
% WandTracking.m
% =========================================================================
% PURPOSE:
%   Automatically detect the two dark circles of the calibration wand in every
%   synchronised image frame, then export detections to a CSV usable by
%   OpenLPT wand calibration.
%
% DETECTION PIPELINE (per image):
%   1. Background subtraction (Inverted: Background - Image for dark targets)
%   2. Gaussian blur to suppress noise
%   3. Percentile-based intensity threshold
%   4. Connected-component labelling + area filtering
%   5. Sub-pixel centroid refinement
%   6. Pairing: exactly two blobs expected per image; ranked by Area (Large vs Small)
% =========================================================================
clear; clc;

% =========================================================================
%  CONFIGURATION
% =========================================================================
% --- Paths ---
wandDir    = 'wand_images';   % parent folder containing cam1/, cam2/, ...
resultsDir = 'results';       % output folder (must exist)

% --- Camera / frame settings ---
nCams      = 4;       % number of cameras
startFrame = 1;       % skip frames before wand appears (set to 1 for no skip)

% --- Dark Ball detection parameters ---
detCfg.gaussSigma    = 1.5;   % Gaussian blur sigma (px)
detCfg.intensityPct  = 95;    % threshold percentile of non-zero difference pixels
detCfg.minArea       = 30;    % min blob area (px²)
detCfg.maxArea       = 1200;  % max blob area (px²)

% --- Post-detection filters ---
minWandPx = 10;   % min wand pixel length (px)
maxWandPx = 800;  % max wand pixel length (px)

% --- Per-camera minimum wand length overrides ---
minWandPxPerCam = containers.Map({1,2,3,4}, {minWandPx, minWandPx, minWandPx, 60});

% --- CSV export parameters ---
METRIC     = 0.5;    % placeholder metric written to CSV
OUTPUT_CSV = fullfile(resultsDir, 'wand_points.csv');

fprintf('=== Dark Wand Tracking ===\n\n');

% -------------------------------------------------------------------------
% 1. Discover synchronised frames
% -------------------------------------------------------------------------
frameSets = cell(nCams, 1);
camDirs   = cell(nCams, 1);
for c = 1:nCams
    camDirs{c} = fullfile(wandDir, sprintf('cam%d', c));
    if ~exist(camDirs{c}, 'dir')
        error('Wand image folder not found: %s', camDirs{c});
    end
    files = dir(fullfile(camDirs{c}, '*.png'));
    if isempty(files), files = dir(fullfile(camDirs{c}, '*.tif')); end
    if isempty(files), files = dir(fullfile(camDirs{c}, '*.jpg')); end
    frameSets{c} = {files.name};
end

nFrames = min(cellfun(@numel, frameSets));
if startFrame > 1
    for c = 1:nCams
        frameSets{c} = frameSets{c}(startFrame:end);
    end
    nFrames = min(cellfun(@numel, frameSets));
    fprintf('Skipping first %d frames. Processing %d frames.\n', startFrame-1, nFrames);
end
fprintf('Processing %d synchronised frames across %d cameras\n\n', nFrames, nCams);

% -------------------------------------------------------------------------
% 2. Compute background model (median filter over initial frames)
% -------------------------------------------------------------------------
fprintf('Computing background models...\n');
bgFrames = min(20, nFrames);
bg = cell(nCams, 1);
for c = 1:nCams
    stack = [];
    for i = 1:bgFrames
        fpath = fullfile(camDirs{c}, frameSets{c}{i});
        im    = imread(fpath);
        if size(im,3)==3, im = rgb2gray(im); end
        stack(:,:,i) = double(im); 
    end
    bg{c} = median(stack, 3);
    fprintf('  Camera %d background computed from %d frames\n', c, bgFrames);
end

% -------------------------------------------------------------------------
% 3. Detect Dark Balls in every frame
% -------------------------------------------------------------------------
detections = cell(nCams, 1);
for c = 1:nCams
    detections{c} = repmat(struct('pts',[],'radii',[],'areas',[],'valid',false,'frame',0), nFrames, 1);
end

fprintf('\nTracking targets...\n');
for f = 1:nFrames
    for c = 1:nCams
        if f <= bgFrames
            detections{c}(f).frame = f;
            continue;
        end
              
        fpath = fullfile(camDirs{c}, frameSets{c}{f});
        im    = double(imread(fpath));
        if size(im,3)==3, im = mean(im,3); end
        
        % CRITICAL CHANGE: Invert subtraction because target features are DARKER than background
        diff_im = max(bg{c} - im, 0);
        
        % Gaussian blur
        sigma   = detCfg.gaussSigma;
        h       = fspecial('gaussian', ceil(6*sigma+1), sigma);
        blurred = imfilter(diff_im, h, 'replicate');
        
        % Threshold dynamic range
        vals = blurred(blurred > 0);
        if isempty(vals)
            detections{c}(f).frame = f;
            continue;
        end
        thresh = prctile(vals, detCfg.intensityPct);
        mask   = blurred >= thresh;
        
        % Morphological cleanup
        mask = bwareaopen(mask, detCfg.minArea);
        mask = imclose(mask, strel('disk', 4));
        
        % Feature Extraction
        cc    = bwconncomp(mask);
        props = regionprops(cc, blurred, 'Area','WeightedCentroid','PixelIdxList','MinorAxisLength');
        
        % Filter components by size boundaries
        areas = [props.Area];
        keep  = areas >= detCfg.minArea & areas <= detCfg.maxArea;
        props = props(keep);
        
        if numel(props) < 2
            detections{c}(f).frame = f;
            continue;
        end
        
        % CRITICAL CHANGE: Sort descending by Area instead of Intensity to separate Large vs Small
        [~, order] = sort([props.Area], 'descend');
        props      = props(order(1:2));
        
        % Sub-pixel centroid configuration
        pts   = zeros(2,2);
        radii = zeros(2,1);
        for p = 1:2
            [rows, cols] = ind2sub(size(blurred), props(p).PixelIdxList);
            weights      = blurred(props(p).PixelIdxList);
            weights      = weights / sum(weights);
            pts(p,1)     = sum(cols .* weights);   % u (x)
            pts(p,2)     = sum(rows .* weights);   % v (y)
            radii(p)     = props(p).MinorAxisLength / 2; 
        end
        
        detections{c}(f).pts    = pts;
        detections{c}(f).radii  = radii;  % radii(1) = Large, radii(2) = Small
        detections{c}(f).areas  = [props(1).Area; props(2).Area];
        detections{c}(f).valid  = true;
        detections{c}(f).frame  = f;
    end
    if mod(f,50)==0
        fprintf('  Frame %d / %d\n', f, nFrames);
    end
end

% -------------------------------------------------------------------------
% 4. Post-Detection Wand Geometry Validation
% -------------------------------------------------------------------------
fprintf('\n--- Wand pixel length validation ---\n');
for c = 1:nCams
    nDropped = 0;
    for f = 1:nFrames
        if ~detections{c}(f).valid, continue; end
        d = norm(detections{c}(f).pts(1,:) - detections{c}(f).pts(2,:));
        if d < minWandPxPerCam(c) || d > maxWandPx
            detections{c}(f).valid = false;
            nDropped = nDropped + 1;
        end
    end
    fprintf('  Cam %d: dropped %d frames failing geometric length restrictions\n', c, nDropped);
end

% Summary Metrics
fprintf('\n--- Final valid tracking counts ---\n');
for c = 1:nCams
    nValid = sum([detections{c}.valid]);
    fprintf('  Camera %d: %d / %d valid frames (%.1f%%)\n', c, nValid, nFrames, 100*nValid/nFrames);
end

% -------------------------------------------------------------------------
% 5. Save Workspace and Export OpenLPT CSV
% -------------------------------------------------------------------------
outFile = fullfile(resultsDir, 'wand_detections.mat');
save(outFile, 'detections', 'nFrames', 'nCams', 'camDirs', 'frameSets');

fprintf('\n--- Exporting OpenLPT CSV ---\n');
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
fprintf('  Frames fully intersecting across all %d cameras: %d\n', nCams, numel(common_frames));

fid = fopen(OUTPUT_CSV, 'w');
fprintf(fid, 'Frame,Camera,Status,PointIdx,X,Y,Radius,Metric\n');
for fi = 1:numel(common_frames)
    f = common_frames(fi);
    for c = 1:nCams
        pts = detections{c}(f).pts;   % row1 = Large, row2 = Small
        x_large = pts(1,1);  y_large = pts(1,2);
        x_small = pts(2,1);  y_small = pts(2,2);
        cam_id  = c - 1;              % 0-indexed base conversion
        r_large = detections{c}(f).radii(1);
        r_small = detections{c}(f).radii(2);
        
        % Write Small target first, followed by Large target
        fprintf(fid, '%d,%d,Filtered_Small,0,%.6f,%.6f,%.4f,%.1f\n', f, cam_id, x_small, y_small, r_small, METRIC);
        fprintf(fid, '%d,%d,Filtered_Large,1,%.6f,%.6f,%.4f,%.1f\n', f, cam_id, x_large, y_large, r_large, METRIC);
    end
end
fclose(fid);
fprintf('  Wrote parsing records to %s\n', OUTPUT_CSV);

% Generate verification diagnostics
generateDetectionOverview(detections, camDirs, frameSets, nCams, nFrames, resultsDir, startFrame);
fprintf('\n=== Process Complete ===\n');

% =========================================================================
%  HELPER VISUALIZATION FUNCTION
% =========================================================================
function generateDetectionOverview(detections, camDirs, frameSets, nCams, nFrames, resultsDir, startFrame)
    nSamples = 5;
    fig = figure('Visible','off', 'Position', [0 0 1800 1000]);
    for c = 1:nCams
        validIdx = find([detections{c}.valid]);
        if isempty(validIdx), continue; end
        pickIdx = round(linspace(1, numel(validIdx), min(nSamples, numel(validIdx))));
        pick    = validIdx(pickIdx);
        for fi = 1:numel(pick)
            f  = pick(fi);
            im = imread(fullfile(camDirs{c}, frameSets{c}{f}));
            if size(im,3)==3, im = rgb2gray(im); end
            
            im_disp = im;
            hi = prctile(double(im_disp(:)), 99.5);
            if hi > 0, im_disp = uint8(double(im_disp) * (255/hi)); end
            
            ax = subplot(nCams, nSamples, (c-1)*nSamples + fi);
            imagesc(im_disp, 'Parent', ax);
            colormap(ax, gray); caxis(ax, [0 255]);
            axis(ax, 'image', 'off');
            hold(ax, 'on');
            
            det = detections{c}(f);
            % Large Ball — Red circle
            plot(ax, det.pts(1,1), det.pts(1,2), 'ro', 'MarkerSize', 14, 'LineWidth', 2);
            % Small Ball — Green circle
            plot(ax, det.pts(2,1), det.pts(2,2), 'go', 'MarkerSize', 14, 'LineWidth', 2);
            % Linking Vector
            plot(ax, det.pts(:,1), det.pts(:,2), 'y-', 'LineWidth', 1.5);
            
            actualFrame = f + startFrame - 1;
            if fi == 1
                title(ax, sprintf('Cam%d  f=%d', c-1, actualFrame), 'FontSize', 7);
            else
                title(ax, sprintf('f=%d', actualFrame), 'FontSize', 7);
            end
            hold(ax, 'off');
        end
    end
    sgtitle('Detections per camera  |  red = Large Ball   green = Small Ball', 'FontSize', 10);
    saveas(fig, fullfile(resultsDir, 'detection_overview.png'));
    close(fig);
    fprintf('[WandTracking] Diagnostic overview montage compiled successfully.\n');
end