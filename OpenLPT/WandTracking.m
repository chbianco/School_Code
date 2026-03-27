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

clear; clc;

% =========================================================================
%  CONFIGURATION — edit this section to match your setup
% =========================================================================

% --- Paths ---
wandDir    = 'wand_images';   % parent folder containing cam1/, cam2/, ...
resultsDir = 'results';       % output folder (must exist)

% --- Camera / frame settings ---
nCams      = 4;       % number of cameras
startFrame = 1;     % skip frames before wand appears (set to 1 for no skip)

% --- LED detection parameters ---
detCfg.gaussSigma    = 1.5;   % Gaussian blur sigma (px)
detCfg.intensityPct  = 93;    % threshold at this percentile of non-zero pixels
detCfg.minArea       = 20;    % min blob area (px²)
detCfg.maxArea       = 800;  % max blob area (px²)

% --- Post-detection filters ---
minWandPx = 20;   % min wand pixel length (px)
maxWandPx = 700;   % max wand pixel length (px)
maxJump   = 150;    % max per-frame displacement for temporal filter (px)
% --- Per-camera minimum wand length overrides ---
minWandPxPerCam = containers.Map({1,2,3,4}, {minWandPx, minWandPx, minWandPx, 60});

% --- CSV export parameters ---
%RADIUS     = 10.0;   % placeholder radius written to CSV
METRIC     = 0.5;    % placeholder metric written to CSV
OUTPUT_CSV = fullfile(resultsDir, 'wand_points.csv');

% --- Brightness-based identification ---
% Large = brighter LED, Small = dimmer LED (assigned per-frame by intensity).
% Requires one LED to be physically dimmer than the other.
% minBrightnessRatio: if the two LEDs are closer in brightness than this
% ratio (bright/dim < threshold), the frame is flagged as ambiguous and
% marked invalid.  Set to 1.0 to disable.
minBrightnessRatio = 1.05;
minBrightnessRatioPerCam = containers.Map({1,2,3,4}, {minBrightnessRatio, 1.15, minBrightnessRatio, 1.5});

% --- Diagnostics ---

% =========================================================================
%  END CONFIGURATION
% =========================================================================

fprintf('=== Wand Tracking ===\n\n');

% -------------------------------------------------------------------------
% 1.  Discover synchronised frames
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

% Warn if frame counts differ across cameras
nFrames = numel(frameSets{1});
for c = 2:nCams
    if numel(frameSets{c}) ~= nFrames
        warning('Camera %d has %d frames; Camera 1 has %d. Using min.', ...
            c, numel(frameSets{c}), nFrames);
    end
end
nFrames = min(cellfun(@numel, frameSets));

% Skip background frames before wand appears
if startFrame > 1
    for c = 1:nCams
        frameSets{c} = frameSets{c}(startFrame:end);
    end
    nFrames = min(cellfun(@numel, frameSets));
    fprintf('Skipping first %d background frames. Processing %d frames.\n', startFrame-1, nFrames);
end

fprintf('Processing %d synchronised frames across %d cameras\n\n', nFrames, nCams);

% -------------------------------------------------------------------------
% 2.  Compute background model (median of first 20 frames or all if fewer)
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
% 3.  Detect LEDs in every frame
% -------------------------------------------------------------------------
% detections{c}(f) = struct with fields:
%   .pts      [2x2] sub-pixel centroid positions (row = point, col = [u v])
%             row 1 = brighter / larger ball  -> Filtered_Large in CSV
%             row 2 = dimmer  / smaller ball  -> Filtered_Small in CSV
%   .valid    logical — true if exactly 2 blobs found
%   .frame    frame index

detections = cell(nCams, 1);
for c = 1:nCams
    detections{c} = repmat(struct('pts',[],'radii',[],'intens',[],'valid',false,'frame',0), nFrames, 1);
end

fprintf('\nDetecting LEDs...\n');
for f = 1:nFrames
    for c = 1:nCams

         % Skip background frames — no wand present
               if f <= bgFrames
                   detections{c}(f).frame = f;
                   continue;
               end
              
        fpath = fullfile(camDirs{c}, frameSets{c}{f});
        im    = double(imread(fpath));
        if size(im,3)==3, im = mean(im,3); end

        % Background subtract + clip
        diff_im = max(im - bg{c}, 0);

        % Gaussian blur
        sigma   = detCfg.gaussSigma;
        h       = fspecial('gaussian', ceil(6*sigma+1), sigma);
        blurred = imfilter(diff_im, h, 'replicate');

        % Threshold at Nth percentile of non-zero pixels
        vals = blurred(blurred > 0);
        if isempty(vals)
            detections{c}(f).frame = f;
            continue;
        end
        thresh = prctile(vals, detCfg.intensityPct);
        mask   = blurred >= thresh;

        % Morphological clean-up
        mask = bwareaopen(mask, detCfg.minArea);
        mask = imclose(mask, strel('disk', 4));

        % Connected components
        cc    = bwconncomp(mask);
        props = regionprops(cc, blurred, ...
            'Area','WeightedCentroid','MeanIntensity','PixelIdxList','MinorAxisLength');

        % Filter by area
        areas = [props.Area];
        keep  = areas >= detCfg.minArea & areas <= detCfg.maxArea;
        props = props(keep);

        if numel(props) < 2
            detections{c}(f).frame = f;
            continue;
        end

        % Sort by mean intensity descending, take top 2
        [~, order] = sort([props.MeanIntensity], 'descend');
        props      = props(order(1:2));

        % Sub-pixel centroid (intensity-weighted)
        pts   = zeros(2,2);
        radii = zeros(2,1);
        for p = 1:2
            [rows, cols] = ind2sub(size(blurred), props(p).PixelIdxList);
            weights      = blurred(props(p).PixelIdxList);
            weights      = weights / sum(weights);
            pts(p,1)     = sum(cols .* weights);        % u (x)
            pts(p,2)     = sum(rows .* weights);        % v (y)
            radii(p) = props(p).MinorAxisLength / 2;  % radius from ellipse minor axis (robust to viewing angle)
        end

        detections{c}(f).pts    = pts;
        detections{c}(f).radii  = radii;            % radii(1)=Large, radii(2)=Small
        detections{c}(f).intens = [props(1).MeanIntensity; props(2).MeanIntensity];
        detections{c}(f).valid  = true;
        detections{c}(f).frame  = f;
    end

    if mod(f,50)==0
        fprintf('  Frame %d / %d\n', f, nFrames);
    end
end

fprintf('\n--- Raw detection counts ---\n');
for c = 1:nCams
    nValid = sum([detections{c}.valid]);
    fprintf('  Camera %d: %d / %d frames with valid detection\n', c, nValid, nFrames);
end

% -------------------------------------------------------------------------
% 3b. POST-DETECTION FILTERS
% -------------------------------------------------------------------------

% ---- Wand pixel length filter ----
fprintf('\n--- Wand pixel length filter [%d, %d] px ---\n', minWandPx, maxWandPx);
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
    fprintf('  Cam %d: dropped %d frames (bad wand length)\n', c, nDropped);
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
    fprintf('  Cam %d: dropped %d frames (temporal jump)\n', c, nDropped);
end

% -------------------------------------------------------------------------
% 3c. VALIDATE BRIGHTNESS SEPARATION & REJECT AMBIGUOUS FRAMES
% -------------------------------------------------------------------------
% Large/Small assignment is already done by the intensity sort in the
% detection loop (row 1 = brighter = Large, row 2 = dimmer = Small).
% Here we reject frames where the two LEDs are too similar in brightness
% to be reliably distinguished, and report per-camera statistics.
fprintf('\n--- Brightness-based Large/Small validation ---\n');

for c = 1:nCams
    ratios   = [];
    nDropped = 0;
    for f = 1:nFrames
        if ~detections{c}(f).valid, continue; end
        intens = detections{c}(f).intens;   % [bright, dim]
        if intens(2) < 1e-6
            detections{c}(f).valid = false;
            nDropped = nDropped + 1;
            continue;
        end
        ratio = intens(1) / intens(2);
        ratios(end+1) = ratio;
        if ratio < minBrightnessRatioPerCam(c)
            detections{c}(f).valid = false;
            nDropped = nDropped + 1;
        end
    end
    if isempty(ratios)
        fprintf('  Camera %d: no valid frames\n', c);
    else
        fprintf('  Camera %d: brightness ratio bright/dim = %.2f +/- %.2f (min %.2f)  |  dropped %d ambiguous frames\n', ...
            c, mean(ratios), std(ratios), min(ratios), nDropped);
    end
end

% ---- Final counts ----
fprintf('\n--- Final detection counts ---\n');
for c = 1:nCams
    nValid = sum([detections{c}.valid]);
    fprintf('  Camera %d: %d / %d valid frames (%.1f%%)\n', ...
        c, nValid, nFrames, 100*nValid/nFrames);
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
            c, mean(lens), std(lens), min(lens), max(lens));
    end
end

% -------------------------------------------------------------------------
% 4.  Save .mat
% -------------------------------------------------------------------------
outFile = fullfile(resultsDir, 'wand_detections.mat');
save(outFile, 'detections', 'nFrames', 'nCams', 'camDirs', 'frameSets');
fprintf('\n[WandTracking] Detections saved to %s\n', outFile);

% -------------------------------------------------------------------------
% 5.  Export OpenLPT CSV
% -------------------------------------------------------------------------
fprintf('\n--- Exporting OpenLPT CSV ---\n');

% Find frames valid in ALL cameras
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
        pts = detections{c}(f).pts;   % row1 = Large, row2 = Small

        x_large = pts(1,1);  y_large = pts(1,2);
        x_small = pts(2,1);  y_small = pts(2,2);

        cam_id = c - 1;   % 0-indexed for OpenLPT

        r_large = detections{c}(f).radii(1);
        r_small = detections{c}(f).radii(2);

        % Small written first (insertion order matters in OpenLPT loader)
        fprintf(fid, '%d,%d,Filtered_Small,0,%.6f,%.6f,%.4f,%.1f\n', ...
            f, cam_id, x_small, y_small, r_small, METRIC);
        fprintf(fid, '%d,%d,Filtered_Large,1,%.6f,%.6f,%.4f,%.1f\n', ...
            f, cam_id, x_large, y_large, r_large, METRIC);
    end
end

fclose(fid);
fprintf('  Wrote %d frames x %d cameras to %s\n', numel(common_frames), nCams, OUTPUT_CSV);

% -------------------------------------------------------------------------
% 6.  Detection overview montage
% -------------------------------------------------------------------------
generateDetectionOverview(detections, camDirs, frameSets, nCams, nFrames, resultsDir, startFrame);

fprintf('\n=== WandTracking complete ===\n');


% =========================================================================
%  HELPER: Detection overview figure
% =========================================================================
function generateDetectionOverview(detections, camDirs, frameSets, nCams, nFrames, resultsDir, startFrame)
    % Show 5 evenly-spaced valid frames for each camera independently.
    % Red circle = Large (bright) LED, green circle = Small (dim) LED.
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

            % Stretch contrast so dark images are visible
            im_disp = im;
            hi = prctile(double(im_disp(:)), 99.5);
            if hi > 0, im_disp = uint8(double(im_disp) * (255/hi)); end

            ax = subplot(nCams, nSamples, (c-1)*nSamples + fi);
            imagesc(im_disp, 'Parent', ax);
            colormap(ax, gray); caxis(ax, [0 255]);
            axis(ax, 'image', 'off');
            hold(ax, 'on');

            det = detections{c}(f);
            % Large LED — red circle
            plot(ax, det.pts(1,1), det.pts(1,2), 'ro', ...
                'MarkerSize', 14, 'LineWidth', 2);
            % Small LED — green circle
            plot(ax, det.pts(2,1), det.pts(2,2), 'go', ...
                'MarkerSize', 14, 'LineWidth', 2);
            % Line connecting them
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

    sgtitle('Detections per camera  |  red = Large (bright)   green = Small (dim)', 'FontSize', 10);
    saveas(fig, fullfile(resultsDir, 'detection_overview.png'));
    close(fig);
    fprintf('[WandTracking] Detection overview saved.\n');
end
