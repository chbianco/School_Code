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
detCfg.intensityPct  = 97;    % threshold at this percentile of non-zero pixels
detCfg.minArea       = 20;    % min blob area (px²)
detCfg.maxArea       = 800;  % max blob area (px²)

% --- Post-detection filters ---
minWandPx = 80;   % min wand pixel length (px)
maxWandPx = 900;   % max wand pixel length (px)
maxJump   = 150;    % max per-frame displacement for temporal filter (px)

% --- CSV export parameters ---
RADIUS     = 10.0;   % placeholder radius written to CSV
METRIC     = 0.5;    % placeholder metric written to CSV
OUTPUT_CSV = fullfile(resultsDir, 'wand_points.csv');

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
    detections{c} = repmat(struct('pts',[],'valid',false,'frame',0), nFrames, 1);
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
            'Area','WeightedCentroid','MeanIntensity','PixelIdxList');

        % Filter by area
        areas = [props.Area];
        keep  = areas >= detCfg.minArea & areas <= detCfg.maxArea;
        props = props(keep);

        if numel(props) < 2
            detections{c}(f).frame = f;
            continue;
        end

        % Sort by mean intensity descending, take top 2
        % props(1) = brighter = Large ball, props(2) = dimmer = Small ball
        [~, order] = sort([props.MeanIntensity], 'descend');
        props      = props(order(1:2));

        % Sub-pixel centroid (intensity-weighted)
        pts = zeros(2,2);
        for p = 1:2
            [rows, cols] = ind2sub(size(blurred), props(p).PixelIdxList);
            weights      = blurred(props(p).PixelIdxList);
            weights      = weights / sum(weights);
            pts(p,1)     = sum(cols .* weights);   % u (x)
            pts(p,2)     = sum(rows .* weights);   % v (y)
        end

        detections{c}(f).pts   = pts;
        detections{c}(f).valid = true;
        detections{c}(f).frame = f;
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
        if d < minWandPx || d > maxWandPx
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
  % 3c. ENFORCE LARGE/SMALL TEMPORAL CONSISTENCY VIA TRACKING
  % -------------------------------------------------------------------------
  fprintf('\n--- Enforcing Large/Small temporal consistency ---\n');
  for c = 1:nCams
      prevValid = 0;
      nSwapped  = 0;
      for f = 1:nFrames
          if ~detections{c}(f).valid, continue; end

          pts = detections{c}(f).pts;  % row 1 = current Large, row 2 = current Small

          if prevValid == 0
              % First valid frame: initialize so that Large = blob with higher Y
              % (lower in image = same physical LED for all upright cameras)
              if pts(1,2) < pts(2,2)   % row 1 is actually higher (smaller Y)
                  detections{c}(f).pts = pts([2,1], :);  % swap
                  nSwapped = nSwapped + 1;
              end
          else
              % Subsequent frames: track which blob is closest to previous Large
              pts_prev = detections{c}(prevValid).pts;
              d_keep = norm(pts(1,:) - pts_prev(1,:));  % curr Large → prev Large
              d_swap = norm(pts(2,:) - pts_prev(1,:));  % curr Small → prev Large

              if d_swap < d_keep
                  detections{c}(f).pts = pts([2,1], :);  % swap to maintain identity
                  nSwapped = nSwapped + 1;
              end
          end

          prevValid = f;
      end
      fprintf('  Camera %d: %d frames re-assigned for consistency\n', c, nSwapped);
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

        % Small written first (insertion order matters in OpenLPT loader)
        fprintf(fid, '%d,%d,Filtered_Small,0,%.6f,%.6f,%.1f,%.1f\n', ...
            f, cam_id, x_small, y_small, RADIUS, METRIC);
        fprintf(fid, '%d,%d,Filtered_Large,1,%.6f,%.6f,%.1f,%.1f\n', ...
            f, cam_id, x_large, y_large, RADIUS, METRIC);
    end
end

fclose(fid);
fprintf('  Wrote %d frames x %d cameras to %s\n', numel(common_frames), nCams, OUTPUT_CSV);

% -------------------------------------------------------------------------
% 6.  Detection overview montage
% -------------------------------------------------------------------------
generateDetectionOverview(detections, camDirs, frameSets, nCams, nFrames, resultsDir);

fprintf('\n=== WandTracking complete ===\n');


% =========================================================================
%  HELPER: Detection overview figure
% =========================================================================
function generateDetectionOverview(detections, camDirs, frameSets, nCams, nFrames, resultsDir)
    % Pick up to 8 evenly spaced frames valid across all cameras
    validFrames = [];
    for f = 1:nFrames
        allValid = all(arrayfun(@(c) detections{c}(f).valid, 1:nCams));
        if allValid, validFrames(end+1) = f; end %#ok<AGROW>
    end
    if isempty(validFrames), return; end

    pick = round(linspace(1, numel(validFrames), min(8, numel(validFrames))));
    pick = validFrames(pick);

    fig = figure('Visible','off','Position',[0 0 1600 900]);
    nPick = numel(pick);
    for fi = 1:nPick
        f  = pick(fi);
        im = imread(fullfile(camDirs{1}, frameSets{1}{f}));
        if size(im,3)==3, im = rgb2gray(im); end
        ax = subplot(2, ceil(nPick/2), fi);
        imagesc(im,'Parent',ax); colormap(ax,gray); axis(ax,'image','off');
        hold(ax,'on');
        det = detections{1}(f);
        if det.valid
            plot(ax, det.pts(1,1), det.pts(1,2), 'r+', 'MarkerSize',12,'LineWidth',2);  % Large
            plot(ax, det.pts(2,1), det.pts(2,2), 'g+', 'MarkerSize',12,'LineWidth',2);  % Small
        end
        title(ax, sprintf('Frame %d', f));
        hold(ax,'off');
    end
    sgtitle('Wand Detections — Camera 1 Sample (red=Large, green=Small)');
    saveas(fig, fullfile(resultsDir, 'detection_overview.png'));
    close(fig);
    fprintf('[WandTracking] Detection overview saved.\n');
end
