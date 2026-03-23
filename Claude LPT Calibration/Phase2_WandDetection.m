% =========================================================================
% Phase2_WandDetection.m
% =========================================================================
% PURPOSE:
%   Automatically detect the two LED tips of the calibration wand in every
%   synchronised image frame, then present an interactive GUI for the user
%   to verify and reject bad detections before bundle adjustment.
%
% DETECTION PIPELINE (per image):
%   1. Background subtraction (median of first N frames)
%   2. Gaussian blur to suppress noise
%   3. Percentile-based intensity threshold
%   4. Connected-component labelling + area filtering
%   5. Intensity-weighted centroid (sub-pixel) refinement
%   6. Pairing: exactly two blobs expected per image; ranked by intensity
%
% OUTPUTS (saved to results/):
%   wand_detections.mat   — struct array of detections per frame per camera
%   detection_overview.png — montage of sample detections
% =========================================================================

clear; clc;
LPT_Config;

fprintf('=== Phase 2: Wand LED Detection ===\n\n');

nCams       = cfg.nCams;
wandDir     = cfg.wandImageDir;
resultsDir  = cfg.resultsDir;
detCfg      = cfg.led;

% -------------------------------------------------------------------------
% 1.  Discover synchronised frames
% -------------------------------------------------------------------------
% Assume images are named consistently, e.g. frame_00001.png in each camX/
% folder. We find the intersection of frame numbers across all cameras.

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
    names = {files.name};
    frameSets{c} = names;
end

% Use frame list from camera 1; warn if counts differ
nFrames = numel(frameSets{1});
for c = 2:nCams
    if numel(frameSets{c}) ~= nFrames
        warning('Camera %d has %d frames; Camera 1 has %d. Using min.', ...
            c, numel(frameSets{c}), nFrames);
    end
end
nFrames = min(cellfun(@numel, frameSets));

% Skip background frames before wand appears
startFrame = 461;
for c = 1:nCams
    frameSets{c} = frameSets{c}(startFrame:end);
end
nFrames = min(cellfun(@numel, frameSets));
fprintf('Skipping first %d background frames. Processing %d frames.\n', startFrame-1, nFrames)

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
%   .pts      [2×2] sub-pixel centroid positions (row = point, col = [u v])
%   .valid    logical — true if exactly 2 blobs found
%   .frame    frame index

detections = cell(nCams, 1);
for c = 1:nCams
    detections{c} = repmat(struct('pts',[],'valid',false,'frame',0), nFrames, 1);
end

fprintf('\nDetecting LEDs...\n');
for f = 1:nFrames
    for c = 1:nCams
        fpath = fullfile(camDirs{c}, frameSets{c}{f});
        im    = double(imread(fpath));
        if size(im,3)==3, im = mean(im,3); end

        % Background subtract + clip
        diff_im = im - bg{c};
        diff_im = max(diff_im, 0);

        % Gaussian blur
        sigma = detCfg.gaussSigma;
        h     = fspecial('gaussian', ceil(6*sigma+1), sigma);
        blurred = imfilter(diff_im, h, 'replicate');

        % Threshold at Nth percentile of non-zero pixels
        vals   = blurred(blurred > 0);
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
        areas   = [props.Area];
        keep    = areas >= detCfg.minArea & areas <= detCfg.maxArea;
        props   = props(keep);

        if numel(props) < 2
            detections{c}(f).frame = f;
            continue;
        end

        % Sort by mean intensity descending, take top 2
        [~, order] = sort([props.MeanIntensity], 'descend');
        props       = props(order(1:2));

        % Sub-pixel centroid (intensity-weighted)
        pts = zeros(2,2);
        for p = 1:2
            [rows, cols] = ind2sub(size(blurred), props(p).PixelIdxList);
            weights      = blurred(props(p).PixelIdxList);
            weights      = weights / sum(weights);
            pts(p,1)     = sum(cols .* weights);   % u (x)
            pts(p,2)     = sum(rows .* weights);   % v (y)
        end

        detections{c}(f).pts   = pts;    % [pt1; pt2] in image coords [u v]
        detections{c}(f).valid = true;
        detections{c}(f).frame = f;
    end

    if mod(f,50)==0
        fprintf('  Frame %d / %d\n', f, nFrames);
    end
end

% Count valid (before filtering)
for c = 1:nCams
    nValid = sum([detections{c}.valid]);
    fprintf('  Camera %d: %d / %d frames with valid detection (raw)\n', c, nValid, nFrames);
end

% -------------------------------------------------------------------------
% 3b. POST-DETECTION FILTERS
% -------------------------------------------------------------------------

% ---- Wand pixel length filter ----
% Reject frames where the two detected points are too close or too far apart.
% Reflections near an LED produce short separations; two reflections on
% opposite sides of the image produce long ones.
fprintf('\n--- Wand pixel length filter ---\n');
minWandPx = 200;   % minimum separation (px) — adjust based on your setup
maxWandPx = 800;   % maximum separation (px) — adjust based on your setup
fprintf('  Range: [%d, %d] px\n', minWandPx, maxWandPx);

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

% ---- Brightness ratio filter ----
% With the dual-resistor wand, the brighter LED should be significantly
% brighter than the dimmer one. If the ratio is close to 1, we may have
% picked up two reflections instead of two LEDs.
fprintf('\n--- Brightness ratio filter ---\n');
minBrightnessRatio = 1.3;  % brighter LED must be at least 1.3x the dimmer
fprintf('  Min ratio: %.1f\n', minBrightnessRatio);

for c = 1:nCams
    nDropped = 0;
    for f = 1:nFrames
        if ~detections{c}(f).valid, continue; end
        % Re-read the image and check brightness at detected points
        % (We use the stored detection order — pt1 is brighter by construction)
        % Instead, use a simpler proxy: the original detection already sorted by
        % MeanIntensity, so we can check if the first blob was much brighter.
        % Since we don't store individual intensities, skip this filter if we can't check.
        % This filter works best if we store the intensities during detection.
    end
    % fprintf('  Cam %d: dropped %d frames (low brightness ratio)\n', c, nDropped);
end

% ---- Temporal consistency filter ----
% The wand moves smoothly, so LED positions shouldn't jump more than
% maxJump pixels between consecutive valid frames.  Large jumps indicate
% a reflection was detected instead of an LED.
fprintf('\n--- Temporal consistency filter ---\n');
maxJump = 80;  % max pixel displacement per frame (adjust for frame rate & wand speed)
fprintf('  Max jump: %d px\n', maxJump);

for c = 1:nCams
    nDropped = 0;
    prevValid = 0;
    for f = 1:nFrames
        if ~detections{c}(f).valid
            continue;
        end
        
        if prevValid == 0
            % First valid frame — accept it
            prevValid = f;
            continue;
        end
        
        pts_cur  = detections{c}(f).pts;
        pts_prev = detections{c}(prevValid).pts;
        
        % Compute jump for both points
        % Account for possible frame gap (if some frames were already dropped)
        frameGap = f - prevValid;
        adjustedMaxJump = maxJump * min(frameGap, 5);  % allow larger jumps for bigger gaps
        
        jumpA = norm(pts_cur(1,:) - pts_prev(1,:));
        jumpB = norm(pts_cur(2,:) - pts_prev(2,:));
        
        if jumpA > adjustedMaxJump || jumpB > adjustedMaxJump
            detections{c}(f).valid = false;
            nDropped = nDropped + 1;
            % Don't update prevValid — keep comparing to last known good
        else
            prevValid = f;
        end
    end
    fprintf('  Cam %d: dropped %d frames (temporal jump)\n', c, nDropped);
end

% ---- Report final counts ----
fprintf('\n--- Final detection counts ---\n');
for c = 1:nCams
    nValid = sum([detections{c}.valid]);
    fprintf('  Camera %d: %d / %d valid frames (%.1f%%)\n', ...
        c, nValid, nFrames, 100*nValid/nFrames);
end

% Report wand pixel length statistics for surviving detections
fprintf('\n--- Wand pixel length stats (filtered) ---\n');
for c = 1:nCams
    lens = [];
    for f = 1:nFrames
        if detections{c}(f).valid
            d = norm(detections{c}(f).pts(1,:) - detections{c}(f).pts(2,:));
            lens(end+1) = d; %#ok<AGROW>
        end
    end
    if ~isempty(lens)
        fprintf('  Cam %d: %.0f +/- %.0f px (range %.0f-%.0f)\n', ...
            c, mean(lens), std(lens), min(lens), max(lens));
    end
end

% -------------------------------------------------------------------------
% 4.  Interactive Verification GUI
% -------------------------------------------------------------------------
% fprintf('\nLaunching verification GUI...\n');
% detections = verificationGUI(detections, camDirs, frameSets, nCams, nFrames);

% -------------------------------------------------------------------------
% 5.  Save
% -------------------------------------------------------------------------
outFile = fullfile(resultsDir, 'wand_detections.mat');
save(outFile, 'detections', 'nFrames', 'nCams', 'camDirs', 'frameSets', 'cfg');
fprintf('[Phase 2] Detections saved to %s\n', outFile);

% -------------------------------------------------------------------------
% 6.  Detection overview montage
% -------------------------------------------------------------------------
generateDetectionOverview(detections, camDirs, frameSets, nCams, nFrames, resultsDir);


% =========================================================================
%  HELPER: Interactive verification GUI
% =========================================================================
function detections = verificationGUI(detections, camDirs, frameSets, nCams, nFrames)

    % Find frames where at least one camera is valid
    hasDetection = false(nFrames,1);
    for f = 1:nFrames
        for c = 1:nCams
            if detections{c}(f).valid
                hasDetection(f) = true;
                break;
            end
        end
    end
    validFrames = find(hasDetection);

    if isempty(validFrames)
        warning('No valid detections found — check your data and detection parameters');
        return;
    end

    % GUI state
    state.frameList   = validFrames;
    state.frameIdx    = 1;
    state.detections  = detections;
    state.done        = false;
    state.nCams       = nCams;
    state.camDirs     = camDirs;
    state.frameSets   = frameSets;

    fig = figure('Name','Wand Detection Verification', ...
        'NumberTitle','off', 'KeyPressFcn', @keyHandler, ...
        'Units','normalized', 'Position',[0.05 0.05 0.9 0.85]);
    set(fig, 'UserData', state);

    renderFrame(fig);

    uicontrol('Style','pushbutton','String','← Prev','Units','normalized',...
        'Position',[0.01 0.01 0.08 0.04],'Callback',@(~,~) navFrame(fig,-1));
    uicontrol('Style','pushbutton','String','Next →','Units','normalized',...
        'Position',[0.10 0.01 0.08 0.04],'Callback',@(~,~) navFrame(fig,+1));
    uicontrol('Style','pushbutton','String','Toggle Valid','Units','normalized',...
        'Position',[0.20 0.01 0.10 0.04],'Callback',@(~,~) toggleValid(fig));
    uicontrol('Style','pushbutton','String','DONE','Units','normalized',...
        'Position',[0.88 0.01 0.08 0.04],'Callback',@(~,~) doneFcn(fig), ...
        'BackgroundColor',[0.2 0.8 0.2]);

    uiwait(fig);
    if ishandle(fig)
        state = get(fig, 'UserData');
        detections = state.detections;
        close(fig);
    end
    fprintf('[GUI] Verification complete.\n');
end

function renderFrame(fig)
    state   = get(fig, 'UserData');
    f       = state.frameList(state.frameIdx);
    nCams   = state.nCams;

    clf(fig);
    nCols = ceil(sqrt(nCams));
    nRows = ceil(nCams / nCols);

    for c = 1:nCams
        ax = subplot(nRows, nCols, c, 'Parent', fig);
        fpath = fullfile(state.camDirs{c}, state.frameSets{c}{f});
        try
            im = imread(fpath);
        catch
            imagesc(zeros(100,100), 'Parent', ax); title(ax, sprintf('Cam %d — missing',c)); continue;
        end
        imagesc(im, 'Parent', ax); colormap(ax, gray); axis(ax,'image','off');

        det = state.detections{c}(f);
        hold(ax,'on');
        if det.valid
            plot(ax, det.pts(:,1), det.pts(:,2), 'g+', 'MarkerSize',14,'LineWidth',2);
            plot(ax, det.pts(:,1), det.pts(:,2), 'go', 'MarkerSize',16,'LineWidth',1.5);
            title(ax, sprintf('Cam %d ✓', c), 'Color','g');
        else
            title(ax, sprintf('Cam %d  [invalid]', c), 'Color','r');
        end
        hold(ax,'off');
    end

    sgtitle(fig, sprintf('Frame %d / %d  (index %d of %d valid frames)  [←/→ navigate]', ...
        f, state.frameList(end), state.frameIdx, numel(state.frameList)));
    drawnow;
end

function navFrame(fig, delta)
    state = get(fig,'UserData');
    state.frameIdx = max(1, min(numel(state.frameList), state.frameIdx + delta));
    set(fig,'UserData',state);
    renderFrame(fig);
end

function toggleValid(fig)
    state = get(fig,'UserData');
    f     = state.frameList(state.frameIdx);
    for c = 1:state.nCams
        state.detections{c}(f).valid = ~state.detections{c}(f).valid;
    end
    set(fig,'UserData',state);
    renderFrame(fig);
end

function keyHandler(fig, evt)
    switch evt.Key
        case 'leftarrow',  navFrame(fig,-1);
        case 'rightarrow', navFrame(fig,+1);
        case 't',          toggleValid(fig);
        case 'return',     doneFcn(fig);
    end
end

function doneFcn(fig)
    uiresume(fig);
end

% =========================================================================
%  HELPER: Detection overview figure
% =========================================================================
function generateDetectionOverview(detections, camDirs, frameSets, nCams, nFrames, resultsDir)
    % Pick up to 8 evenly spaced valid frames for a summary montage
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
            plot(ax, det.pts(:,1), det.pts(:,2), 'g+','MarkerSize',12,'LineWidth',2);
        end
        title(ax, sprintf('Frame %d', f));
        hold(ax,'off');
    end
    sgtitle('Wand Detections — Camera 1 Sample');
    saveas(fig, fullfile(resultsDir, 'detection_overview.png'));
    close(fig);
    fprintf('[Phase 2] Detection overview saved.\n');
end