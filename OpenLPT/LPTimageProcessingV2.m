%% LPTimageProcessing.m
%
% Preprocessing pipeline for multi-camera particle-tracking image sequences:
%   1. Background subtraction (stationary features -> 0)
%   2. Tracer (bright point) identification, with sub-pixel centroids and
%      temporal hysteresis thresholding (reduces frame-to-frame flicker)
%   3. Assign each tracer a user-prescribed peak intensity
%   4. Re-render each tracer as a Gaussian blob of user-prescribed radius
%   5. Write the resulting images + OpenLPT-style image list to disk
%
% Expects a parent folder (selected via dialog) containing one
% subfolder per camera: cam0/, cam1/, ..., camN-1/, each full of .tif
% frames. Output frames + image lists are written under <outputDir_folder>.
%
% Requires: Image Processing Toolbox (imregionalmax)
%
% Author: Christopher Bianco

clear; clc; close all;

%% ======================= USER PARAMETERS =======================

% --- I/O ---
outputDir_folder = 'imgFile';   % parent folder to write processed frames to
filePattern       = '*.tif';    % pattern to match input files
num_cams          = 4;          % number of cameras (expects cam0 ... camN-1)

% --- Background subtraction ---
bgMethod = 'median';   % 'median' | 'mean' | 'min' | 'running average'
                        % 'median' is recommended: robust to the tracers
                        % themselves, since any single pixel is only
                        % "hit" by a particle in a small fraction of frames.
                        % 'running average' should be better for moving
                        % particles that dwell near one location for a
                        % large fraction of the sequence.
runningAvgAlpha = 0.05; % learning rate for 'running average' (0 = static
                         % background, 1 = instant update). Only used if
                         % bgMethod = 'running average'.

% --- Tracer detection ---
threshScope  = 'global';   % 'global'   -> one threshold for the whole stack
                            %               (recommended: consistent particle
                            %               count/behavior frame-to-frame)
                            % 'perframe' -> recompute threshold each frame
threshMethod = 'meanstd';  % 'meanstd'    -> thresh = mean + threshParam*std
                            % 'percentile' -> thresh = prctile(data,threshParam)
                            % 'absolute'   -> thresh = threshParam (raw intensity)
threshParam  = 2.6;        % "confident" threshold parameter (see threshMethod)

useHysteresis  = true;   % apply temporal hysteresis thresholding: a weak
                          % detection is kept if it's spatially close to a
                          % CONFIDENT detection in a neighboring frame. This
                          % directly targets frame-to-frame "flicker" caused
                          % by a particle's peak intensity hovering near a
                          % single static threshold, where noise alone pushes
                          % it above/below the line from one frame to the next.
threshParamLow = 1.2;     % "weak"/candidate threshold parameter (same units
                          % as threshParam/threshMethod). Must be less strict
                          % than threshParam. If flicker persists even with
                          % hysteresis on, it usually means real dips are
                          % falling below threshParamLow entirely, so there's
                          % no candidate at all for hysteresis to rescue that
                          % frame. Check the "N candidates (M confident)"
                          % console output from Stage 1 -- if N barely
                          % exceeds M, threshParamLow is still too strict.
linkRadius     = 4;       % [px] max allowed frame-to-frame displacement for
                          % hysteresis linking. If this is smaller than your
                          % actual per-frame particle displacement (flow
                          % speed x dt x calibration), fast-moving particles
                          % will keep failing to link even when they're
                          % real. Check your real displacement/frame and set
                          % this with margin above it -- too large risks
                          % linking two different nearby particles in dense
                          % regions instead of confirming the same one.
hysteresisFrames = 2;    % how many neighboring frames on each side to check
                          % for a confident anchor. 2 bridges slightly longer
                          % dropouts (a particle dipping below threshold for
                          % 2 consecutive frames) than the previous default
                          % of 1. Larger values increase false-linking risk
                          % in dense regions -- don't push this much further
                          % without checking results.

rejectStaticArtifacts = true;  % flags and discards detections that sit at
                                % (nearly) the same location across an
                                % implausible fraction of frames. Real
                                % tracers shouldn't do this -- persistent
                                % fixed-location detections are almost always
                                % background-subtraction residue, debris on
                                % the glass, or a sensor hot pixel. NOTE:
                                % hysteresis linking makes these MORE likely
                                % to show up as apparently "solid" rather
                                % than flickering, since a non-moving
                                % artifact trivially satisfies the linkRadius
                                % proximity check every frame -- this filter
                                % specifically targets that failure mode.
staticRejectRadius  = 3;       % [px] radius used to group detections into
                                % the same "location" for the occupancy count.
                                % Raised from 1 -> 3: with radius=1, normal
                                % per-frame subpixel jitter (especially with
                                % gauss3pt, which trades some variance for
                                % lower bias) can split a truly-static
                                % artifact's occupancy across 2-3 neighboring
                                % pixels, so no single pixel individually
                                % crosses the flagging threshold even though
                                % the artifact is clearly present almost
                                % every frame. Confirmed via the peak-
                                % occupancy diagnostic: radius=1 kept known
                                % artifacts under threshold; radius=3 caught
                                % them cleanly (~755/800 frames, 94%).
staticOccupancyFrac = 0.8;     % if a location has a KEPT detection (confident
                                % or hysteresis-rescued -- i.e. whatever
                                % actually gets rendered) in more than this
                                % fraction of ALL frames, it's flagged as a
                                % static artifact and excluded everywhere.
                                % Lower this (e.g. 0.5-0.6) if artifacts still
                                % get through; raise it (e.g. 0.9+) if you
                                % have genuinely near-stagnant real particles
                                % (e.g. in recirculation zones) that you
                                % don't want to risk discarding.

minPeakSep     = 2;     % [px] minimum allowed separation between two
                         % detected tracer centers (non-max suppression radius)
localMaxWindow = 3;     % [px] neighborhood size used to define a "local max"
                         % (must be odd, >=3)

subpixelMethod = 'gauss3pt';  % 'gauss3pt'         -> 3-point Gaussian
                               %                       (log-parabola) peak
                               %                       estimator, applied
                               %                       separately along x
                               %                       and y using only the
                               %                       immediate neighbor
                               %                       pixels. Lower bias,
                               %                       higher variance --
                               %                       preferred for
                               %                       individual-particle
                               %                       tracking, since
                               %                       trajectory smoothing
                               %                       removes noise but
                               %                       not systematic bias.
                               % 'weightedCentroid' -> intensity-weighted
                               %                       center of mass over a
                               %                       window of radius
                               %                       centroidWindowRadius.
                               %                       Lower variance,
                               %                       higher bias
                               %                       (peak-locking) for
                               %                       small Gaussian-PSF
                               %                       particles like yours.
                               % 'none'             -> keep integer pixel
                               %                       location (no subpixel
                               %                       refinement)
                               % Use check_pixel_locking_bias.m and
                               % subpixel_estimator_validation.m (companion
                               % scripts) to verify which is actually best
                               % for your data.

centroidWindowRadius  = 2;     % [px] half-width of the window used by
                                % 'weightedCentroid' (only used if selected)
centroidSubtractFloor = true;  % if true, subtract the window's local minimum
                                % intensity before weighting (reduces bias
                                % from the non-zero pedestal under the tails)

% --- Output rendering ---
I_out           = 100;      % prescribed peak brightness of each rendered tracer
sigma           = 2;        % [px] Gaussian standard deviation ("radius")
renderHalfWidth = ceil(3*sigma);   % rendering window half-width (>=3 sigma
                                    % captures ~99% of the Gaussian energy)
outputBitDepth  = 'same';   % 'same' -> match input image class
                             % or explicitly: 'uint8' | 'uint16' | 'double'

% --- Misc ---
saveCentroids = false;   % if true, also saves detected (sub-pixel) centroid
                          % coordinates for each frame to a .mat file
                          % (useful for cross-checking counts against your
                          % downstream OpenLPT 2D detections)

% --- Image list file naming (per camera; e.g. for OpenLPT image input) ---
listNamePrefix = 'img';    % filename prefix
listNumDigits  = 6;        % zero-padding width
listExt        = '.tif';

%% ======================= SETUP: SELECT FOLDER =======================

mainDir = uigetdir(pwd, 'Select parent OpenLPT folder');

if mainDir ~= 0
    cd(mainDir);
    fprintf('Working directory successfully changed to:\n%s\n', mainDir);
else
    disp('Folder selection cancelled. Working directory remains unchanged.');
end

%% ======================= SETUP: PARALLEL POOL =======================

useParallel = true;   % if true, starts (or reuses) a parallel pool and the
                       % per-frame loops below run as `parfor`. Requires
                       % Parallel Computing Toolbox -- if that toolbox isn't
                       % installed, `parfor` automatically falls back to an
                       % ordinary serial loop, so this is safe to leave on
                       % either way.
numWorkers  = [];     % [] = MATLAB's default pool size (normally the number
                       % of physical cores). Set an integer to override --
                       % e.g. useful to leave headroom if you're running
                       % other memory-heavy processes (OpenLPT, etc.)
                       % alongside this.

if useParallel
    poolObj = gcp('nocreate');
    if isempty(poolObj)
        if isempty(numWorkers)
            parpool('local');
        else
            parpool('local', numWorkers);
        end
    end
end

%% ======================= PER-CAMERA PROCESSING =======================

for n = 0:num_cams - 1

    cam = ['cam', num2str(n)];
    cam_outputDir = [outputDir_folder, '/', cam];
    if ~exist(cam_outputDir, 'dir')
        mkdir(cam_outputDir);
    end

    inputDir = cam;

    fileList = getSortedFileList(inputDir, filePattern);
    nFrames  = numel(fileList);
    if nFrames == 0
        error('No files matching "%s" found in "%s".', filePattern, inputDir);
    end
    fprintf('\n=== %s: found %d frames ===\n', cam, nFrames);

    % Peek at the first image to get size/class
    info0      = imfinfo(fullfile(inputDir, fileList(1).name));
    imgH       = info0.Height;
    imgW       = info0.Width;
    firstImg   = imread(fullfile(inputDir, fileList(1).name));
    inputClass = class(firstImg);

    if strcmp(outputBitDepth, 'same')
        outputClass = inputClass;
    else
        outputClass = outputBitDepth;
    end

    switch outputClass
        case 'uint8',  outMaxVal = 255;
        case 'uint16', outMaxVal = 65535;
        otherwise,     outMaxVal = 1;   % 'double' -> assume normalized [0,1]
    end

    %% --- STEP 1: BACKGROUND ---
    fprintf('Loading stack for background computation...\n');
    stack = zeros(imgH, imgW, nFrames, 'double');
    parfor k = 1:nFrames
        stack(:,:,k) = double(imread(fullfile(inputDir, fileList(k).name))); %#ok<PFBNS>
    end

    fprintf('Computing background (%s projection)...\n', bgMethod);
    switch bgMethod
        case 'median'
            background = median(stack, 3);
        case 'mean'
            background = mean(stack, 3);
        case 'min'
            background = min(stack, [], 3);
        case 'running average'
            % Recursive filter -- inherently sequential, cannot be
            % parallelized across frames (each step depends on the last).
            % Typically cheap relative to detection, so this isn't a
            % priority target anyway.
            background = stack(:,:,1);   % initialize with first frame
            for f = 2:size(stack, 3)
                background = (1 - runningAvgAlpha) * background + runningAvgAlpha * stack(:,:,f);
            end
        otherwise
            error('Unknown bgMethod: %s', bgMethod);
    end

    % Background-subtracted stack (stationary features -> ~0, negatives clipped)
    bgsubStack = stack - background;          % implicit expansion over dim 3
    bgsubStack(bgsubStack < 0) = 0;
    clear stack   % no longer needed; frees ~imgH*imgW*nFrames*8 bytes before
                   % the (parallelized) detection passes below

    %% --- GLOBAL THRESHOLDS (if requested) ---
    if strcmp(threshScope, 'global')
        globalThreshHigh = computeThreshold(bgsubStack(:), threshMethod, threshParam);
        fprintf('Global confident threshold = %.3f\n', globalThreshHigh);
        if useHysteresis
            globalThreshLow = computeThreshold(bgsubStack(:), threshMethod, threshParamLow);
            fprintf('Global weak/candidate threshold = %.3f\n', globalThreshLow);
        end
    end

    %% --- STEPS 2-5: DETECT, RENDER, SAVE ---
    if saveCentroids
        allCentroids = cell(nFrames, 1);
    end

    if ~useHysteresis
        % ---------------- Single-pass pipeline (fully parallel) ----------------
        % Every frame is independent here (fixed threshold, no cross-frame
        % reads), so this is a clean parfor with no broadcast concerns.
        fprintf('Processing frames (no hysteresis)...\n');
        parfor k = 1:nFrames
            bgsub = bgsubStack(:,:,k); %#ok<PFBNS>

            if strcmp(threshScope, 'perframe')
                thresh = computeThreshold(bgsub(:), threshMethod, threshParam);
            else
                thresh = globalThreshHigh;
            end

            centroids = detectTracers(bgsub, thresh, localMaxWindow, minPeakSep, ...
                subpixelMethod, centroidWindowRadius, centroidSubtractFloor);

            if saveCentroids
                allCentroids{k} = centroids;
            end

            outFrame = renderGaussianTracers(imgH, imgW, centroids, I_out, sigma, renderHalfWidth);
            outFrame = min(outFrame, outMaxVal);
            outFrame = castToClass(outFrame, outputClass);

            outName = fullfile(cam_outputDir, fileList(k).name);
            imwrite(outFrame, outName);

            if mod(k, 50) == 0 || k == nFrames
                fprintf('  frame %d/%d: %d tracers detected\n', k, nFrames, size(centroids,1));
            end
        end

    else
        % ---------------- Three-stage temporal hysteresis pipeline ----------------
        % Split into: (1) parallel detection [expensive, frame-independent],
        % (2) sequential hysteresis linking [cheap, needs neighboring
        % frames' candidate lists], (3) parallel render+save [expensive,
        % frame-independent once stage 2 has resolved the final centroid
        % list per frame]. This avoids broadcasting the full candidate
        % lists to every worker -- only stage 2 (which is NOT parallelized)
        % ever reads across frames, so there's no large-array broadcast
        % cost for the parallel stages.

        % --- Stage 1/3 (parallel): detect every candidate at the LOW
        %     (permissive) threshold, keeping each candidate's peak
        %     intensity so we can classify it as confident/weak after. ---
        fprintf('Stage 1/3: detecting candidates at low threshold (parallel)...\n');
        candCentroids = cell(nFrames, 1);
        candVals      = cell(nFrames, 1);
        isConfident   = cell(nFrames, 1);

        parfor k = 1:nFrames
            bgsub = bgsubStack(:,:,k); %#ok<PFBNS>

            if strcmp(threshScope, 'perframe')
                threshLow  = computeThreshold(bgsub(:), threshMethod, threshParamLow);
                threshHigh = computeThreshold(bgsub(:), threshMethod, threshParam);
            else
                threshLow  = globalThreshLow;
                threshHigh = globalThreshHigh;
            end

            [c, v] = detectTracers(bgsub, threshLow, localMaxWindow, minPeakSep, ...
                subpixelMethod, centroidWindowRadius, centroidSubtractFloor);

            candCentroids{k} = c;
            candVals{k}      = v;
            isConfident{k}   = v >= threshHigh;

            if mod(k, 50) == 0 || k == nFrames
                fprintf('  frame %d/%d: %d candidates (%d confident)\n', ...
                    k, nFrames, numel(v), nnz(isConfident{k}));
            end
        end

        clear bgsubStack   % no longer needed; frees memory before stage 3

        % --- Stage 2a/3 (sequential, cheap): hysteresis linking only.
        %     Produces a PRELIMINARY per-frame kept set (confident +
        %     hysteresis-rescued weak candidates). Static-artifact
        %     detection below deliberately runs AFTER this, not before --
        %     see note. ---
        fprintf('Stage 2a/3: hysteresis linking (sequential)...\n');

        [dxg, dyg] = meshgrid(-linkRadius:linkRadius, -linkRadius:linkRadius);
        diskMask = (dxg.^2 + dyg.^2) <= linkRadius^2;
        [ddx, ddy] = deal(dxg(diskMask), dyg(diskMask));

        preliminaryKeep = cell(nFrames, 1);

        for k = 1:nFrames
            conf = isConfident{k};
            keep = conf;   % confident detections are always kept

            if any(~conf)
                c = candCentroids{k};
                % Build an occupancy mask from CONFIDENT detections in
                % neighboring frames, then confirm weak candidates that fall
                % within linkRadius of it.
                neighborMask = false(imgH, imgW);
                for offset = [-hysteresisFrames:-1, 1:hysteresisFrames]
                    kn = k + offset;
                    if kn >= 1 && kn <= nFrames
                        cn = candCentroids{kn}(isConfident{kn}, :);
                        for i = 1:size(cn, 1)
                            xi = round(cn(i,1)) + ddx;
                            yi = round(cn(i,2)) + ddy;
                            valid = xi >= 1 & xi <= imgW & yi >= 1 & yi <= imgH;
                            lin = sub2ind([imgH, imgW], yi(valid), xi(valid));
                            neighborMask(lin) = true;
                        end
                    end
                end

                weakIdx = find(~conf);
                for wi = weakIdx'
                    xi = round(c(wi,1)); yi = round(c(wi,2));
                    if xi >= 1 && xi <= imgW && yi >= 1 && yi <= imgH && neighborMask(yi, xi)
                        keep(wi) = true;
                    end
                end
            end

            preliminaryKeep{k} = keep;
        end

        % --- Stage 2b/3 (sequential, cheap): flag static (non-moving)
        %     false detections. IMPORTANT: this is built from the
        %     PRELIMINARY POST-HYSTERESIS kept set above, not from
        %     isConfident alone. A non-moving artifact trivially satisfies
        %     the linkRadius proximity check against its own neighboring
        %     frames every time, so hysteresis will rescue it on almost
        %     every weak frame -- meaning it can be rendered in ~100% of
        %     frames while only being "confident" a fraction of the time.
        %     Computing occupancy from isConfident alone would systematically
        %     UNDERCOUNT exactly these cases and let them slip through the
        %     filter. Using the post-hysteresis kept set closes that gap. ---
        if rejectStaticArtifacts
            fprintf('Stage 2b/3: flagging static (non-moving) false detections...\n');

            [dxs, dys] = meshgrid(-staticRejectRadius:staticRejectRadius, -staticRejectRadius:staticRejectRadius);
            staticDiskMask = (dxs.^2 + dys.^2) <= staticRejectRadius^2;
            [sddx, sddy] = deal(dxs(staticDiskMask), dys(staticDiskMask));

            occupancyCount = zeros(imgH, imgW, 'uint32');
            for k = 1:nFrames
                keptPts = candCentroids{k}(preliminaryKeep{k}, :);
                if isempty(keptPts)
                    continue;
                end
                frameMask = false(imgH, imgW);
                for i = 1:size(keptPts, 1)
                    xi = round(keptPts(i,1)) + sddx;
                    yi = round(keptPts(i,2)) + sddy;
                    valid = xi >= 1 & xi <= imgW & yi >= 1 & yi <= imgH;
                    lin = sub2ind([imgH, imgW], yi(valid), xi(valid));
                    frameMask(lin) = true;
                end
                occupancyCount = occupancyCount + uint32(frameMask);
            end

            % --- DIAGNOSTIC (temporary): report the single most-occupied
            %     pixel, even if it doesn't cross the flagging threshold.
            %     Helps distinguish "no artifacts present" (small max) from
            %     "occupancy is being split across neighboring pixels"
            %     (large max that still doesn't individually clear
            %     staticOccupancyFrac*nFrames). Safe to remove once resolved.
            [maxOcc, maxIdx] = max(occupancyCount(:));
            [maxY, maxX] = ind2sub([imgH, imgW], maxIdx);
            fprintf('  peak occupancy: %d/%d frames at (x=%d, y=%d)\n', maxOcc, nFrames, maxX, maxY);

            staticArtifactMask = occupancyCount > (staticOccupancyFrac * nFrames);
            rawFlaggedCount = nnz(staticArtifactMask);

            % Dilate by staticRejectRadius before using this mask for
            % rejection. occupancyCount was built by spreading each
            % detection across a radius-staticRejectRadius disk, so a
            % jittery-but-static artifact's occupancy is smeared across a
            % small neighborhood and only its highest-count pixel(s) cross
            % the threshold. Stage 2c below checks a detection's EXACT
            % rounded pixel against this mask with no radius tolerance --
            % without dilating first, positions belonging to the same
            % artifact that jitter a pixel or two away from that peak
            % would slip through untouched (confirmed empirically: the
            % exact flagged peak pixel dropped from ~94% to ~41% presence
            % after rejection, while an unflagged neighbor 1px away
            % stayed at ~69%).
            [flagY, flagX] = find(staticArtifactMask);
            dilatedMask = false(imgH, imgW);
            for i = 1:numel(flagX)
                xi = flagX(i) + sddx;
                yi = flagY(i) + sddy;
                valid = xi >= 1 & xi <= imgW & yi >= 1 & yi <= imgH;
                lin = sub2ind([imgH, imgW], yi(valid), xi(valid));
                dilatedMask(lin) = true;
            end
            staticArtifactMask = dilatedMask;

            fprintf('  flagged %d px as static artifacts (present in >%.0f%% of frames), %d after dilation\n', ...
                rawFlaggedCount, staticOccupancyFrac*100, nnz(staticArtifactMask));
        else
            staticArtifactMask = false(imgH, imgW);
        end

        % --- Stage 2c/3 (sequential, cheap): apply static-artifact
        %     rejection to the preliminary kept set -> finalCentroids. ---
        finalCentroids = cell(nFrames, 1);

        for k = 1:nFrames
            c = candCentroids{k};
            keep = preliminaryKeep{k};

            if rejectStaticArtifacts && any(keep)
                keepIdx = find(keep);
                for ki = keepIdx'
                    xi = round(c(ki,1)); yi = round(c(ki,2));
                    if xi >= 1 && xi <= imgW && yi >= 1 && yi <= imgH && staticArtifactMask(yi, xi)
                        keep(ki) = false;
                    end
                end
            end

            finalCentroids{k} = c(keep, :);
        end

        clear candCentroids candVals isConfident preliminaryKeep   % no longer needed

        % --- Stage 3/3 (parallel): render + save. Only reads finalCentroids{k}
        %     (sliced, one frame at a time) so this is a clean parfor. ---
        fprintf('Stage 3/3: rendering and saving (parallel)...\n');
        parfor k = 1:nFrames
            centroids = finalCentroids{k}; %#ok<PFBNS>

            if saveCentroids
                allCentroids{k} = centroids;
            end

            outFrame = renderGaussianTracers(imgH, imgW, centroids, I_out, sigma, renderHalfWidth);
            outFrame = min(outFrame, outMaxVal);
            outFrame = castToClass(outFrame, outputClass);

            outName = fullfile(cam_outputDir, fileList(k).name);
            imwrite(outFrame, outName);

            if mod(k, 50) == 0 || k == nFrames
                fprintf('  frame %d/%d: %d tracers rendered\n', k, nFrames, size(centroids,1));
            end
        end
    end

    if saveCentroids
        save(fullfile(cam_outputDir, 'centroids.mat'), 'allCentroids', 'fileList');
    end

    % --- Image list file (e.g. for OpenLPT camera image input) ---
    camFolder     = [outputDir_folder, '/', cam, '/'];   % virtual path prefix written into the list
    imageListFile = fullfile(outputDir_folder, [cam, 'ImageNames.txt']);

    fid = fopen(imageListFile, 'w');
    if fid == -1
        error('Could not open %s for writing.', imageListFile);
    end
    for k = 1:nFrames
        frameNum  = sprintf('%0*d', listNumDigits, k-1);          % 0, 1, 2, ...
        frameName = [camFolder, listNamePrefix, frameNum, listExt];
        fprintf(fid, '%s\n', frameName);
    end
    fclose(fid);
    fprintf('Wrote image list to: %s\n', imageListFile);

    fprintf('Done with %s. Output written to: %s\n', cam, cam_outputDir);

end

fprintf('\nAll cameras complete.\n');


%% ======================= LOCAL FUNCTIONS =======================

function fileList = getSortedFileList(folder, pattern)
% Returns dir() struct array sorted in natural numeric order based on any
% digits found in the filename (robust to inconsistent zero-padding).
    fileList = dir(fullfile(folder, pattern));
    if isempty(fileList)
        return;
    end
    nums = zeros(numel(fileList), 1);
    for i = 1:numel(fileList)
        d = regexp(fileList(i).name, '\d+', 'match');
        if isempty(d)
            nums(i) = i;   % no digits found: fall back to original order
        else
            nums(i) = str2double(d{end});   % use the last numeric group
        end
    end
    [~, order] = sort(nums);
    fileList = fileList(order);
end


function thresh = computeThreshold(data, method, param)
% Computes a scalar intensity threshold from a vector of pixel values.
    data = data(data > 0);   % ignore hard zeros (clipped background)
    switch method
        case 'meanstd'
            thresh = mean(data) + param * std(data);
        case 'percentile'
            thresh = simplePrctile(data, param);
        case 'absolute'
            thresh = param;
        otherwise
            error('Unknown threshMethod: %s', method);
    end
end


function p = simplePrctile(data, pct)
% Minimal percentile implementation (linear interpolation on sorted data),
% avoids a Statistics and Machine Learning Toolbox dependency.
    data = sort(data(:));
    n = numel(data);
    if n == 0
        p = 0;
        return;
    end
    ranks = (pct/100) * (n - 1) + 1;
    lo = floor(ranks); hi = ceil(ranks);
    lo = min(max(lo,1),n); hi = min(max(hi,1),n);
    w = ranks - lo;
    p = (1-w)*data(lo) + w*data(hi);
end


function [centroids, keepVal] = detectTracers(bgsub, thresh, winSize, minSep, ...
    subpixelMethod, centroidWindowRadius, centroidSubtractFloor)
% Detects bright tracer points in a background-subtracted image.
%   bgsub          - 2D background-subtracted image (double)
%   thresh         - intensity threshold
%   winSize        - odd integer, size of local-max neighborhood
%   minSep         - minimum allowed separation between kept peaks [px]
%   subpixelMethod - 'gauss3pt' | 'weightedCentroid' | 'none'
%   centroidWindowRadius  - window radius for 'weightedCentroid' [px]
%   centroidSubtractFloor - subtract local window minimum before weighting
%
% Returns:
%   centroids - [N x 2] = [x, y] in image coordinates (x = column, y = row),
%               sub-pixel if subpixelMethod ~= 'none'
%   keepVal   - [N x 1] peak (pre-subpixel) intensity at each detection,
%               row-aligned with centroids -- used upstream for hysteresis
%               confident/weak classification

    [H, W] = size(bgsub);

    % --- candidate local maxima ---
    localMaxMask = imregionalmax(bgsub, 8);
    % Additionally enforce that the "flat" regionalmax blobs are true
    % single-pixel peaks by intersecting with a moving-window max filter
    % (guards against imregionalmax returning multi-pixel plateaus)
    winMax = ordfilt2(bgsub, winSize^2, true(winSize), 'symmetric');
    localMaxMask = localMaxMask & (bgsub >= winMax) & (bgsub > thresh);

    [ys, xs] = find(localMaxMask);
    vals     = bgsub(localMaxMask);

    if isempty(xs)
        centroids = zeros(0, 2);
        keepVal = zeros(0, 1);
        return;
    end

    % --- greedy non-maximum suppression by descending intensity ---
    [~, order] = sort(vals, 'descend');
    xs = xs(order); ys = ys(order); vals = vals(order);

    excludeMask = false(H, W);
    keepX = zeros(numel(xs), 1);
    keepY = zeros(numel(xs), 1);
    keepValAll = zeros(numel(xs), 1);
    nKept = 0;

    % Precompute disk offsets for marking the exclusion zone
    [dxg, dyg] = meshgrid(-minSep:minSep, -minSep:minSep);
    diskMask = (dxg.^2 + dyg.^2) <= minSep^2;
    [ddx, ddy] = deal(dxg(diskMask), dyg(diskMask));

    for i = 1:numel(xs)
        x = xs(i); y = ys(i);
        if ~excludeMask(y, x)
            nKept = nKept + 1;
            keepX(nKept) = x;
            keepY(nKept) = y;
            keepValAll(nKept) = vals(i);

            % mark exclusion disk around this peak
            xi = x + ddx; yi = y + ddy;
            valid = xi >= 1 & xi <= W & yi >= 1 & yi <= H;
            lin = sub2ind([H, W], yi(valid), xi(valid));
            excludeMask(lin) = true;
        end
    end
    keepX = keepX(1:nKept);
    keepY = keepY(1:nKept);
    keepVal = keepValAll(1:nKept);

    % --- sub-pixel refinement ---
    switch subpixelMethod
        case 'none'
            centroids = [keepX, keepY];
            return;
        case 'gauss3pt'
            xSub = zeros(nKept, 1);
            ySub = zeros(nKept, 1);
            for i = 1:nKept
                x = keepX(i); y = keepY(i);
                if x > 1 && x < W && y > 1 && y < H
                    Ixm = bgsub(y, x-1); Ix0 = bgsub(y, x); Ixp = bgsub(y, x+1);
                    Iym = bgsub(y-1, x); Iy0 = bgsub(y, x); Iyp = bgsub(y+1, x);

                    dx = gauss3ptOffset(Ixm, Ix0, Ixp);
                    dy = gauss3ptOffset(Iym, Iy0, Iyp);

                    xSub(i) = x + dx;
                    ySub(i) = y + dy;
                else
                    xSub(i) = x;
                    ySub(i) = y;
                end
            end
            centroids = [xSub, ySub];
        case 'weightedCentroid'
            xSub = zeros(nKept, 1);
            ySub = zeros(nKept, 1);
            r = centroidWindowRadius;
            for i = 1:nKept
                x = keepX(i); y = keepY(i);
                [xc, yc] = weightedCentroidOffset(bgsub, x, y, r, W, H, centroidSubtractFloor);
                xSub(i) = xc;
                ySub(i) = yc;
            end
            centroids = [xSub, ySub];
        otherwise
            error('Unknown subpixelMethod: %s', subpixelMethod);
    end
end


function offset = gauss3ptOffset(Im, I0, Ip)
% 3-point Gaussian sub-pixel peak estimator along one axis.
% Im, I0, Ip = intensities at (peak-1), (peak), (peak+1).
% Falls back to 0 offset if inputs are non-positive or degenerate
% (avoids log(0) / division-by-zero).
    if Im <= 0 || I0 <= 0 || Ip <= 0
        offset = 0;
        return;
    end
    lm = log(Im); l0 = log(I0); lp = log(Ip);
    denom = (lm - 2*l0 + lp);
    if denom == 0
        offset = 0;
    else
        offset = 0.5 * (lm - lp) / denom;
        offset = max(min(offset, 0.5), -0.5);   % sanity clamp
    end
end


function [xc, yc] = weightedCentroidOffset(bgsub, x, y, r, W, H, subtractFloor)
% Intensity-weighted centroid (center of mass) computed over a
% (2r+1)x(2r+1) window centered on the integer peak (x,y).
%   subtractFloor - if true, subtract the window's local minimum before
%                   weighting. This reduces (but does not eliminate) the
%                   bias caused by the finite window truncating the tails
%                   of the particle's intensity profile asymmetrically.
    xlo = max(x-r, 1); xhi = min(x+r, W);
    ylo = max(y-r, 1); yhi = min(y+r, H);

    patch = bgsub(ylo:yhi, xlo:xhi);
    [XX, YY] = meshgrid(xlo:xhi, ylo:yhi);

    w = patch;
    if subtractFloor
        w = w - min(w(:));
    end
    w(w < 0) = 0;

    total = sum(w(:));
    if total <= 0
        xc = x; yc = y;   % degenerate patch: fall back to integer peak
        return;
    end

    xc = sum(XX(:) .* w(:)) / total;
    yc = sum(YY(:) .* w(:)) / total;
end


function outFrame = renderGaussianTracers(H, W, centroids, I0, sigma, halfWidth)
% Renders each tracer centroid as an additive 2D Gaussian blob of peak
% amplitude I0 and standard deviation sigma onto an (H x W) canvas.
    outFrame = zeros(H, W, 'double');
    if isempty(centroids)
        return;
    end

    [wx, wy] = meshgrid(-halfWidth:halfWidth, -halfWidth:halfWidth);

    for i = 1:size(centroids, 1)
        xc = centroids(i, 1);
        yc = centroids(i, 2);

        xi = round(xc); yi = round(yc);
        % sub-pixel remainder shifts the kernel grid, not just its amplitude
        fracX = xc - xi; fracY = yc - yi;

        xs = xi + (-halfWidth:halfWidth);
        ys = yi + (-halfWidth:halfWidth);

        validX = xs >= 1 & xs <= W;
        validY = ys >= 1 & ys <= H;

        if ~any(validX) || ~any(validY)
            continue;
        end

        % recompute kernel with true sub-pixel offset for this tracer
        localKernel = exp(-((wx - fracX).^2 + (wy - fracY).^2) / (2*sigma^2));

        outFrame(ys(validY), xs(validX)) = outFrame(ys(validY), xs(validX)) + ...
            I0 * localKernel(validY, validX);
    end
end


function out = castToClass(img, className)
    switch className
        case 'uint8'
            out = uint8(round(img));
        case 'uint16'
            out = uint16(round(img));
        case 'double'
            out = img;
        otherwise
            error('Unsupported outputBitDepth: %s', className);
    end
end