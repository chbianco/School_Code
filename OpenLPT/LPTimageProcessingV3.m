%% LPTimageProcessing_simple.m
%
% A SIMPLE, LESS-AGGRESSIVE preprocessing pipeline for
% multi-camera particle-tracking image sequences, intended as a clean
% baseline to feed OpenLPT/STB.
%
% Design philosophy
% -----------------
% Do only what STB cannot do for itself, and make no irreversible,
% single-camera keep/drop decision beyond one conservative threshold.
% Concretely, relative to the earlier pipeline this script:

%       Uses a spatial SPATIAL BANDPASS filter. 
%       Glare off the flume floor is spatially extended,
%       low-frequency, and slowly wandering, so a per-pixel stationarity
%       test misses it whenever it fluctuates enough that no single pixel
%       clears the occupancy fraction (this is what leaked through on cam2).
%       A spatial bandpass removes contamination by SCALE, not by temporal
%       persistence, so it catches glare regardless of whether it holds
%       still. The temporal median still handles fixed-pattern noise and
%       hot pixels, so the occupancy filter's only remaining job is already
%       covered and it can go.
%
%       DEFAULTS to writing the CLEANED REAL IMAGE rather than re-rendering
%       every particle as an identical Gaussian. Re-rendering discards the
%       per-particle intensity and PSF variation that STB's OTF and shake-
%       residual machinery use to disambiguate, and it makes every 2D
%       keep/drop decision final. Feeding cleaned real images lets OpenLPT
%       detect on true intensities with nothing irreversible upstream. A
%       'gaussian' render mode is retained (opt-in) but note that even
%       in 'gaussian' mode this defaults to the MEASURED peak intensity,
%       not a uniform value, so intensity information survives.
%
%
% Expects a parent folder (selected via dialog) containing one subfolder
% per camera: cam0/, cam1/, ..., camN-1/, each full of .tif frames.
%
% Requires: Image Processing Toolbox (imregionalmax, ordfilt2, imgaussfilt)
% Parallel Computing Toolbox optional (parfor falls back to serial).
%
% Author: Christopher Bianco

clear; clc; close all;

%% ======================= USER PARAMETERS =======================

% --- I/O ---
outputDir_folder = 'imgFile';   % parent folder to write processed frames to
filePattern      = '*.tif';     % pattern to match input files
num_cams         = 4;           % number of cameras (expects cam0 ... camN-1)

% --- Temporal background (removes stationary structure + hot pixels) ---
bgMethod = 'median';   % 'median' (robust, recommended) | 'mean' | 'none'

% --- Spatial bandpass (removes smooth glare, per frame, by scale) ---
useSpatialBandpass = true;   % subtract a large-scale background estimate so
                             % smooth glare is removed whether or not it is
                             % stationary. Complementary to the temporal median
sigmaBG = 15;                % [px] std of the Gaussian used to ESTIMATE the
                             % smooth background. MUST be several times the
                             % particle radius, otherwise you subtract a
                             % blurred copy of the particles themselves and
                             % attenuate real signal. With ~5 px particles
                             % (sigma~2), 12-20 px is a safe range.
useNoiseSmoothing = true;    % gently smooth at the particle scale to raise
                             % detection SNR and suppress single-pixel shot
                             % noise. Set NEAR the particle sigma so it acts as
                             % a matched filter: averages down the fine 1-px
                             % noise carpet while preserving particle-sized
                             % blobs, lifting dim particles above noise. Do not
                             % exceed the particle sigma or neighbours blur.
sigmaNoise = 1.5;            % [px] std of the noise-smoothing Gaussian.
                             % should be less than the particle size 
                          

% --- Detection  ---
thresholdK     = 0.75;     % threshold = thresholdK * (robust noise sigma).
                        % 0.75 worked well for current setup, 4-6 might be
                        % a better goal if we get more lifhting 
absMinThreshold = 15;   % absolute detection-threshold floor, in input-count
                        % units, applied as thresh = max(thresholdK*sigma,
                        % absMinThreshold). This stops a low-noise camera from
                        % detecting arbitrarily deep into faint speckle: with a
                        % global MAD-based threshold, a camera whose quiet
                        % regions dominate the noise estimate gets a
                        % proportionally tiny threshold and lights up on residual
                        % floor texture. Set it ABOVE your speckle level and
                        % BELOW your dim-particle peaks--12-25 is the useful band. 
                        % Watch the per-camera printout to see where it binds, 
                        % and raise it until per-camera matches are smilar 
                        % Set to 0 to disable. SCALE this for non-8-bit data.
                       
localMaxWindow = 3;     % [px] neighborhood defining a local max (odd, >=3)
minPeakSep     = 2;     % [px] min separation between kept peaks (NMS radius)
subpixelMethod = 'gauss3pt';   % 'gauss3pt' (low bias, recommended) | 'none'

% --- Output mode ---
renderMode = 'cleaned';   % 'cleaned'  -> write the cleaned real image
                          %               (RECOMMENDED: least aggressive,
                          %                preserves true intensity/PSF, lets
                          %                OpenLPT do its own detection)
                          % 'gaussian' -> re-render detected centroids as
                          %               Gaussian blobs. Probably too
                          %               agressive

% Parameters used only when renderMode = 'gaussian':
peakIntensity   = 'measured';  % 'measured' -> use each particle's detected
                               %               peak (preserves intensity
                               %               variation for STB). Or give a
                               %               NUMBER for a fixed amplitude --
                               %               if you do, scale it to your
                               %               output bit depth (e.g. ~0.4x
                               %               the max, NOT a raw 100 into a
                               %               16-bit image, which would be
                               %               nearly black and get
                               %               thresholded away).
sigma           = 2;           % [px] rendered Gaussian std
renderHalfWidth = ceil(3*sigma);

% Parameters used only when renderMode = 'cleaned':
cleanedGain = 4;   % fallback gain, used ONLY if brightRef cannot be estimated.
                   % The cleaned output now uses a PER-CAMERA gain (camGain,
                   % computed below) instead of this fixed value.
noiseFloorK  = 2.5;  % subtract this many robust-noise-sigmas as a soft floor
                     % before gain -> removes the sub-noise "carpet" per camera
                     % while leaving real particles (well above noise) intact.
targetBright = 180;  % per-camera gain maps the typical bright particle (~p90 of
                     % detected peaks) to this level, equalizing brightness
                     % across cameras and avoiding 255 saturation.

% --- Output bit depth ---
outputBitDepth = 'same';   % 'same' | 'uint8' | 'uint16' | 'double'

% --- Misc ---
saveCentroids = false;   % also save detected sub-pixel centroids per frame
                         % (useful to cross-check counts vs OpenLPT's own 2D
                         % detections). Detection always runs for the console
                         % count regardless of renderMode.

% --- Image list naming (must match saved filenames; drop-in for OpenLPT) ---
listNamePrefix = 'img';
listNumDigits  = 6;
listExt        = '.tif';

% --- Parallel ---
useParallel = true;
numWorkers  = [];   % [] = default pool size

%% ======================= SETUP =======================

mainDir = uigetdir(pwd, 'Select parent OpenLPT folder');
if mainDir ~= 0
    cd(mainDir);
    fprintf('Working directory changed to:\n%s\n', mainDir);
else
    disp('Folder selection cancelled. Working directory unchanged.');
end

if useParallel && isempty(gcp('nocreate'))
    if isempty(numWorkers), parpool('local'); else, parpool('local', numWorkers); end
end

%% ======================= PER-CAMERA PROCESSING =======================

for n = 0:num_cams - 1

    cam           = ['cam', num2str(n)];
    cam_outputDir = fullfile(outputDir_folder, cam);
    if ~exist(cam_outputDir, 'dir'), mkdir(cam_outputDir); end
    inputDir = cam;

    fileList = getSortedFileList(inputDir, filePattern);
    nFrames  = numel(fileList);
    if nFrames == 0
        error('No files matching "%s" found in "%s".', filePattern, inputDir);
    end
    fprintf('\n=== %s: found %d frames ===\n', cam, nFrames);

    % Size/class from first frame
    info0      = imfinfo(fullfile(inputDir, fileList(1).name));
    imgH       = info0.Height;
    imgW       = info0.Width;
    inputClass = class(imread(fullfile(inputDir, fileList(1).name)));

    if strcmp(outputBitDepth, 'same'), outputClass = inputClass;
    else, outputClass = outputBitDepth; end
    switch outputClass
        case 'uint8',  outMaxVal = 255;
        case 'uint16', outMaxVal = 65535;
        otherwise,     outMaxVal = 1;     % 'double' assumed normalized [0,1]
    end

    % --- Load stack in NATIVE class (low memory) for the temporal background ---
    fprintf('Loading stack (native class) for temporal background...\n');
    stack = zeros(imgH, imgW, nFrames, inputClass);
    parfor k = 1:nFrames
        stack(:,:,k) = imread(fullfile(inputDir, fileList(k).name)); %#ok<PFBNS>
    end

    fprintf('Computing temporal background (%s)...\n', bgMethod);
    switch bgMethod
        case 'median', background = double(median(stack, 3));
        case 'mean',   background = mean(double(stack), 3);
        case 'none',   background = zeros(imgH, imgW);
        otherwise,     error('Unknown bgMethod: %s', bgMethod);
    end

    % --- Global threshold from a SAMPLED, robust noise estimate ---
    % (no full-stack flatten). Sample evenly spaced frames, process them the
    % same way the detector will, and take a MAD-based sigma of the signed
    % residual. Particles are sparse outliers, so the MAD reflects the noise.
    fprintf('Estimating global threshold (sampled, robust)...\n');
    nSampFrames = min(40, nFrames);
    sampIdx     = round(linspace(1, nFrames, nSampFrames));
    pxPerFrame  = min(1e5, imgH*imgW);
    samples     = zeros(pxPerFrame, nSampFrames);
    peakSamples = [];   % detected peak intensities across sample frames (for gain)
    for j = 1:nSampFrames
        k   = sampIdx(j);
        res = double(stack(:,:,k)) - background;
        bp  = spatialBandpass(res, useSpatialBandpass, sigmaBG, ...
                              useNoiseSmoothing, sigmaNoise);   % signed
        idx = randperm(imgH*imgW, pxPerFrame);
        samples(:, j) = bp(idx);
        [~, pv] = detectTracers(bp, absMinThreshold, localMaxWindow, ...
                                minPeakSep, subpixelMethod);
        peakSamples = [peakSamples; pv]; %#ok<AGROW>
    end
    samples      = samples(:);
    medSamp      = median(samples);
    noiseSig     = 1.4826 * median(abs(samples - medSamp));   % robust sigma
    robustThresh = thresholdK * noiseSig;
    thresh       = max(robustThresh, absMinThreshold);        % absolute floor
    if thresh > robustThresh
        fprintf(['  robust noise sigma = %.4f -> k*sigma = %.4f; ', ...
                 'absolute floor (%.2f) binds -> threshold = %.4f\n'], ...
                noiseSig, robustThresh, absMinThreshold, thresh);
    else
        fprintf('  robust noise sigma = %.4f -> threshold = %.4f (floor %.2f inactive)\n', ...
                noiseSig, thresh, absMinThreshold);
    end

    % --- Per-camera bright reference + equalization gain (cleaned mode) ---
    % Map this camera's typical bright particle to targetBright so all cameras
    % end up at comparable intensity (fixes the weak-camera cross-view mismatch).
    if ~isempty(peakSamples)
        sp        = sort(peakSamples);
        brightRef = sp(max(1, round(0.90 * numel(sp))));   % ~p90, no toolbox needed
    else
        brightRef = targetBright / cleanedGain;            % fallback -> old fixed gain
    end
    camGain = targetBright / max(brightRef, eps);
    fprintf('  brightRef(p90) = %.2f -> per-camera gain = %.3f\n', brightRef, camGain);

    if saveCentroids, allCentroids = cell(nFrames, 1); end

    % --- Main pass: every frame is independent -> clean parfor ---
    fprintf('Processing frames (%s output)...\n', renderMode);
    parfor k = 1:nFrames
        res = double(stack(:,:,k)) - background; %#ok<PFBNS>
        bp  = spatialBandpass(res, useSpatialBandpass, sigmaBG, ...
                              useNoiseSmoothing, sigmaNoise);   % signed

        [centroids, vals] = detectTracers(bp, thresh, localMaxWindow, ...
                                          minPeakSep, subpixelMethod);

        if saveCentroids, allCentroids{k} = centroids; end

        switch renderMode
            case 'cleaned'
                outFrame = max(bp - noiseFloorK*noiseSig, 0) * camGain;
            case 'gaussian'
                if ischar(peakIntensity) || isstring(peakIntensity)
                    amp = vals;                       % 'measured'
                else
                    amp = peakIntensity * ones(size(vals));
                end
                outFrame = renderGaussianTracers(imgH, imgW, centroids, amp, ...
                                                 sigma, renderHalfWidth);
            otherwise
                error('Unknown renderMode: %s', renderMode);
        end

        outFrame = min(outFrame, outMaxVal);
        outFrame = castToClass(outFrame, outputClass);

        frameNum = sprintf('%0*d', listNumDigits, k-1);
        outName  = fullfile(cam_outputDir, [listNamePrefix, frameNum, listExt]);
        imwrite(outFrame, outName);

        if mod(k, 50) == 0 || k == nFrames
            fprintf('  frame %d/%d: %d tracers detected\n', k, nFrames, size(centroids,1));
        end
    end

    if saveCentroids
        save(fullfile(cam_outputDir, 'centroids.mat'), 'allCentroids', 'fileList');
    end

    % --- OpenLPT image list ---
    camFolder     = [outputDir_folder, '/', cam, '/'];
    imageListFile = fullfile(outputDir_folder, [cam, 'ImageNames.txt']);
    fid = fopen(imageListFile, 'w');
    if fid == -1, error('Could not open %s for writing.', imageListFile); end
    for k = 1:nFrames
        frameNum = sprintf('%0*d', listNumDigits, k-1);
        fprintf(fid, '%s\n', [camFolder, listNamePrefix, frameNum, listExt]);
    end
    fclose(fid);
    fprintf('Wrote image list to: %s\n', imageListFile);
    fprintf('Done with %s. Output in: %s\n', cam, cam_outputDir);

    clear stack   % free before next camera
end

fprintf('\nAll cameras complete.\n');


%% ======================= LOCAL FUNCTIONS =======================

function fileList = getSortedFileList(folder, pattern)
% dir() sorted in natural numeric order on the last digit-group in each name.
    fileList = dir(fullfile(folder, pattern));
    if isempty(fileList), return; end
    nums = zeros(numel(fileList), 1);
    for i = 1:numel(fileList)
        d = regexp(fileList(i).name, '\d+', 'match');
        if isempty(d), nums(i) = i; else, nums(i) = str2double(d{end}); end
    end
    [~, order] = sort(nums);
    fileList = fileList(order);
end


function bp = spatialBandpass(res, doBandpass, sigmaBG, doNoise, sigmaNoise)
% Per-frame spatial bandpass on a (median-subtracted) residual image.
%   - Subtracts a large-scale Gaussian background estimate to remove smooth
%     glare, regardless of whether that glare is temporally stationary.
%   - Optionally smooths at the particle scale to raise detection SNR.
% Returns the SIGNED high-pass result (no negative clipping) so that a
% robust noise sigma can be estimated from it upstream; callers clip where
% appropriate. sigmaBG must be >> particle radius or particles are attenuated.
    if doBandpass
        bg = imgaussfilt(res, sigmaBG, 'Padding', 'replicate');
        hp = res - bg;
    else
        hp = res;
    end
    if doNoise && sigmaNoise > 0
        hp = imgaussfilt(hp, sigmaNoise, 'Padding', 'replicate');
    end
    bp = hp;
end


function [centroids, keepVal] = detectTracers(img, thresh, winSize, minSep, subpixelMethod)
% Detect bright tracer points in a (signed) bandpassed image.
%   centroids - [N x 2] = [x, y] = [col, row], sub-pixel unless 'none'
%   keepVal   - [N x 1] peak intensity at each detection (row-aligned)
    [H, W] = size(img);

    localMaxMask = imregionalmax(img, 8);
    winMax       = ordfilt2(img, winSize^2, true(winSize), 'symmetric');
    localMaxMask = localMaxMask & (img >= winMax) & (img > thresh);

    [ys, xs] = find(localMaxMask);
    vals     = img(localMaxMask);
    if isempty(xs)
        centroids = zeros(0, 2); keepVal = zeros(0, 1); return;
    end

    % Greedy non-max suppression by descending intensity
    [~, order] = sort(vals, 'descend');
    xs = xs(order); ys = ys(order); vals = vals(order);

    excludeMask = false(H, W);
    keepX = zeros(numel(xs),1); keepY = zeros(numel(xs),1); keepValAll = zeros(numel(xs),1);
    nKept = 0;
    [dxg, dyg] = meshgrid(-minSep:minSep, -minSep:minSep);
    diskMask   = (dxg.^2 + dyg.^2) <= minSep^2;
    [ddx, ddy] = deal(dxg(diskMask), dyg(diskMask));

    for i = 1:numel(xs)
        x = xs(i); y = ys(i);
        if ~excludeMask(y, x)
            nKept = nKept + 1;
            keepX(nKept) = x; keepY(nKept) = y; keepValAll(nKept) = vals(i);
            xi = x + ddx; yi = y + ddy;
            valid = xi >= 1 & xi <= W & yi >= 1 & yi <= H;
            excludeMask(sub2ind([H,W], yi(valid), xi(valid))) = true;
        end
    end
    keepX = keepX(1:nKept); keepY = keepY(1:nKept); keepVal = keepValAll(1:nKept);

    switch subpixelMethod
        case 'none'
            centroids = [keepX, keepY];
        case 'gauss3pt'
            xSub = keepX; ySub = keepY;
            for i = 1:nKept
                x = keepX(i); y = keepY(i);
                if x > 1 && x < W && y > 1 && y < H
                    dx = gauss3ptOffset(img(y,x-1), img(y,x), img(y,x+1));
                    dy = gauss3ptOffset(img(y-1,x), img(y,x), img(y+1,x));
                    xSub(i) = x + dx; ySub(i) = y + dy;
                end
            end
            centroids = [xSub, ySub];
        otherwise
            error('Unknown subpixelMethod: %s', subpixelMethod);
    end
end


function offset = gauss3ptOffset(Im, I0, Ip)
% 3-point Gaussian (log-parabola) sub-pixel peak estimator along one axis.
    if Im <= 0 || I0 <= 0 || Ip <= 0, offset = 0; return; end
    lm = log(Im); l0 = log(I0); lp = log(Ip);
    denom = (lm - 2*l0 + lp);
    if denom == 0
        offset = 0;
    else
        offset = 0.5 * (lm - lp) / denom;
        offset = max(min(offset, 0.5), -0.5);
    end
end


function outFrame = renderGaussianTracers(H, W, centroids, amp, sigma, halfWidth)
% Render each centroid as an additive 2D Gaussian, preserving the sub-pixel
% location via a fractional kernel shift. amp is a per-particle amplitude
% vector (measured peak) or a constant-per-particle vector.
    outFrame = zeros(H, W);
    if isempty(centroids), return; end
    [wx, wy] = meshgrid(-halfWidth:halfWidth, -halfWidth:halfWidth);
    for i = 1:size(centroids, 1)
        xc = centroids(i,1); yc = centroids(i,2);
        xi = round(xc); yi = round(yc);
        fracX = xc - xi; fracY = yc - yi;
        xs = xi + (-halfWidth:halfWidth);
        ys = yi + (-halfWidth:halfWidth);
        validX = xs >= 1 & xs <= W; validY = ys >= 1 & ys <= H;
        if ~any(validX) || ~any(validY), continue; end
        kernel = exp(-((wx - fracX).^2 + (wy - fracY).^2) / (2*sigma^2));
        outFrame(ys(validY), xs(validX)) = outFrame(ys(validY), xs(validX)) + ...
            amp(i) * kernel(validY, validX);
    end
end


function out = castToClass(img, className)
    switch className
        case 'uint8',  out = uint8(round(img));
        case 'uint16', out = uint16(round(img));
        case 'double', out = img;
        otherwise,     error('Unsupported outputBitDepth: %s', className);
    end
end