%% LPTimageProcessing.m
%
% Preprocessing pipeline for particle-tracking image sequences:
%   1. Background subtraction (stationary features -> 0)
%   2. Tracer (bright point) identification, with sub-pixel centroids
%   3. Assign each tracer a user-prescribed peak intensity
%   4. Re-render each tracer as a Gaussian blob of user-prescribed radius
%   5. Write the resulting images to disk
%
% Requires: Image Processing Toolbox (imregionalmax)
%
% Author: Christopher Bianco

clear; clc; close all;

%% ======================= USER PARAMETERS =======================

% --- I/O ---
outputDir_folder   = 'myProcessed';       % parent folder to write processed frames to
filePattern = '*.tif';               % pattern to match input files
num_cams = 4;   %Number of cameras 

% --- Background subtraction ---
bgMethod = 'median';   % 'median' | 'mean' | 'min'
                        % 'median' is recommended: robust to the tracers
                        % themselves, since any single pixel is only
                        % "hit" by a particle in a small fraction of frames.

% --- Tracer detection ---
threshScope  = 'global';   % 'global'   -> one threshold for the whole stack
                            %               (recommended: consistent particle
                            %               count/behavior frame-to-frame)
                            % 'perframe' -> recompute threshold each frame
threshMethod = 'meanstd';  % 'meanstd'    -> thresh = mean + threshParam*std
                            % 'percentile' -> thresh = prctile(data,threshParam)
                            % 'absolute'   -> thresh = threshParam (raw intensity)
threshParam  = 3;          % meaning depends on threshMethod (see above)

minPeakSep     = 2;     % [px] minimum allowed separation between two
                         % detected tracer centers (non-max suppression radius)
localMaxWindow = 3;     % [px] neighborhood size used to define a "local max"
                         % (must be odd, >=3)
subpixelRefine = true;  % true -> 3-point Gaussian estimator for sub-pixel
                         % centroid localization

% --- Output rendering ---
I_out           = 100;      % prescribed peak brightness of each rendered tracer
sigma           = 2;      % [px] Gaussian standard deviation ("radius")
renderHalfWidth = ceil(3*sigma);   % rendering window half-width (>=3 sigma
                                    % captures ~99% of the Gaussian energy)
outputBitDepth  = 'same';   % 'same' -> match input image class
                             % or explicitly: 'uint8' | 'uint16' | 'double'

% --- Misc ---
saveCentroids = false;   % if true, also saves detected (sub-pixel) centroid
                          % coordinates for each frame to a .mat file
                          % (useful for cross-checking counts against your
                          % downstream OpenLPT 2D detections)

%% ======================= SETUP =======================

mainDir = uigetdir(pwd, 'Select folder parent OpenLPT folder');

if mainDir ~= 0
    cd(mainDir);
    fprintf('Working directory successfully changed to:\n%s\n', mainDir);
else
    disp('Folder selection cancelled. Working directory remains unchanged.');
end

for n = 0:num_cams - 1

cam = ['cam', num2str(n)];
cam_outputDir = [outputDir_folder, '/', cam];
if ~exist(cam_outputDir , 'dir')
    mkdir(cam_outputDir);
end

inputDir = cam;

fileList = getSortedFileList(inputDir, filePattern);
nFrames  = numel(fileList);
if nFrames == 0
    error('No files matching "%s" found in "%s".', filePattern, inputDir);
end
fprintf('Found %d frames.\n', nFrames);

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

%% ======================= STEP 1: BACKGROUND =======================
% Load the full stack (as double, for accurate stats) and compute the
% per-pixel background statistic. Stationary features will subtract to
% (approximately) zero in every frame.

fprintf('Loading stack for background computation...\n');
stack = zeros(imgH, imgW, nFrames, 'double');
for k = 1:nFrames
    stack(:,:,k) = double(imread(fullfile(inputDir, fileList(k).name)));
end

fprintf('Computing background (%s projection)...\n', bgMethod);
switch bgMethod
    case 'median'
        background = median(stack, 3);
    case 'mean'
        background = mean(stack, 3);
    case 'min'
        background = min(stack, [], 3);
    otherwise
        error('Unknown bgMethod: %s', bgMethod);
end

% Background-subtracted stack (stationary features -> ~0, negatives clipped)
bgsubStack = stack - background;          % implicit expansion over dim 3
bgsubStack(bgsubStack < 0) = 0;

%% ======================= GLOBAL THRESHOLD (if requested) =======================

if strcmp(threshScope, 'global')
    globalThresh = computeThreshold(bgsubStack(:), threshMethod, threshParam);
    fprintf('Global detection threshold = %.3f\n', globalThresh);
end

%% ======================= STEPS 2-5: DETECT, RENDER, SAVE =======================

if saveCentroids
    allCentroids = cell(nFrames, 1);
end

fprintf('Processing frames...\n');
for k = 1:nFrames

    bgsub = bgsubStack(:,:,k);

    % --- threshold for this frame ---
    if strcmp(threshScope, 'perframe')
        thresh = computeThreshold(bgsub(:), threshMethod, threshParam);
    else
        thresh = globalThresh;
    end

    % --- Step 2: identify tracers (local maxima + NMS + subpixel) ---
    centroids = detectTracers(bgsub, thresh, localMaxWindow, minPeakSep, subpixelRefine);
    % centroids is [N x 2] with columns [x, y] (sub-pixel, image coordinates)

    if saveCentroids
        allCentroids{k} = centroids;
    end

    % --- Steps 3-4: assign fixed intensity + render as Gaussian blobs ---
    outFrame = renderGaussianTracers(imgH, imgW, centroids, I_out, sigma, renderHalfWidth);

    % --- clip to valid range and cast to output class ---
    outFrame = min(outFrame, outMaxVal);
    outFrame = castToClass(outFrame, outputClass);

    % --- Step 5: write output ---
    outName = fullfile(cam_outputDir, fileList(k).name);
    imwrite(outFrame, outName);

    if mod(k, 50) == 0 || k == nFrames
        fprintf('  frame %d/%d: %d tracers detected\n', k, nFrames, size(centroids,1));
    end
end

if saveCentroids
    save(fullfile(cam_outputDir, 'centroids.mat'), 'allCentroids', 'fileList');
end

fprintf('Done. Output written to: %s\n', cam_outputDir);

end

fprintf('All cameras complete')


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


function centroids = detectTracers(bgsub, thresh, winSize, minSep, doSubpixel)
% Detects bright tracer points in a background-subtracted image.
%   bgsub      - 2D background-subtracted image (double)
%   thresh     - intensity threshold
%   winSize    - odd integer, size of local-max neighborhood
%   minSep     - minimum allowed separation between kept peaks [px]
%   doSubpixel - if true, refine each centroid with a 3-point Gaussian fit
%
% Returns centroids as [N x 2] = [x, y] in image coordinates
% (x = column, y = row), sub-pixel if doSubpixel is true.

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
        return;
    end

    % --- greedy non-maximum suppression by descending intensity ---
    [~, order] = sort(vals, 'descend');
    xs = xs(order); ys = ys(order); vals = vals(order);

    excludeMask = false(H, W);
    keepX = zeros(numel(xs), 1);
    keepY = zeros(numel(xs), 1);
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

            % mark exclusion disk around this peak
            xi = x + ddx; yi = y + ddy;
            valid = xi >= 1 & xi <= W & yi >= 1 & yi <= H;
            lin = sub2ind([H, W], yi(valid), xi(valid));
            excludeMask(lin) = true;
        end
    end
    keepX = keepX(1:nKept);
    keepY = keepY(1:nKept);

    % --- sub-pixel refinement (3-point Gaussian estimator) ---
    if doSubpixel
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
    else
        centroids = [keepX, keepY];
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