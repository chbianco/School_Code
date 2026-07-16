%% estimateSTBparameters.m
%
% Estimate OpenLPT/STB parameters from processed (cleaned) image frames plus
% the physical facts of your setup.
%
% WHAT THIS CAN AND CANNOT DO -- read this first
% ----------------------------------------------
% The STB parameters split into two kinds, and only one is derivable from an
% image:
%
%   * SPACING is a density measurement. Count real particles N in a frame,
%     divide the true illuminated volume V, take the cube root. The frame
%     gives N; YOU supply V. This is the number that goes in
%     "Convergence Avg Spacing (vox)".
%
%   * SEARCH RADII depend on per-frame DISPLACEMENT, which a single frame
%     cannot provide (no motion in one frame). Displacement = U * dt / voxel
%     size, so this script takes U (velocity), fps, and voxelToMm as INPUTS
%     and computes it. It does not read motion from the image.
%
%   * TOLERANCES depend on calibration REPROJECTION error, also an input.
%
% So this is "measure N and particle size from the frame, take the physics as
% inputs, turn both into recommended parameters with sanity ratios."
%
% CAVEATS THAT MATTER
%   - Use CLEAN frames (properly background-subtracted AND spatially
%     bandpassed, cam3 speckle and cam1 band removed). Contamination inflates
%     N and shrinks the estimated spacing. The script reports N variability
%     across frames as a contamination flag.
%   - Use the REAL illuminated slab for V, not an inflated reconstruction
%     volume. Spacing ~ V^(1/3), so even a rough slab depth is fine (2x error
%     in V -> only 26% in spacing), but a 6x-too-deep V matters.
%   - The gold-standard spacing is the mean nearest-neighbor distance of a
%     3D RECONSTRUCTED field. This density estimate is the best you can do
%     BEFORE a trustworthy reconstruction exists; once you have one, measure
%     spacing from it directly and prefer that.
%   - U should be the FASTEST velocity you intend to track (freestream), since
%     that sets the largest displacement and therefore the radii.
%
% Author: Christopher Bianco

clear; clc;

%% ======================= USER INPUTS =======================

% --- Frames to analyze (clean, single camera, the ones fed to OpenLPT) ---
frameGlob    = 'imgFile/cam0/img*.tif';
nFramesToUse = 20;     % average N over this many evenly spaced frames

% --- Acquisition / calibration ---
voxelToMm   = 0.27;    % [mm/vox] your reconstruction grid scale
fps         = 800;     % [1/s] frame rate
U_track     = 0.50;    % [m/s] FASTEST velocity to track (freestream). This
                       % sets the largest per-frame displacement.
pxToMm      = 0.27;    % [mm/px] object-space size of one IMAGE pixel. Only
                       % used to cross-check the density spacing against the
                       % measured 2D nearest-neighbor distance. If unknown,
                       % leaving it equal to voxelToMm assumes voxel ~ pixel.

% --- Real illuminated measurement volume [mm] (NOT the inflated bounds) ---
Lx_mm = 270;           % streamwise / X extent of the seeded, imaged slab
Ly_mm = 320;           % Y extent
Lz_mm = 200;            % <-- DEPTH of the illuminated slab. This is the
                       %     uncertain one; set it to the real lit thickness,
                       %     not the OpenLPT view-volume Z range.

% --- Calibration reprojection RMS error per camera [px] ---
reprojRMS_px = [0.9 0.9 0.9 1.4];

% --- Detection (match your preprocessing so N is consistent) ---
countThreshold = 25;   % same as absMinThreshold in the preprocessing script
localMaxWindow = 3;    % [px] local-max neighborhood (odd, >=3)
minPeakSep     = 2;    % [px] non-max suppression radius

% --- Optional reconstruction efficiency (fraction of real particles that
%     actually reconstruct in 3D). 1.0 = assume every imaged particle
%     reconstructs. Lowering it raises the implied density / lowers spacing.
reconEfficiency = 1.0;

% --- Recommendation multipliers (documented, tune to taste) ---
k_2Dtol     = 1.8;     % 2D tolerance = k_2Dtol * worst-camera reproj error
predictMult = 2.0;     % predict-field search radius = predictMult * displacement
initMult    = 3.0;     % initial-phase search radius = initMult * displacement

%% ======================= DETECT + MEASURE =======================
mainDir = uigetdir('C:\Users\FlumePIV\Desktop\', 'Select Your Project Folder');
cd(mainDir);

files = dir(frameGlob);
if isempty(files)
    error('No frames matched "%s".', frameGlob);
end
[~, ord] = sort({files.name});
files = files(ord);
folder = files(1).folder;

useIdx = unique(round(linspace(1, numel(files), min(nFramesToUse, numel(files)))));
counts = zeros(numel(useIdx), 1);

diamSamples = [];      % particle-image sigma samples (from the first frame)
nnSample    = [];      % 2D nearest-neighbor distances (from the first frame)

fprintf('Analyzing %d frames from %s\n', numel(useIdx), folder);
for j = 1:numel(useIdx)
    img = double(imread(fullfile(folder, files(useIdx(j)).name)));
    c   = detectTracers(img, countThreshold, localMaxWindow, minPeakSep);
    counts(j) = size(c, 1);

    if j == 1 && ~isempty(c)
        nnSample    = nearestNeighborDist(c);
        diamSamples = particleSigmas(img, c, 6);   % window radius 6 px
    end
end

Nmean = mean(counts) * reconEfficiency;
Nstd  = std(counts);
Ncv   = Nstd / mean(counts);          % coefficient of variation (contamination flag)

nn2D_med_px  = median(nnSample);
nn2D_med_vox = nn2D_med_px * pxToMm / voxelToMm;

sigma_px = median(diamSamples);
fwhm_px  = 2.3548 * sigma_px;         % particle image diameter (FWHM)

%% ======================= DERIVE PARAMETERS =======================

V_mm3   = Lx_mm * Ly_mm * Lz_mm;
n_mm3   = Nmean / V_mm3;               % number density [particles/mm^3]
s_mm    = (V_mm3 / Nmean)^(1/3);       % mean spacing [mm]
s_vox   = s_mm / voxelToMm;            % mean spacing [vox]  <-- the answer

% Per-frame displacement (physics, not image)
dt        = 1 / fps;                   % [s]
disp_mm   = U_track * dt * 1000;       % [mm]
disp_vox  = disp_mm / voxelToMm;       % [vox]

% Recommendations
predictR_vox = max(2, predictMult * disp_vox);
initR_vox    = min(initMult * disp_vox, 0.7 * s_vox);
initR_vox    = max(initR_vox, 3);
tol2D_px     = k_2Dtol * max(reprojRMS_px);
tol3D_vox    = tol2D_px * (pxToMm / voxelToMm);   % heuristic: reproj err in vox

% Sanity ratios
trackRatio   = disp_vox / s_vox;       % want << 1 (comfortably < ~0.33)
predRatio    = predictR_vox / s_vox;
initRatio    = initR_vox / s_vox;

%% ======================= REPORT =======================

fprintf('\n================ MEASURED FROM FRAMES ================\n');
fprintf('  particles/frame N  : %.0f  (std %.0f, CV %.2f)\n', Nmean, Nstd, Ncv);
fprintf('  particle image     : sigma %.2f px  ->  FWHM diam %.2f px\n', sigma_px, fwhm_px);
fprintf('  2D nearest-neighbor: median %.1f px  (~%.1f vox at pxToMm=%.3f)\n', ...
        nn2D_med_px, nn2D_med_vox, pxToMm);

fprintf('\n================ SUPPLIED PHYSICS ====================\n');
fprintf('  slab volume        : %.0f x %.0f x %.0f mm = %.3g mm^3\n', Lx_mm, Ly_mm, Lz_mm, V_mm3);
fprintf('  voxelToMm          : %.4f mm/vox\n', voxelToMm);
fprintf('  fps / U_track      : %.0f / %.3f m/s  ->  displacement %.2f vox/frame\n', ...
        fps, U_track, disp_vox);

fprintf('\n================ DERIVED =============================\n');
fprintf('  number density     : %.4g /mm^3\n', n_mm3);
fprintf('  mean spacing s     : %.2f mm  =  %.1f vox   (density route)\n', s_mm, s_vox);
fprintf('  cross-check (2D NN): %.1f vox  -- a LOWER BOUND on 3D spacing\n', nn2D_med_vox);
if nn2D_med_vox > s_vox
    fprintf('  [!] 2D-NN exceeds the density spacing -> V or N is off, or field is contaminated.\n');
end

fprintf('\n================ RECOMMENDED STB PARAMETERS ==========\n');
fprintf('  Convergence Avg Spacing (vox) : %.1f      [= density spacing s]\n', s_vox);
fprintf('  Predict Field Search Radius   : %.1f vox  [= %.1f x displacement]\n', predictR_vox, predictMult);
fprintf('  Initial Phase Search Radius   : %.1f vox  [~ %.1f x displacement, capped < s]\n', initR_vox, initMult);
fprintf('  2D Tolerance (px)             : %.2f     [= %.1f x worst reproj %.2f, cam%d]\n', ...
        tol2D_px, k_2Dtol, max(reprojRMS_px), find(reprojRMS_px==max(reprojRMS_px),1)-1);
fprintf('  3D Tolerance (vox)            : %.2f     [HEURISTIC -- refine vs IPR residuals]\n', tol3D_vox);

fprintf('\n================ SANITY CHECKS =======================\n');
printCheck('displacement / spacing', trackRatio, 0.33, ...
    'trackable: per-frame motion is a small fraction of spacing', ...
    'AMBIGUOUS: motion approaches spacing -> raise fps or lower seeding; no radius fixes this');
printCheck('predict radius / spacing', predRatio, 0.30, ...
    'predict radius safely below spacing', ...
    'predict radius too large vs spacing -> multiple candidates per sphere');
printCheck('initial radius / spacing', initRatio, 0.60, ...
    'initial radius below spacing', ...
    'initial radius near/over spacing -> ambiguous initialization');

if Ncv > 0.15
    fprintf('  [!] N varies %.0f%% across frames -- likely contamination (glare/speckle).\n', 100*Ncv);
    fprintf('      Clean the frames before trusting the spacing estimate.\n');
end

fprintf('\nNote: not shown (algorithmic defaults, not derivable from a frame):\n');
fprintf('  IPR Loops (~4), Reduced Loops (~2), Cameras to Reduce (1 = 3-of-4),\n');
fprintf('  Initial Phase Frames (~4), Shake Width (~0.4 vox), Shake Loops (~4),\n');
fprintf('  Ghost Threshold (revisit AFTER radii/spacing are correct).\n');
fprintf('  Predict Field Grid: set so each cell spans a few spacings (~%.0f vox).\n', 3*s_vox);


%% ======================= LOCAL FUNCTIONS =======================

function centroids = detectTracers(img, thresh, winSize, minSep)
% Local-maxima detection matching the preprocessing pipeline (integer peaks).
    [H, W] = size(img);
    lm  = imregionalmax(img, 8);
    wmx = ordfilt2(img, winSize^2, true(winSize), 'symmetric');
    lm  = lm & (img >= wmx) & (img > thresh);
    [ys, xs] = find(lm);
    vals = img(lm);
    if isempty(xs), centroids = zeros(0,2); return; end
    [~, o] = sort(vals, 'descend'); xs = xs(o); ys = ys(o);
    excl = false(H, W);
    keep = false(numel(xs),1);
    [dxg,dyg] = meshgrid(-minSep:minSep, -minSep:minSep);
    m = (dxg.^2+dyg.^2) <= minSep^2; ddx = dxg(m); ddy = dyg(m);
    for i = 1:numel(xs)
        if ~excl(ys(i), xs(i))
            keep(i) = true;
            xi = xs(i)+ddx; yi = ys(i)+ddy;
            v = xi>=1 & xi<=W & yi>=1 & yi<=H;
            excl(sub2ind([H,W], yi(v), xi(v))) = true;
        end
    end
    centroids = [xs(keep), ys(keep)];
end


function nn = nearestNeighborDist(pts)
% Brute-force 2D nearest-neighbor distance (no toolbox needed).
    n = size(pts,1); nn = zeros(n,1);
    for i = 1:n
        d2 = (pts(:,1)-pts(i,1)).^2 + (pts(:,2)-pts(i,2)).^2;
        d2(i) = inf;
        nn(i) = sqrt(min(d2));
    end
end


function sig = particleSigmas(img, pts, r)
% Intensity-weighted radial sigma of the brightest detections (particle size).
    [H, W] = size(img);
    vals = img(sub2ind([H,W], round(pts(:,2)), round(pts(:,1))));
    [~, o] = sort(vals, 'descend');
    take = o(1:min(200, numel(o)));
    sig = zeros(numel(take),1); k = 0;
    for t = take'
        x = round(pts(t,1)); y = round(pts(t,2));
        if x>r && x<=W-r && y>r && y<=H-r
            patch = img(y-r:y+r, x-r:x+r);
            patch = patch - min(patch(:)); patch(patch<0) = 0;
            [XX,YY] = meshgrid(-r:r, -r:r);
            tot = sum(patch(:));
            if tot > 0
                cx = sum(XX(:).*patch(:))/tot; cy = sum(YY(:).*patch(:))/tot;
                varr = sum(((XX(:)-cx).^2 + (YY(:)-cy).^2).*patch(:))/tot;
                k = k+1; sig(k) = sqrt(varr/2);   % per-axis sigma
            end
        end
    end
    sig = sig(1:k);
    if isempty(sig), sig = NaN; end
end


function printCheck(name, ratio, limit, okMsg, badMsg)
    if ratio <= limit
        fprintf('  [OK]   %-26s = %.2f  (<= %.2f)  %s\n', name, ratio, limit, okMsg);
    else
        fprintf('  [WARN] %-26s = %.2f  ( > %.2f)  %s\n', name, ratio, limit, badMsg);
    end
end