% =========================================================================
% SyntheticDataGenerator.m
% =========================================================================
% PURPOSE:
%   Generate synthetic wand images + ground-truth data to validate the
%   entire calibration pipeline before running on real data.
%
%   This is the FIRST thing you should run to verify the pipeline works
%   end-to-end with known ground truth.
%
% WHAT IT DOES:
%   1. Generates random (but physically realistic) camera poses
%   2. Simulates a wand being waved through the measurement volume
%   3. Projects wand endpoints through the full refractive model
%   4. Adds Gaussian image noise (simulating real detection uncertainty)
%   5. Saves synthetic observations in the same format as Phase 2 output
%   6. Saves ground-truth poses for comparison after calibration
%
% USAGE:
%   1. Run this script
%   2. Run Phase4_BundleAdjustment.m  (using synthetic detections)
%   3. Run EvaluateCalibration.m      (compares result to ground truth)
% =========================================================================

clear; clc;
LPT_Config;

fprintf('=== Synthetic Data Generator ===\n\n');

rng(42);   % fixed seed for reproducibility

resultsDir = cfg.resultsDir;
nCams      = cfg.nCams;
wL         = cfg.wandLength;
n          = struct('air', cfg.n_air, 'glass', cfg.n_glass, 'water', cfg.n_water);

% -------------------------------------------------------------------------
% 1.  GENERATE GROUND-TRUTH CAMERA POSES
% -------------------------------------------------------------------------
% Place cameras around a measurement volume centred at origin.
% Volume assumed to be ~0.2 × 0.2 × 0.2 m inside the flume.

fprintf('Generating ground-truth camera poses...\n');

gt_cameras = struct();

% Camera positions (in air, outside the glass walls)
camPositions = [
    0.35,  0.00,  0.10;    % Cam 1 — +X side
    0.35,  0.05, -0.05;    % Cam 2 — +X side, offset
   -0.35,  0.00,  0.05;    % Cam 3 — -X side
    0.00,  0.00, -0.35;    % Cam 4 — bottom (floor)
];

% Aim each camera at the measurement volume centre + small random offset
volCentre = [0; 0; 0];
imageNoise_sigma = 1.0;   % px — simulated detection noise

for c = 1:nCams
    C  = camPositions(c,:)';
    aim = volCentre + 0.01*randn(3,1);  % slight aim error
    zAxis = aim - C; zAxis = zAxis/norm(zAxis);

    % Build orthonormal frame
    up = [0;0;1];
    if abs(dot(zAxis,up)) > 0.9, up = [0;1;0]; end
    xAxis = cross(zAxis, up); xAxis = xAxis/norm(xAxis);
    yAxis = cross(zAxis, xAxis);

    R  = [xAxis, yAxis, zAxis]';
    t  = -R * C;

    % Slightly randomise intrinsics
    f_x = 800 + 50*randn;
    f_y = f_x + 5*randn;
    cx  = cfg.imageSize(1)/2 + 10*randn;
    cy  = cfg.imageSize(2)/2 + 8*randn;
    K   = [f_x 0 cx; 0 f_y cy; 0 0 1];
    dc  = [0.05 -0.03 0.001 0.002 0.01] + 0.01*randn(1,5);

    gt_cameras(c).R          = R;
    gt_cameras(c).t          = t;
    gt_cameras(c).K          = K;
    gt_cameras(c).distCoeffs = dc;
    gt_cameras(c).C          = C;
    gt_cameras(c).wall       = cfg.walls(c);

    fprintf('  Camera %d: C=[%.3f %.3f %.3f]\n', c, C(1), C(2), C(3));
end

% -------------------------------------------------------------------------
% 2.  SIMULATE WAND TRAJECTORIES
% -------------------------------------------------------------------------
nFrames  = 100;
fprintf('\nSimulating %d wand frames...\n', nFrames);

gt_wand = cell(nFrames, 1);

for f = 1:nFrames
    % Random midpoint inside the measurement volume
    mid = 0.08 * (2*rand(3,1) - 1);

    % Random wand orientation
    phi   = 2*pi * rand;
    theta = pi * rand;
    wdir  = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];

    gt_wand{f} = struct('midpoint', mid, 'direction', wdir);
end

% -------------------------------------------------------------------------
% 3.  PROJECT WAND POINTS THROUGH REFRACTIVE MODEL
% -------------------------------------------------------------------------
fprintf('Projecting through refractive model + adding noise...\n');

detections = cell(nCams,1);
for c = 1:nCams
    detections{c} = repmat(struct('pts',[],'valid',false,'frame',0), nFrames,1);
end

nInvisible = 0;
for f = 1:nFrames
    mid  = gt_wand{f}.midpoint;
    wdir = gt_wand{f}.direction;
    XA   = mid - (wL/2)*wdir;
    XB   = mid + (wL/2)*wdir;

    for c = 1:nCams
        pose = struct('R', gt_cameras(c).R, 't', gt_cameras(c).t);
        K    = gt_cameras(c).K;
        dc   = gt_cameras(c).distCoeffs;
        wg   = struct('normal',    gt_cameras(c).wall.normal(:), ...
                       'point',     gt_cameras(c).wall.point(:), ...
                       'thickness', gt_cameras(c).wall.thickness);

        uvA = projectPointRefractive(XA, pose, K, dc, wg, n);
        uvB = projectPointRefractive(XB, pose, K, dc, wg, n);

        % Check visibility (within image bounds + in front of camera)
        imgW = cfg.imageSize(1); imgH = cfg.imageSize(2);
        inBoundsA = ~any(isnan(uvA)) && uvA(1)>=1 && uvA(1)<=imgW && uvA(2)>=1 && uvA(2)<=imgH;
        inBoundsB = ~any(isnan(uvB)) && uvB(1)>=1 && uvB(1)<=imgW && uvB(2)>=1 && uvB(2)<=imgH;
        
        % Randomly drop some detections to simulate occlusion
        dropProb = 0.15;  % 15% chance of missing detection per camera per frame
        isDropped = rand < dropProb;

        if inBoundsA && inBoundsB && ~isDropped
            uvA_noisy = uvA + imageNoise_sigma * randn(2,1);
            uvB_noisy = uvB + imageNoise_sigma * randn(2,1);

            detections{c}(f).pts   = [uvA_noisy'; uvB_noisy'];
            detections{c}(f).valid = true;
            detections{c}(f).frame = f;
        else
            nInvisible = nInvisible + 1;
        end
    end
end

fprintf('  Invisible detections: %d / %d\n', nInvisible, nFrames*nCams);

% -------------------------------------------------------------------------
% 4.  GENERATE SYNTHETIC INTRINSICS (as if from Phase 1)
% -------------------------------------------------------------------------
fprintf('\nGenerating synthetic intrinsics...\n');
intrinsics = cell(nCams,1);
calib      = cell(nCams,1);

for c = 1:nCams
    K  = gt_cameras(c).K;
    dc = gt_cameras(c).distCoeffs;

    % Create a cameraIntrinsics-like struct for compatibility
    % (Since we can't call estimateCameraParameters, we build the struct manually)
    intrinsics{c} = struct( ...
        'K',             K, ...
        'FocalLength',   [K(1,1), K(2,2)], ...
        'PrincipalPoint',[K(1,3), K(2,3)], ...
        'ImageSize',     cfg.imageSize);

    calib{c} = struct( ...
        'params',        struct('RadialDistortion', dc(1:3), ...
                                'TangentialDistortion', dc(4:5), ...
                                'Intrinsics', intrinsics{c}), ...
        'meanReprojErr', imageNoise_sigma * 0.8 + 0.05*rand);  % simulated error

    fprintf('  Camera %d: f=[%.1f %.1f] cx=%.1f cy=%.1f\n', c, ...
        K(1,1), K(2,2), K(1,3), K(2,3));
end

% -------------------------------------------------------------------------
% 5.  SAVE — in same format as Phase 1 & 2
% -------------------------------------------------------------------------
% Save as Phase 1 would have
save(fullfile(resultsDir, 'phase1_report.mat'), 'intrinsics', 'calib', 'cfg');

% Save per-camera intrinsics
for c = 1:nCams
    K          = gt_cameras(c).K;
    distCoeffs = gt_cameras(c).distCoeffs;
    imageSize  = cfg.imageSize;
    camParams  = calib{c}.params;
    save(fullfile(resultsDir, sprintf('intrinsics_cam%d.mat',c)), ...
         'K','distCoeffs','imageSize','camParams');
end

% Save as Phase 2 would have
nCams_save = nCams;
camDirs    = cell(nCams,1);
frameSets  = cell(nCams,1);
for c = 1:nCams
    camDirs{c}  = sprintf('synthetic_cam%d', c);
    frameSets{c} = arrayfun(@(i) sprintf('frame_%05d.png',i), 1:nFrames, 'UniformOutput',false);
end
save(fullfile(resultsDir, 'wand_detections.mat'), ...
     'detections', 'nFrames', 'nCams', 'camDirs', 'frameSets', 'cfg');

% Save ground truth for later evaluation
save(fullfile(resultsDir, 'ground_truth.mat'), 'gt_cameras', 'gt_wand', ...
     'imageNoise_sigma', 'nFrames', 'nCams', 'cfg');

fprintf('\n[SyntheticDataGenerator] Saved:\n');
fprintf('  phase1_report.mat    (synthetic intrinsics)\n');
fprintf('  wand_detections.mat  (synthetic observations)\n');
fprintf('  ground_truth.mat     (for evaluation)\n');
fprintf('\nNow run Phase4_BundleAdjustment.m, then EvaluateCalibration.m\n');
