% =========================================================================
% Phase4_BundleAdjustment.m
% =========================================================================
% PURPOSE:
%   Given the Phase-1 intrinsics and Phase-2 wand detections, jointly
%   optimise:
%     • Camera extrinsic poses  (R_c, t_c)  × nCams
%     • 3-D wand positions      (mid-point + orientation) × nFrames
%     • (optionally) wall normal perturbations
%
%   The cost function is the sum of squared reprojection errors computed
%   through the full refractive ray model (Phase 3).
%
%   INITIALISATION STRATEGY:
%     1. Build an initial camera network using the wand constraint:
%        |P1 - P2|_world = wandLength  (scale constraint)
%     2. Use RANSAC + DLT (ignoring refraction) for a rough pose estimate
%     3. Refine with full refractive bundle adjustment via lsqnonlin
%
% OUTPUTS (saved to results/):
%   extrinsics.mat         — R, t per camera in world coords
%   bundle_adjustment.mat  — full BA result struct
%   ba_convergence.png     — cost vs iteration plot
% =========================================================================

clear; clc;
LPT_Config;

fprintf('=== Phase 4: Bundle Adjustment ===\n\n');

resultsDir = cfg.resultsDir;

% -------------------------------------------------------------------------
% Load prerequisite data
% -------------------------------------------------------------------------
phase1 = load(fullfile(resultsDir, 'phase1_report.mat'));
phase2 = load(fullfile(resultsDir, 'wand_detections.mat'));

intrinsics  = phase1.intrinsics;    % cell array of cameraIntrinsics objects
detections  = phase2.detections;
nCams       = cfg.nCams;
nFrames     = phase2.nFrames;

% Unpack intrinsics into K and distCoeffs cells
Ks          = cell(nCams,1);
distCoeffs  = cell(nCams,1);


for c = 1:nCams
    Ks{c}         = intrinsics{c}.K;
    dc            = [phase1.calib{c}.params.RadialDistortion, ...
                     phase1.calib{c}.params.TangentialDistortion];
    if numel(dc) < 5, dc(end+1:5) = 0; end
    distCoeffs{c} = dc;
end

% Pack refractive index struct
n = struct('air', cfg.n_air, 'glass', cfg.n_glass, 'water', cfg.n_water);

% -------------------------------------------------------------------------
% Build per-frame, per-camera observation matrix
% obs{f}(c,1:4) = [u1 v1 u2 v2] for wand point A and B in camera c, frame f
% NaN if camera c did not detect wand in frame f
% -------------------------------------------------------------------------
fprintf('Building observation matrix...\n');
obs = cell(nFrames,1);

for f = 1:nFrames
    obs{f} = NaN(nCams, 4);
    for c = 1:nCams
        det = detections{c}(f);
        if det.valid
            % pts is [2×2]: row1 = point A, row2 = point B
            % We use brightness ordering consistently — point A = brighter
            obs{f}(c,:) = [det.pts(1,:), det.pts(2,:)];
        end
    end
end

% Keep only frames seen by >= 2 cameras
frameValid = false(nFrames,1);
for f = 1:nFrames
    nSeen = sum(~any(isnan(obs{f}),2));
    frameValid(f) = nSeen >= 2;
end
validFrames = find(frameValid);

% Subsample frames for faster BA
maxFrames = 400;
if numel(validFrames) > maxFrames
    idx = round(linspace(1, numel(validFrames), maxFrames));
    validFrames = validFrames(idx);
    fprintf('  Subsampled to %d frames for BA\n', numel(validFrames));
end

fprintf('  Frames used in BA: %d / %d\n', numel(validFrames), nFrames);

% -------------------------------------------------------------------------
% INITIALISATION
% -------------------------------------------------------------------------
fprintf('\nInitialising camera poses...\n');
[poses0, wandPts0, validFrames] = initialisePoses(obs, validFrames, Ks, distCoeffs, cfg);

% -------------------------------------------------------------------------
% BUNDLE ADJUSTMENT  — lsqnonlin with analytical structure
% -------------------------------------------------------------------------
fprintf('\nRunning bundle adjustment...\n');
fprintf('  Free parameters:\n');
fprintf('    Camera poses:  %d × 6 = %d params\n', nCams, nCams*6);
fprintf('    Wand poses:    %d × 6 = %d params\n', numel(validFrames), numel(validFrames)*6);

% Parameter vector layout:
%   [cam_params (nCams×6) | wand_params (nValid×6)]
%
%   Camera:  [r1 r2 r3 tx ty tz]  (axis-angle + translation)
%   Wand:    [r1 r2 r3 mx my mz]  (axis-angle of wand orientation + midpoint)

x0 = packParams(poses0, wandPts0, nCams, numel(validFrames));

costHistory = [];
iterCount   = 0;

clear recordCost

opts = optimoptions('lsqnonlin', ...
    'Algorithm',            'levenberg-marquardt', ...
    'MaxIterations',        cfg.ba.maxIter, ...
    'FunctionTolerance',    cfg.ba.fTol, ...
    'StepTolerance',        cfg.ba.xTol, ...
    'Display',              'iter-detailed', ...
    'OutputFcn',            @recordCost, ...
    'SpecifyObjectiveGradient', false);

costFun = @(x) bundleAdjustmentCost(x, obs, validFrames, Ks, distCoeffs, ...
                                     cfg.walls, n, cfg.wandLength, nCams);

[~, initErr] = computePerCameraError(x0, obs, validFrames, Ks, ...
                distCoeffs, cfg.walls, n, cfg.wandLength, nCams);
fprintf('Initial per-camera reprojection errors:\n');
for c = 1:nCams
    fprintf('  Cam %d: %.2f px\n', c, initErr(c));
end

% Bounds: camera translations within ±2m, rotations unbounded
nParams = numel(x0);
nCamParams = nCams * 6;
lb = -inf(nParams, 1);
ub =  inf(nParams, 1);
% Bound camera translations (params 4,5,6 of each camera block)
for c = 1:nCams
    idx = (c-1)*6 + (4:6);
    lb(idx) = -5;
    ub(idx) =  5;
end

[x_opt, ~, residuals, exitflag, output, ~, J_sparse] = ...
    lsqnonlin(costFun, x0, lb, ub, opts);

fprintf('\n  Exit flag: %d  (%s)\n', exitflag, output.message);
fprintf('  Final RMS reprojection error: %.4f px\n', rms(residuals) / sqrt(2));

% -------------------------------------------------------------------------
% UNPACK & COMPUTE UNCERTAINTIES
% -------------------------------------------------------------------------
[poses_opt, wandPts_opt] = unpackParams(x_opt, nCams, numel(validFrames));

% Covariance estimate from Jacobian:  Cov ≈ s² (J'J)^{-1}
s2       = sum(residuals.^2) / (numel(residuals) - numel(x_opt));
JtJ      = full(J_sparse' * J_sparse);
covMatrix = s2 * pinv(JtJ);
paramStd = sqrt(max(0, diag(covMatrix)));

% Per-camera pose uncertainty (first nCams*6 params)
poseStd = reshape(paramStd(1:nCams*6), 6, nCams)';
fprintf('\nCamera pose uncertainties (1σ):\n');
fprintf('%-8s %-18s %-18s\n', 'Camera', 'Rotation (rad)', 'Translation (m)');
for c = 1:nCams
    fprintf('  Cam %d   [%.2e %.2e %.2e]  [%.2e %.2e %.2e]\n', c, ...
        poseStd(c,1), poseStd(c,2), poseStd(c,3), ...
        poseStd(c,4), poseStd(c,5), poseStd(c,6));
end

% Per-camera reprojection error
[~, perCamErr] = computePerCameraError(x_opt, obs, validFrames, Ks, ...
                                        distCoeffs, cfg.walls, n, cfg.wandLength, nCams);
fprintf('\nPer-camera RMS reprojection errors:\n');
for c = 1:nCams
    fprintf('  Camera %d: %.4f px\n', c, perCamErr(c));
end

% -------------------------------------------------------------------------
% OUTLIER DETECTION  (mark frames with high reprojection error)
% -------------------------------------------------------------------------
[~, frameErrors] = computePerFrameError(x_opt, obs, validFrames, Ks, ...
                                         distCoeffs, cfg.walls, n, cfg.wandLength, nCams);
outlierFrames = validFrames(frameErrors > cfg.ba.reprojThresh);
fprintf('\nFrames with reprojection error > %.1f px: %d / %d\n', ...
    cfg.ba.reprojThresh, numel(outlierFrames), numel(validFrames));

% -------------------------------------------------------------------------
% SAVE RESULTS
% -------------------------------------------------------------------------
% Convert axis-angle poses to R, t matrices
cameras = struct();
for c = 1:nCams
    cameras(c).R           = axang2rotm(poses_opt(c,1:3));
    cameras(c).t           = poses_opt(c,4:6)';
    cameras(c).K           = Ks{c};
    cameras(c).distCoeffs  = distCoeffs{c};
    cameras(c).wall        = cfg.walls(c);
    cameras(c).reprojRMS   = perCamErr(c);
    cameras(c).poseStd     = poseStd(c,:);
end

save(fullfile(resultsDir, 'extrinsics.mat'), 'cameras', 'cfg');
save(fullfile(resultsDir, 'bundle_adjustment.mat'), ...
    'cameras', 'wandPts_opt', 'residuals', 'covMatrix', ...
    'frameErrors', 'outlierFrames', 'costHistory', 'cfg');

fprintf('\n[Phase 4] Results saved to %s\n', resultsDir);

% -------------------------------------------------------------------------
% FIGURES
% -------------------------------------------------------------------------
% Convergence plot
fig = figure('Visible','off');
plot(costHistory, 'b-', 'LineWidth', 1.5);
xlabel('Iteration'); ylabel('Total cost (sum sq. residuals)');
title('Bundle Adjustment Convergence'); grid on;
saveas(fig, fullfile(resultsDir, 'ba_convergence.png')); close(fig);

% Per-frame error histogram
fig2 = figure('Visible','off');
histogram(frameErrors, 40, 'FaceColor',[0.2 0.5 0.8]);
xline(cfg.ba.reprojThresh, 'r--', 'LineWidth', 2, 'Label', 'Outlier thresh');
xlabel('RMS reprojection error (px)'); ylabel('Frame count');
title('Per-frame Reprojection Errors After BA'); grid on;
saveas(fig2, fullfile(resultsDir, 'ba_frame_errors.png')); close(fig2);


% =========================================================================
%  OUTPUT FUNCTION — records cost each iteration
% =========================================================================
function stop = recordCost(~, optimValues, state)
    persistent hist
    if strcmp(state, 'init')
        hist = [];
    elseif strcmp(state, 'iter')
        hist(end+1) = optimValues.resnorm;
    elseif strcmp(state, 'done')
        assignin('base', 'costHistory', hist);
    end
    stop = false;
end


% =========================================================================
%  INITIALISATION
% =========================================================================
function [poses0, wandPts0, validFrames] = initialisePoses(obs, validFrames, Ks, distCoeffs, cfg)
% Initialise camera poses from approximate positions specified in cfg.initPoses.
% This is far more robust than essential matrix initialisation for refractive
% systems where the pinhole approximation breaks down at oblique angles.

    nCams  = cfg.nCams;
    wL     = cfg.wandLength;
    n      = struct('air',cfg.n_air,'glass',cfg.n_glass,'water',cfg.n_water);

    fprintf('  Using approximate pose initialisation from LPT_Config\n');

    % -----------------------------------------------------------------------
    % Step 1: Build initial R, t for each camera from C_approx and target
    % -----------------------------------------------------------------------
    R_world = cell(nCams,1);
    t_world = cell(nCams,1);

    for c = 1:nCams
        C   = cfg.initPoses(c).C_approx(:);
        tgt = cfg.initPoses(c).target(:);

        zAxis = tgt - C;
        if norm(zAxis) < 1e-6
            error('Camera %d: C_approx and target are the same point.', c);
        end
        zAxis = zAxis / norm(zAxis);

        up = [0;0;1];
        if abs(dot(zAxis, up)) > 0.9
            up = [0;1;0];
        end
        xAxis = cross(zAxis, up); xAxis = xAxis / norm(xAxis);
        yAxis = cross(zAxis, xAxis);

        R_world{c} = [xAxis, yAxis, zAxis]';
        t_world{c} = -R_world{c} * C;

        fprintf('  Camera %d: C=[%.3f %.3f %.3f] aimed at [%.3f %.3f %.3f]\n', ...
            c, C(1), C(2), C(3), tgt(1), tgt(2), tgt(3));
    end

    % -----------------------------------------------------------------------
    % Step 2: Estimate scale from wand length constraint
    % Uses refraction-aware ray tracing to triangulate wand endpoints
    % and compare to known wand length
    % -----------------------------------------------------------------------
    fprintf('  Estimating scale from wand constraint...\n');

    wallGeoms = cell(nCams,1);
    for c = 1:nCams
        wallGeoms{c} = struct('normal',    cfg.walls(c).normal(:), ...
                               'point',     cfg.walls(c).point(:), ...
                               'thickness', cfg.walls(c).thickness);
    end

    % Find best camera pair for scale estimation (most shared frames)
    pairCount = zeros(nCams);
    for f = validFrames'
        seen = find(~any(isnan(obs{f}(:,1:2)), 2));
        for i = 1:numel(seen)
            for j = i+1:numel(seen)
                pairCount(seen(i),seen(j)) = pairCount(seen(i),seen(j)) + 1;
            end
        end
    end
    [~, bestIdx] = max(pairCount(:));
    [c1, c2]     = ind2sub([nCams nCams], bestIdx);
    fprintf('  Scale estimation using camera pair %d-%d\n', c1, c2);

    % Triangulate wand endpoints using approximate poses (no refraction for speed)
    % and compare reconstructed length to known wL
    scales = [];
    for f = validFrames'
        if any(isnan(obs{f}(c1,:))) || any(isnan(obs{f}(c2,:))), continue; end

        pose1 = struct('R', R_world{c1}, 't', t_world{c1});
        pose2 = struct('R', R_world{c2}, 't', t_world{c2});

        XA = dltTriangulate(obs{f}(c1,1:2)', obs{f}(c2,1:2)', ...
                             Ks{c1}, Ks{c2}, pose1, pose2);
        XB = dltTriangulate(obs{f}(c1,3:4)', obs{f}(c2,3:4)', ...
                             Ks{c1}, Ks{c2}, pose1, pose2);
        d = norm(XA - XB);
        if d > 1e-6 && d < 10  % sanity check
            scales(end+1) = wL / d; 
        end
    end

    if isempty(scales)
        warning('Could not estimate scale — using scale=1. Check approximate poses.');
        scale = 1;
    else
        scale = median(scales);
        fprintf('  Scale factor: %.6f (median of %d estimates)\n', scale, numel(scales));
        fprintf('  Scale min/max: %.6f / %.6f\n', min(scales), max(scales));
        fprintf('  First few scale estimates: '); disp(scales(1:min(5,end)))
    end

    % Recompute translations directly from approximate positions (no scaling)
    for c = 1:nCams
        C = cfg.initPoses(c).C_approx(:);
        t_world{c} = -R_world{c} * C;
    end
    % -----------------------------------------------------------------------
    % Step 3: Pack into axis-angle + translation
    % -----------------------------------------------------------------------
    poses0 = zeros(nCams, 6);
    for c = 1:nCams
        aa = rotm2axang(R_world{c});
        poses0(c,:) = [aa(1:3)*aa(4), t_world{c}'];
    end

    % -----------------------------------------------------------------------
    % Step 4: Initial wand point estimates via refraction-aware triangulation
    % -----------------------------------------------------------------------
    nValid   = numel(validFrames);
    wandPts0 = zeros(nValid, 6);

    camPoses = cell(nCams,1);
    for c = 1:nCams
        camPoses{c} = struct('R', R_world{c}, 't', t_world{c});
    end

    fprintf('  Initialising %d wand poses...\n', nValid);
    goodFrames = true(nValid, 1);
    for fi = 1:nValid
        f    = validFrames(fi);
        seen = find(~any(isnan(obs{f}(:,1:2)), 2));
        
        % Try all camera pairs, use best triangulation
        bestXA = []; bestXB = []; bestDist = inf;
        for pi = 1:numel(seen)
            for pj = pi+1:numel(seen)
                ca = seen(pi); cb = seen(pj);
                % Skip same-side camera pairs (poor geometry)
                if cfg.walls(ca).normal(1) == cfg.walls(cb).normal(1)
                    continue;
                end
                poseA = struct('R',R_world{ca},'t',t_world{ca});
                poseB = struct('R',R_world{cb},'t',t_world{cb});
                XA_try = dltTriangulate(obs{f}(ca,1:2)', obs{f}(cb,1:2)', ...
                                         Ks{ca}, Ks{cb}, poseA, poseB);
                XB_try = dltTriangulate(obs{f}(ca,3:4)', obs{f}(cb,3:4)', ...
                                         Ks{ca}, Ks{cb}, poseA, poseB);
                d = norm(XA_try - XB_try);
                distFromOrigin = max(norm(XA_try), norm(XB_try));
                if distFromOrigin < 0.5 && d < 0.3 && d > 0.05
                    if abs(d - cfg.wandLength) < bestDist
                        bestDist = abs(d - cfg.wandLength);
                        bestXA = XA_try; bestXB = XB_try;
                    end
                end
            end
        end
        
        if isempty(bestXA)
            goodFrames(fi) = false;
            continue;
        end
        
        mid  = (bestXA + bestXB) / 2;
        wdir = bestXB - bestXA;
        wnorm = norm(wdir);
        if wnorm < 1e-8, wdir = [1;0;0]; else, wdir = wdir/wnorm; end
        ax  = cross([1;0;0], wdir);
        axn = norm(ax);
        if axn < 1e-8
            aa_wand = [0 0 0];
        else
            ang = atan2(axn, dot([1;0;0], wdir));
            aa_wand = (ax/axn * ang)';
        end
        wandPts0(fi,:) = [aa_wand, mid'];
    end
    
    % Remove bad frames
    validFrames = validFrames(goodFrames);
    wandPts0    = wandPts0(goodFrames,:);
    fprintf('  Good wand initialisations: %d / %d\n', sum(goodFrames), nValid);

    fprintf('  Initialisation complete.\n');
 end



% =========================================================================
%  PARAMETER PACK / UNPACK
% =========================================================================
function x = packParams(poses, wandPts, nCams, nFrames)
    x = [reshape(poses(1:nCams,:)', [], 1); ...
         reshape(wandPts(1:nFrames,:)', [], 1)];
end

function [poses, wandPts] = unpackParams(x, nCams, nFrames)
    poses   = reshape(x(1:nCams*6), 6, nCams)';
    wandPts = reshape(x(nCams*6+1:end), 6, nFrames)';
end


% =========================================================================
%  BUNDLE ADJUSTMENT COST
% =========================================================================
function r = bundleAdjustmentCost(x, obs, validFrames, Ks, distCoeffs, walls, n, wL, nCams)
% Returns residual vector for lsqnonlin.

    nValid = numel(validFrames);
    [poses, wandPts] = unpackParams(x, nCams, nValid);
    
    % Pre-build pose structs
    camPoses = cell(nCams,1);
    for c = 1:nCams
        aa = poses(c,1:3);
        ang = norm(aa);
        if ang < 1e-10
            R = eye(3);
        else
            R = axang2rotm([aa/ang, ang]);
        end
        camPoses{c} = struct('R', R, 't', poses(c,4:6)');
    end

    % Build wall geometry structs per camera (allow perturbation if enabled)
    wallGeoms = cell(nCams,1);
    for c = 1:nCams
        wallGeoms{c} = struct('normal',    walls(c).normal(:), ...
                               'point',     walls(c).point(:), ...
                               'thickness', walls(c).thickness);
    end

    residuals = [];
    for fi = 1:nValid
        f = validFrames(fi);

        % Decode wand midpoint and direction
        aa_w  = wandPts(fi,1:3);
        mid   = wandPts(fi,4:6)';
        ang   = norm(aa_w);
        if ang < 1e-10
            wdir = [1;0;0];
        else
            wdir = axang2rotm([aa_w/ang, ang]) * [1;0;0];
        end
        wdir = wdir / norm(wdir);

        XA = mid - (wL/2) * wdir;
        XB = mid + (wL/2) * wdir;

        for c = 1:nCams
            if any(isnan(obs{f}(c,:))), continue; end

            uvA_obs = obs{f}(c,1:2)';
            uvB_obs = obs{f}(c,3:4)';

            uvA_proj = projectPointRefractive(XA, camPoses{c}, Ks{c}, ...
                                               distCoeffs{c}, wallGeoms{c}, n);
            uvB_proj = projectPointRefractive(XB, camPoses{c}, Ks{c}, ...
                                               distCoeffs{c}, wallGeoms{c}, n);

            if any(isnan(uvA_proj)) || any(isnan(uvB_proj))
                % Penalty that grows with distance from measurement volume
                % This gives the optimiser a gradient to pull cameras back
                C_c = -camPoses{c}.R' * camPoses{c}.t;
                dist_penalty = 1000 * max(0, norm(C_c) - 1.0);
                residuals = [residuals; 100+dist_penalty; 100+dist_penalty; ...
                                        100+dist_penalty; 100+dist_penalty];
                continue;
            end

            residuals = [residuals; uvA_proj - uvA_obs; uvB_proj - uvB_obs]; %#ok<AGROW>
        end
    end
    r = residuals;
end


% =========================================================================
%  PER-CAMERA AND PER-FRAME ERROR
% =========================================================================
function [rms_total, perCamRMS] = computePerCameraError(x, obs, validFrames, ...
                                    Ks, distCoeffs, walls, n, wL, nCams)
    nValid = numel(validFrames);
    [poses, wandPts] = unpackParams(x, nCams, nValid);
    camPoses  = buildCamPoses(poses, nCams);
    wallGeoms = buildWallGeoms(walls, nCams);

    errSums = zeros(nCams,1);
    counts  = zeros(nCams,1);

    for fi = 1:nValid
        f = validFrames(fi);
        [XA, XB] = decodeWand(wandPts(fi,:), wL);
        for c = 1:nCams
            if any(isnan(obs{f}(c,:))), continue; end
            uvA = projectPointRefractive(XA, camPoses{c}, Ks{c}, distCoeffs{c}, wallGeoms{c}, n);
            uvB = projectPointRefractive(XB, camPoses{c}, Ks{c}, distCoeffs{c}, wallGeoms{c}, n);
            if any(isnan(uvA)) || any(isnan(uvB)), continue; end
            eA = norm(uvA - obs{f}(c,1:2)');
            eB = norm(uvB - obs{f}(c,3:4)');
            errSums(c) = errSums(c) + eA^2 + eB^2;
            counts(c)  = counts(c)  + 2;
        end
    end
    perCamRMS = sqrt(errSums ./ max(counts,1));
    rms_total = sqrt(sum(errSums) / max(sum(counts),1));
end

function [rms_total, frameRMS] = computePerFrameError(x, obs, validFrames, ...
                                    Ks, distCoeffs, walls, n, wL, nCams)
    nValid = numel(validFrames);
    [poses, wandPts] = unpackParams(x, nCams, nValid);
    camPoses  = buildCamPoses(poses, nCams);
    wallGeoms = buildWallGeoms(walls, nCams);

    frameRMS = zeros(nValid,1);
    for fi = 1:nValid
        f = validFrames(fi);
        [XA, XB] = decodeWand(wandPts(fi,:), wL);
        errSum = 0; cnt = 0;
        for c = 1:nCams
            if any(isnan(obs{f}(c,:))), continue; end
            uvA = projectPointRefractive(XA, camPoses{c}, Ks{c}, distCoeffs{c}, wallGeoms{c}, n);
            uvB = projectPointRefractive(XB, camPoses{c}, Ks{c}, distCoeffs{c}, wallGeoms{c}, n);
            if any(isnan(uvA)) || any(isnan(uvB)), continue; end
            errSum = errSum + norm(uvA-obs{f}(c,1:2)')^2 + norm(uvB-obs{f}(c,3:4)')^2;
            cnt = cnt + 2;
        end
        frameRMS(fi) = sqrt(errSum / max(cnt,1));
    end
    rms_total = rms(frameRMS);
end


% =========================================================================
%  UTILITY FUNCTIONS
% =========================================================================
function camPoses = buildCamPoses(poses, nCams)
    camPoses = cell(nCams,1);
    for c = 1:nCams
        aa = poses(c,1:3); ang = norm(aa);
        if ang < 1e-10
            R = eye(3);
        else
            R = axang2rotm([aa/ang, ang]);
        end
        camPoses{c} = struct('R', R, 't', poses(c,4:6)');
    end
end

function wallGeoms = buildWallGeoms(walls, nCams)
    wallGeoms = cell(nCams,1);
    for c = 1:nCams
        wallGeoms{c} = struct('normal',    walls(c).normal(:), ...
                               'point',     walls(c).point(:), ...
                               'thickness', walls(c).thickness);
    end
end

function [XA, XB] = decodeWand(wp, wL)
    aa = wp(1:3); ang = norm(aa);
    if ang < 1e-10, wdir = [1;0;0];
    else, wdir = axang2rotm([aa/ang, ang]) * [1;0;0]; end
    wdir = wdir / norm(wdir);
    mid  = wp(4:6)';
    XA   = mid - (wL/2)*wdir;
    XB   = mid + (wL/2)*wdir;
end

function R = axang2rotm(aa)
% Axis-angle to rotation matrix. aa = [ax ay az angle] or [ax ay az] (angle=norm).
    if numel(aa)==3
        angle = norm(aa);
        if angle < 1e-10, R = eye(3); return; end
        ax = aa(:)/angle;
    else
        ax = aa(1:3)'; angle = aa(4);
    end
    c=cos(angle); s=sin(angle); t=1-c;
    x=ax(1); y=ax(2); z=ax(3);
    R = [t*x*x+c,   t*x*y-s*z, t*x*z+s*y;
         t*x*y+s*z, t*y*y+c,   t*y*z-s*x;
         t*x*z-s*y, t*y*z+s*x, t*z*z+c  ];
end

function aa = rotm2axang(R)
    angle = acos(max(-1, min(1, (trace(R)-1)/2)));
    if abs(angle) < 1e-10
        aa = [1 0 0 0]; return;
    end
    ax = [(R(3,2)-R(2,3)); (R(1,3)-R(3,1)); (R(2,1)-R(1,2))] / (2*sin(angle));
    aa = [ax(:)', angle];
end

function X = dltTriangulate(uv1, uv2, K1, K2, pose1, pose2)
% Simple DLT triangulation between two cameras (no refraction).
    P1 = K1 * [pose1.R, pose1.t];
    P2 = K2 * [pose2.R, pose2.t];
    A  = [uv1(1)*P1(3,:)-P1(1,:);
          uv1(2)*P1(3,:)-P1(2,:);
          uv2(1)*P2(3,:)-P2(1,:);
          uv2(2)*P2(3,:)-P2(2,:)];
    [~,~,V] = svd(A);
    X = V(1:3,end) / V(4,end);
end

function [R, t] = solvePnP(worldPts, imgPts, K, distCoeffs, imageSize)
% Direct linear PnP without RANSAC — more reliable for clean synthetic data

    n  = size(worldPts, 1);
    A  = zeros(2*n, 12);
    for i = 1:n
        X = worldPts(i,:);
        u = imgPts(i,1);
        v = imgPts(i,2);
        A(2*i-1,:) = [X, 1, 0, 0, 0, 0, -u*X, -u];
        A(2*i,  :) = [0, 0, 0, 0, X, 1, -v*X, -v];
    end

    [~,~,V] = svd(A);
    P = reshape(V(:,end), 4, 3)';

    % Extract R and t from P = K * [R | t]
    M  = K \ P(:,1:3);
    t  = K \ P(:,4);

    % Enforce rotation matrix via SVD
    [U,~,V] = svd(M);
    R = U * V';
    if det(R) < 0, R = -R; t = -t; end

    % Resolve scale
    scale = norm(M) / norm(R);
    t     = t / scale;
end

function [E, inliers] = computeEssentialMatrix(pts1, pts2, K1, K2, imageSize)
% Estimate essential matrix using MATLAB's built-in.
    intr1 = cameraIntrinsics([K1(1,1), K1(2,2)], K1(1:2,3)', fliplr(imageSize));
    intr2 = cameraIntrinsics([K2(1,1), K2(2,2)], K2(1:2,3)', fliplr(imageSize));
    [E, inliers] = estimateEssentialMatrix(pts1, pts2, intr1, intr2, ...
        'Confidence', 99.99, 'MaxDistance', 1.0);
end

function [R, t] = recoverPose(E, pts1, pts2, K)
% Recover R, t from essential matrix.
    [U,~,V] = svd(E);
    W = [0 -1 0; 1 0 0; 0 0 1];
    R1 = U*W*V';   R2 = U*W'*V';
    t1 = U(:,3);
    if det(R1) < 0, R1 = -R1; end
    if det(R2) < 0, R2 = -R2; end
    % Pick solution where points are in front of both cameras
    [R, t] = pickCorrectPose(R1, R2, t1, -t1, pts1, pts2, K);
end

function [R, t] = pickCorrectPose(R1, R2, t1, t2, pts1, pts2, K)
    candidates = {R1,t1; R1,t2; R2,t1; R2,t2};
    bestScore  = -1;
    R = R1; t = t1;
    P1 = K * eye(3,4);
    for i = 1:4
        P2   = K * [candidates{i,1}, candidates{i,2}];
        X    = triangulate(pts1(1,:), pts2(1,:), P1, P2);
        Xc2  = candidates{i,1} * X(1:3)' + candidates{i,2};
        if X(3) > 0 && Xc2(3) > 0
            score = sum(X(:,3) > 0);
            if score > bestScore
                bestScore = score; R = candidates{i,1}; t = candidates{i,2};
            end
        end
    end
end

function scale = estimateScaleFromWand(obs, validFrames, c1, c2, K1, K2, ...
                                        dc1, dc2, R_rel, t_rel, wL)
% Estimate scale by triangulating wand endpoints and comparing to known wL.
    scales = [];
    for f = validFrames'
        if any(isnan(obs{f}(c1,:))) || any(isnan(obs{f}(c2,:))), continue; end
        pose1 = struct('R',eye(3),'t',zeros(3,1));
        pose2 = struct('R',R_rel,'t',t_rel);
        XA = dltTriangulate(obs{f}(c1,1:2)', obs{f}(c2,1:2)', K1, K2, pose1, pose2);
        XB = dltTriangulate(obs{f}(c1,3:4)', obs{f}(c2,3:4)', K1, K2, pose1, pose2);
        d  = norm(XA - XB);
        if d > 1e-6
            scales(end+1) = wL / d; %#ok<AGROW>
        end
    end
    scale = median(scales);
end
