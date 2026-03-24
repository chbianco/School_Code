% =========================================================================
% EvaluateCalibration.m
% =========================================================================
% PURPOSE:
%   Quantitatively evaluate calibration accuracy after bundle adjustment.
%
%   In SYNTHETIC MODE (ground_truth.mat exists):
%     • Compare recovered R, t to ground truth for each camera
%     • Report rotation and translation errors
%     • Report wand length accuracy
%     • Report 3-D triangulation accuracy
%
%   In REAL DATA MODE (no ground truth):
%     • Report per-camera reprojection errors
%     • Check wand length consistency across frames
%     • Check epipolar geometry consistency
%     • Report 3-D uncertainty estimates
% =========================================================================

clear; clc;
LPT_Config;

fprintf('=== Calibration Evaluation ===\n\n');

resultsDir = cfg.resultsDir;
nCams      = cfg.nCams;

cal = load(fullfile(resultsDir, 'extrinsics.mat'));
ba  = load(fullfile(resultsDir, 'bundle_adjustment.mat'));
cameras = cal.cameras;

n   = struct('air',cfg.n_air,'glass',cfg.n_glass,'water',cfg.n_water);

hasSynth = exist(fullfile(resultsDir,'ground_truth.mat'),'file');

% =========================================================================
%  SECTION A: Reprojection errors (always available)
% =========================================================================
fprintf('=== A. Reprojection Error Summary ===\n\n');
fprintf('%-10s %-20s %-20s\n', 'Camera', 'RMS Reproj (px)', 'Max Reproj (px)');
fprintf('%s\n', repmat('-',1,52));
% Load observations and recompute reprojection errors fresh
% Use the BA's own validFrames and obs if available (saved by Phase4)
if isfield(ba, 'validFrames') && isfield(ba, 'obs')
    validFrames_tmp = ba.validFrames;
    obs_tmp         = ba.obs;
    nFrames_tmp     = numel(ba.obs);
    fprintf('  Using BA-saved validFrames (%d frames)\n\n', numel(validFrames_tmp));
else
    % Fallback: reconstruct from Phase 2 detections
    phase2_tmp = load(fullfile(resultsDir,'wand_detections.mat'));
    det_tmp    = phase2_tmp.detections;
    nFrames_tmp = phase2_tmp.nFrames;

    obs_tmp = cell(nFrames_tmp,1);
    for f = 1:nFrames_tmp
        obs_tmp{f} = NaN(nCams,4);
        for c = 1:nCams
            if det_tmp{c}(f).valid
                obs_tmp{f}(c,:) = [det_tmp{c}(f).pts(1,:), det_tmp{c}(f).pts(2,:)];
            end
        end
    end
    % Match to the number of wand points saved
    nWand = size(ba.wandPts_opt, 1);
    allValid = find(cellfun(@(o) sum(~any(isnan(o),2))>=2, obs_tmp));
    if numel(allValid) > nWand
        % Subsample to match BA frame count
        validFrames_tmp = allValid(round(linspace(1, numel(allValid), nWand)));
    else
        validFrames_tmp = allValid(1:min(numel(allValid), nWand));
    end
    fprintf('  Reconstructed validFrames (%d frames, %d wand pts)\n\n', ...
        numel(validFrames_tmp), nWand);
end

% Recompute per-camera errors
Ks_tmp        = {cameras.K};
dc_tmp        = {cameras.distCoeffs};
wallGeoms_tmp = arrayfun(@(c) struct('normal',   cameras(c).wall.normal(:), ...
                                      'point',    cameras(c).wall.point(:), ...
                                      'thickness',cameras(c).wall.thickness), ...
                          1:nCams, 'UniformOutput', false);

errSums = zeros(nCams,1);
counts  = zeros(nCams,1);
for fi = 1:numel(validFrames_tmp)
    f = validFrames_tmp(fi);
    wp = ba.wandPts_opt(fi,:);
    aa = wp(1:3); ang = norm(aa);
    if ang < 1e-10, wdir=[1;0;0];
    else, wdir = axang2rotm_local([aa/ang,ang])*[1;0;0]; end
    wdir = wdir/norm(wdir);
    mid  = wp(4:6)';
    XA   = mid - (cfg.wandLength/2)*wdir;
    XB   = mid + (cfg.wandLength/2)*wdir;

    for c = 1:nCams
        if any(isnan(obs_tmp{f}(c,:))), continue; end
        pose_c = struct('R', cameras(c).R, 't', cameras(c).t);
        uvA = projectPointRefractive(XA, pose_c, Ks_tmp{c}, dc_tmp{c}, wallGeoms_tmp{c}, n);
        uvB = projectPointRefractive(XB, pose_c, Ks_tmp{c}, dc_tmp{c}, wallGeoms_tmp{c}, n);
        if any(isnan(uvA)) || any(isnan(uvB)), continue; end
        errSums(c) = errSums(c) + norm(uvA-obs_tmp{f}(c,1:2)')^2 ...
                                + norm(uvB-obs_tmp{f}(c,3:4)')^2;
        counts(c)  = counts(c) + 2;
    end
end

for c = 1:nCams
    perCamRMS = sqrt(errSums(c) / max(counts(c),1));
    fprintf('  Cam %d     %.4f               %.4f\n', c, perCamRMS, nan);
end

% Overall RMS
nObs     = numel(ba.residuals);
overallRMS = rms(ba.residuals) / sqrt(2);
fprintf('\nOverall RMS (all cameras, all frames): %.4f px\n', overallRMS);

if overallRMS < 0.5
    fprintf('[✓] EXCELLENT — sub-0.5 px RMS is very good for LPT\n');
elseif overallRMS < 1.0
    fprintf('[~] ACCEPTABLE — consider re-running wand collection\n');
else
    fprintf('[✗] POOR — reprojection errors >1 px suggest problems\n');
end

% =========================================================================
%  SECTION B: Wand length consistency
% =========================================================================
fprintf('\n=== B. Wand Length Consistency ===\n\n');

wandPts  = ba.wandPts_opt;
wL_known = cfg.wandLength;
wL_meas  = zeros(size(wandPts,1),1);

for fi = 1:size(wandPts,1)
    wp = wandPts(fi,:);
    aa = wp(1:3); ang = norm(aa);
    if ang < 1e-10, wdir=[1;0;0];
    else, wdir = axang2rotm_local([aa/ang,ang]) * [1;0;0]; end
    wdir = wdir/norm(wdir);
    mid  = wp(4:6)';
    XA   = mid - (wL_known/2)*wdir;
    XB   = mid + (wL_known/2)*wdir;

    % Re-triangulate with current calibration
    obs_f = NaN(nCams,4);
    % (We use the bundle-adjusted 3-D positions directly here)
    wL_meas(fi) = norm(XA - XB);  % will equal wL by construction; more useful below
end

% More meaningful: triangulate from observations and check length
phase2 = load(fullfile(resultsDir,'wand_detections.mat'));
detections = phase2.detections;
frameValid = find(arrayfun(@(f) any(cellfun(@(d) d(f).valid, detections)), 1:phase2.nFrames));

Ks         = {cameras.K};
distCoeffs = {cameras.distCoeffs};
camPoses   = arrayfun(@(c) struct('R',cameras(c).R,'t',cameras(c).t), 1:nCams,'UniformOutput',false);
wallGeoms  = arrayfun(@(c) struct('normal',cameras(c).wall.normal(:),...
                                   'point',cameras(c).wall.point(:),...
                                   'thickness',cameras(c).wall.thickness), ...
                       1:nCams,'UniformOutput',false);

reconLengths = [];
for f = frameValid(1:min(200,end))
    obs_A = NaN(nCams,2); obs_B = NaN(nCams,2);
    for c = 1:nCams
        det = detections{c}(f);
        if det.valid
            obs_A(c,:) = det.pts(1,:);
            obs_B(c,:) = det.pts(2,:);
        end
    end
    rA = triangulateRefractive(obs_A, camPoses, Ks, distCoeffs, wallGeoms, n, cfg.recon);
    rB = triangulateRefractive(obs_B, camPoses, Ks, distCoeffs, wallGeoms, n, cfg.recon);
    if ~any(isnan(rA.pos)) && ~any(isnan(rB.pos))
        reconLengths(end+1) = norm(rA.pos - rB.pos) * 1000;  % mm
    end
end

if ~isempty(reconLengths)
    fprintf('Known wand length:           %.2f mm\n', wL_known*1000);
    fprintf('Reconstructed mean ± std:    %.2f ± %.2f mm\n', mean(reconLengths), std(reconLengths));
    fprintf('Relative error:              %.3f%%\n', ...
        100*abs(mean(reconLengths)-wL_known*1000)/(wL_known*1000));
    fprintf('Coefficient of variation:    %.3f%%\n', 100*std(reconLengths)/mean(reconLengths));
end

% =========================================================================
%  SECTION C: Ground truth comparison (synthetic mode)
% =========================================================================
if hasSynth
    fprintf('\n=== C. Ground Truth Comparison (Synthetic Mode) ===\n\n');
    gt = load(fullfile(resultsDir,'ground_truth.mat'));

    fprintf('%-10s %-20s %-20s\n','Camera','Rot error (deg)','Trans error (mm)');
    fprintf('%s\n', repmat('-',1,52));

    for c = 1:nCams
        R_gt  = gt.gt_cameras(c).R;
        t_gt  = gt.gt_cameras(c).t;
        R_est = cameras(c).R;
        t_est = cameras(c).t;

        % Rotation error: angle of R_est * R_gt^{-1}
        dR    = R_est * R_gt';
        angle = acos(max(-1, min(1, (trace(dR)-1)/2)));
        rot_err_deg = rad2deg(angle);

        % Translation error (in world coords): compare camera centres
        C_gt  = -R_gt'  * t_gt;
        C_est = -R_est' * t_est;
        t_err_mm = norm(C_gt - C_est) * 1000;

        fprintf('  Cam %d     %.4f deg            %.4f mm\n', c, rot_err_deg, t_err_mm);
    end

    % 3-D reconstruction accuracy on synthetic wand
    reconErrors = [];
    for f = 1:min(100, gt.nFrames)
        mid_gt   = gt.gt_wand{f}.midpoint;
        obs_A = NaN(nCams,2); obs_B = NaN(nCams,2);
        for c = 1:nCams
            det = detections{c}(f);
            if det.valid
                obs_A(c,:) = det.pts(1,:);
                obs_B(c,:) = det.pts(2,:);
            end
        end
        rA = triangulateRefractive(obs_A, camPoses, Ks, distCoeffs, wallGeoms, n, cfg.recon);
        rB = triangulateRefractive(obs_B, camPoses, Ks, distCoeffs, wallGeoms, n, cfg.recon);
        if ~any(isnan(rA.pos)) && ~any(isnan(rB.pos))
            mid_est = (rA.pos + rB.pos) / 2;
            reconErrors(end+1) = norm(mid_est - mid_gt) * 1000;  % mm
        end
    end

    if ~isempty(reconErrors)
        fprintf('\n3-D Reconstruction Accuracy (midpoint vs ground truth):\n');
        fprintf('  Mean error: %.4f mm\n', mean(reconErrors));
        fprintf('  Std:        %.4f mm\n', std(reconErrors));
        fprintf('  Max error:  %.4f mm\n', max(reconErrors));
        fprintf('  Added image noise was: %.2f px\n', gt.imageNoise_sigma);
        fprintf('  Expected 3-D error at typical distance:\n');
    end
end

% =========================================================================
%  SECTION D: Epipolar consistency check
% =========================================================================
fprintf('\n=== D. Epipolar Consistency ===\n\n');
fprintf('Checking fundamental matrix consistency for each camera pair...\n\n');

for c1 = 1:nCams
    for c2 = c1+1:nCams
        R1=cameras(c1).R; t1=cameras(c1).t; K1=cameras(c1).K;
        R2=cameras(c2).R; t2=cameras(c2).t; K2=cameras(c2).K;

        % Relative pose
        R12 = R2 * R1';
        t12 = t2 - R2*R1'*t1;

        % Fundamental matrix
        tx  = [0 -t12(3) t12(2); t12(3) 0 -t12(1); -t12(2) t12(1) 0];
        E   = tx * R12;
        F   = inv(K2)' * E * inv(K1);

        % Check epipolar residuals on shared wand observations
        errs = [];
        for f = frameValid(1:min(100,end))
            det1 = detections{c1}(f); det2 = detections{c2}(f);
            if ~det1.valid || ~det2.valid, continue; end
            for pt = 1:2
                p1  = [det1.pts(pt,:), 1]';
                p2  = [det2.pts(pt,:), 1]';
                l2  = F * p1;        % epipolar line in cam2
                d   = abs(p2'*l2) / norm(l2(1:2));   % point-to-line distance
                errs(end+1) = d; 
            end
        end

        if ~isempty(errs)
            fprintf('  Cam %d – Cam %d:  epipolar dist = %.4f ± %.4f px  (%d pts)\n', ...
                c1, c2, mean(errs), std(errs), numel(errs));
        end
    end
end

fprintf('\n[EvaluateCalibration] Done.\n');

% Helper
function R = axang2rotm_local(aa)
    if numel(aa)==3
        angle=norm(aa); if angle<1e-10,R=eye(3);return; end; ax=aa(:)/angle;
    else, ax=aa(1:3)'; angle=aa(4); end
    c=cos(angle);s=sin(angle);t=1-c;x=ax(1);y=ax(2);z=ax(3);
    R=[t*x*x+c,t*x*y-s*z,t*x*z+s*y;t*x*y+s*z,t*y*y+c,t*y*z-s*x;t*x*z-s*y,t*y*z+s*x,t*z*z+c];
end