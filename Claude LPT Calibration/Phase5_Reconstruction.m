% =========================================================================
% Phase5_Reconstruction.m
% =========================================================================
% PURPOSE:
%   Given calibrated cameras (Phase 4), reconstruct 3-D particle positions
%   from multi-camera 2-D observations using refraction-corrected
%   triangulation (Phase 3 model).
%
%   This script handles BOTH:
%     (a) Wand/calibration target reconstruction (for validation)
%     (b) General particle tracking reconstruction
%
%   TRIANGULATION METHOD:
%     For each particle observed in ≥2 cameras:
%     1. Trace a refracted ray from each camera through its glass wall
%     2. Find the 3-D point minimising the sum of squared angular distances
%        to all refracted rays (nonlinear least squares)
%     3. Compute reprojection errors and flag outliers
%     4. Propagate uncertainty via the linearised refractive Jacobian
%
%   INPUT FORMAT:
%     Particle observations should be provided as a cell array or struct:
%       particles{f}(p,c) = [u, v]  image coords of particle p in camera c
%       particles{f}(p,c) = [NaN NaN] if particle p not seen by camera c
%
% OUTPUTS (saved to results/):
%   reconstruction.mat    — 3-D positions + uncertainties per particle/frame
%   reconstruction_QC.png — reprojection error summary figure
% =========================================================================

clear; clc;
LPT_Config;

fprintf('=== Phase 5: 3-D Reconstruction ===\n\n');

resultsDir = cfg.resultsDir;

% -------------------------------------------------------------------------
% Load calibration
% -------------------------------------------------------------------------
cal = load(fullfile(resultsDir, 'extrinsics.mat'));
cameras = cal.cameras;
nCams   = cfg.nCams;

% Repack into cell arrays for the refractive functions
Ks          = cell(nCams,1);
distCoeffs  = cell(nCams,1);
camPoses    = cell(nCams,1);
wallGeoms   = cell(nCams,1);

for c = 1:nCams
    Ks{c}         = cameras(c).K;
    distCoeffs{c} = cameras(c).distCoeffs;
    camPoses{c}   = struct('R', cameras(c).R, 't', cameras(c).t);
    wallGeoms{c}  = struct('normal',    cameras(c).wall.normal(:), ...
                            'point',     cameras(c).wall.point(:), ...
                            'thickness', cameras(c).wall.thickness);
end

n = struct('air', cfg.n_air, 'glass', cfg.n_glass, 'water', cfg.n_water);

% -------------------------------------------------------------------------
% Load particle observations
% -------------------------------------------------------------------------
% Expected file: results/particle_observations.mat
% Variable: particles  — cell array {nFrames} of [nParticles × nCams × 2] arrays
%                        or struct array with fields .uv{cam}

particleFile = fullfile(resultsDir, 'particle_observations.mat');
if ~exist(particleFile, 'file')
    fprintf('[INFO] No particle_observations.mat found.\n');
    fprintf('       Running in VALIDATION MODE using wand detections.\n\n');
    runValidationMode = true;
else
    runValidationMode = false;
    pData = load(particleFile);
    particles = pData.particles;
end

if runValidationMode
    % Use Phase-2 wand detections as test particles to validate reconstruction
    phase2 = load(fullfile(resultsDir, 'wand_detections.mat'));
    detections = phase2.detections;
    nFrames    = phase2.nFrames;
    phase4     = load(fullfile(resultsDir, 'bundle_adjustment.mat'));

    % Convert wand endpoint observations to particle format
    particles  = cell(nFrames,1);
    for f = 1:nFrames
        obs = NaN(2, nCams, 2);   % [nParticles=2, nCams, uv=2]
        for c = 1:nCams
            det = detections{c}(f);
            if det.valid
                obs(1,c,:) = det.pts(1,:);
                obs(2,c,:) = det.pts(2,:);
            end
        end
        particles{f} = obs;
    end
end

nFrames = numel(particles);
fprintf('Reconstructing %d frames...\n\n', nFrames);

% -------------------------------------------------------------------------
% RECONSTRUCT FRAME BY FRAME
% -------------------------------------------------------------------------
recon(nFrames) = struct('frame',0,'positions',[],'reprojErrors',[],...
                        'uncertainty',[],'nCamsSeen',[]);

totalPts   = 0;
goodPts    = 0;
allErrors  = [];

for f = 1:nFrames
    obs3d = particles{f};  % [nPart × nCams × 2]
    if isempty(obs3d)
        recon(f).frame = f;
        continue;
    end
    nPart = size(obs3d, 1);

    positions   = NaN(nPart, 3);
    reproj      = NaN(nPart, nCams);
    uncertainty = NaN(nPart, 3);  % 1-sigma position uncertainty [m]
    nSeen       = zeros(nPart, 1);

    for p = 1:nPart
        % Build observation matrix for this particle
        uvObs = squeeze(obs3d(p,:,:));   % [nCams × 2]

        % Count visible cameras
        visible = ~any(isnan(uvObs), 2);
        nSeen(p) = sum(visible);

        if nSeen(p) < cfg.recon.minCams
            continue;
        end

        % Triangulate
        result = triangulateRefractive(uvObs, camPoses, Ks, distCoeffs, ...
                                        wallGeoms, n, cfg.recon);

        if any(isnan(result.pos))
            continue;
        end

        % Check reprojection errors
        maxErr = max(result.reprojErrors(~isnan(result.reprojErrors)));
        if maxErr > cfg.recon.maxReprojErr
            continue;
        end

        % Compute position uncertainty via Jacobian propagation
        sigma_pos = computePositionUncertainty(result.pos, camPoses, Ks, ...
                                                distCoeffs, wallGeoms, n, ...
                                                result.reprojErrors, visible);

        positions(p,:)   = result.pos';
        reproj(p,:)      = result.reprojErrors';
        uncertainty(p,:) = sigma_pos';
        totalPts         = totalPts + 1;
        goodPts          = goodPts  + 1;
        allErrors        = [allErrors; result.reprojErrors(~isnan(result.reprojErrors))]; %#ok<AGROW>
    end

    recon(f).frame       = f;
    recon(f).positions   = positions;
    recon(f).reprojErrors = reproj;
    recon(f).uncertainty  = uncertainty;
    recon(f).nCamsSeen   = nSeen;
    totalPts = totalPts + nPart;

    if mod(f,100)==0
        fprintf('  Frame %d / %d  |  particles reconstructed: %d\n', f, nFrames, goodPts);
    end
end

fprintf('\nReconstruction complete.\n');
fprintf('  Total valid reconstructions: %d\n', goodPts);
if ~isempty(allErrors)
    fprintf('  Mean reprojection error:     %.4f px\n', mean(allErrors));
    fprintf('  Median reprojection error:   %.4f px\n', median(allErrors));
    fprintf('  Max reprojection error:      %.4f px\n', max(allErrors));
end

% -------------------------------------------------------------------------
% VALIDATION: Compare reconstructed wand length to known value
% -------------------------------------------------------------------------
if runValidationMode
    phase4 = load(fullfile(resultsDir, 'bundle_adjustment.mat'));
    measuredLengths = [];
    for f = 1:nFrames
        if size(recon(f).positions,1) >= 2 && ...
           ~any(isnan(recon(f).positions(1,:))) && ~any(isnan(recon(f).positions(2,:)))
            d = norm(recon(f).positions(1,:) - recon(f).positions(2,:));
            measuredLengths(end+1) = d; %#ok<AGROW>
        end
    end
    if ~isempty(measuredLengths)
        fprintf('\nWand Length Validation:\n');
        fprintf('  Known wand length:       %.4f m\n', cfg.wandLength);
        fprintf('  Measured (mean ± std):   %.4f ± %.4f m\n', ...
            mean(measuredLengths), std(measuredLengths));
        fprintf('  Relative error:          %.3f%%\n', ...
            100*abs(mean(measuredLengths)-cfg.wandLength)/cfg.wandLength);
    end
end

% -------------------------------------------------------------------------
% SAVE
% -------------------------------------------------------------------------
save(fullfile(resultsDir, 'reconstruction.mat'), 'recon', 'cfg');
fprintf('\n[Phase 5] Results saved to %s\n', resultsDir);

% -------------------------------------------------------------------------
% QC FIGURES
% -------------------------------------------------------------------------
if ~isempty(allErrors)
    fig = figure('Visible','off', 'Position',[0 0 1200 400]);

    subplot(1,3,1);
    histogram(allErrors, 50, 'FaceColor',[0.2 0.5 0.8]);
    xline(cfg.recon.maxReprojErr,'r--','LineWidth',2,'Label','Threshold');
    xlabel('Reprojection error (px)'); ylabel('Count');
    title('All reprojection errors'); grid on;

    subplot(1,3,2);
    perCamMean = nanmean(vertcat(recon.reprojErrors), 1);
    bar(perCamMean, 'FaceColor',[0.3 0.7 0.4]);
    xlabel('Camera'); ylabel('Mean reprojection error (px)');
    title('Per-camera mean error'); grid on;

    subplot(1,3,3);
    if runValidationMode && ~isempty(measuredLengths)
        histogram(measuredLengths*1000, 30, 'FaceColor',[0.8 0.4 0.2]);
        xline(cfg.wandLength*1000,'r--','LineWidth',2,'Label','True length');
        xlabel('Reconstructed wand length (mm)'); ylabel('Count');
        title('Wand length distribution (validation)'); grid on;
    end

    sgtitle('Reconstruction Quality Control');
    saveas(fig, fullfile(resultsDir, 'reconstruction_QC.png'));
    close(fig);
    fprintf('[Phase 5] QC figure saved.\n');
end


% =========================================================================
%  POSITION UNCERTAINTY  (linearised error propagation)
% =========================================================================
function sigma_pos = computePositionUncertainty(X, camPoses, Ks, distCoeffs, ...
                                                  wallGeoms, n, reprojErrors, visible)
% Propagate reprojection errors to 3-D position uncertainty.
%
% Uses the pseudo-inverse of the stacked Jacobian:
%   [J1; J2; ...] * dX = [duv1; duv2; ...]
%   Cov(X) = (J'J)^{-1} J' diag(sigma_obs^2) J (J'J)^{-1}

    nCams = numel(camPoses);
    J_stack = [];
    W_stack = [];

    for c = 1:nCams
        if ~visible(c), continue; end
        J_c = refractionJacobian(X, camPoses{c}, Ks{c}, distCoeffs{c}, wallGeoms{c}, n);
        sig_c = max(reprojErrors(c), 0.1);   % floor at 0.1 px
        J_stack = [J_stack; J_c];             %#ok<AGROW>
        W_stack = [W_stack; sig_c; sig_c];    %#ok<AGROW>
    end

    if isempty(J_stack) || rank(J_stack) < 3
        sigma_pos = [Inf; Inf; Inf];
        return;
    end

    W    = diag(W_stack.^2);
    JtWJ = J_stack' * (W \ J_stack);

    try
        Cov_X    = inv(JtWJ);
        sigma_pos = sqrt(max(0, diag(Cov_X)));
    catch
        sigma_pos = [Inf; Inf; Inf];
    end
end
