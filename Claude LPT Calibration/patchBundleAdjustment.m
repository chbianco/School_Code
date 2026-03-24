% PatchBundleAdjustment.m
% Adds validFrames and obs to bundle_adjustment.mat without re-running Phase4.
% Replays the same data preparation (LED swap, detection filter, frame selection).

clear; clc;
LPT_Config;

resultsDir = cfg.resultsDir;
nCams = cfg.nCams;

% Load Phase 2 detections
phase2 = load(fullfile(resultsDir, 'wand_detections.mat'));
detections = phase2.detections;
nFrames = phase2.nFrames;

% Load BA to get the number of wand points
ba = load(fullfile(resultsDir, 'bundle_adjustment.mat'));
nWand = size(ba.wandPts_opt, 1);
fprintf('BA has %d wand points\n', nWand);

% Build observations (same as Phase4)
obs = cell(nFrames,1);
for f = 1:nFrames
    obs{f} = NaN(nCams,4);
    for c = 1:nCams
        det = detections{c}(f);
        if det.valid, obs{f}(c,:) = [det.pts(1,:), det.pts(2,:)]; end
    end
end

% LED swap (must match Phase4 setting)
swapCams = [1, 4];
fprintf('Applying LED swap for cameras [%s]\n', num2str(swapCams));
for f = 1:nFrames
    for c = swapCams
        if ~any(isnan(obs{f}(c,:)))
            obs{f}(c,:) = [obs{f}(c,3:4), obs{f}(c,1:2)];
        end
    end
end

% Camera exclusion (must match Phase4 setting)
excludeCams = [];
if ~isempty(excludeCams)
    for f = 1:nFrames
        for c = excludeCams
            obs{f}(c,:) = NaN(1,4);
        end
    end
end

% Detection quality filter (same as Phase4)
minWandPx = 100; maxWandPx = 800;
for f = 1:nFrames
    for c = 1:nCams
        if ~any(isnan(obs{f}(c,:)))
            wandPxLen = sqrt((obs{f}(c,1)-obs{f}(c,3))^2 + (obs{f}(c,2)-obs{f}(c,4))^2);
            if wandPxLen < minWandPx || wandPxLen > maxWandPx
                obs{f}(c,:) = NaN(1,4);
            end
        end
    end
end

% Find valid frames (≥2 cameras)
allValid = find(arrayfun(@(f) sum(~any(isnan(obs{f}),2))>=2, (1:nFrames)'));
fprintf('Total valid frames: %d\n', numel(allValid));

% Subsample to 400 max (same as Phase4)
maxFrames = 400;
if numel(allValid) > maxFrames
    validFrames = allValid(round(linspace(1, numel(allValid), maxFrames)));
else
    validFrames = allValid;
end
fprintf('After subsampling: %d frames\n', numel(validFrames));

% Now we need to match this to the BA's nWand frames.
% The BA may have dropped additional frames during triangulation filtering.
% We can't perfectly reconstruct which frames were dropped, so we take
% the first nWand frames from validFrames.
if numel(validFrames) > nWand
    fprintf('Trimming validFrames from %d to %d to match BA wand points\n', ...
        numel(validFrames), nWand);
    validFrames = validFrames(1:nWand);
elseif numel(validFrames) < nWand
    fprintf('[WARNING] validFrames (%d) < nWand (%d) — mismatch!\n', ...
        numel(validFrames), nWand);
end

fprintf('Final validFrames: %d (BA wand points: %d)\n', numel(validFrames), nWand);

% Save back to bundle_adjustment.mat
ba.validFrames = validFrames;
ba.obs = obs;
save(fullfile(resultsDir, 'bundle_adjustment.mat'), '-struct', 'ba');
fprintf('Patched bundle_adjustment.mat with validFrames and obs\n');