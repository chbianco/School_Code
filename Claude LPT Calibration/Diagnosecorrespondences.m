% =========================================================================
% DiagnoseSwapCombinations.m
% =========================================================================
% Tests all 16 possible LED A/B swap combinations (each camera can be
% swapped or not) and reports which combination gives the best overall
% epipolar consistency across all camera pairs.
% =========================================================================

clear; clc;
LPT_Config;

fprintf('=== LED Swap Combination Search ===\n\n');

resultsDir = cfg.resultsDir;
phase1 = load(fullfile(resultsDir, 'phase1_report.mat'));
phase2 = load(fullfile(resultsDir, 'wand_detections.mat'));

detections = phase2.detections;
nCams      = cfg.nCams;
nFrames    = phase2.nFrames;

Ks = cell(nCams,1); distCoeffs = cell(nCams,1);
for c = 1:nCams
    Ks{c} = phase1.intrinsics{c}.K;
    dc = [phase1.calib{c}.params.RadialDistortion, ...
          phase1.calib{c}.params.TangentialDistortion];
    if numel(dc)<5, dc(end+1:5)=0; end
    distCoeffs{c} = dc;
end

% Build base observations
obs_base = cell(nFrames,1);
for f = 1:nFrames
    obs_base{f} = NaN(nCams,4);
    for c = 1:nCams
        det = detections{c}(f);
        if det.valid
            obs_base{f}(c,:) = [det.pts(1,:), det.pts(2,:)];
        end
    end
end

% Identify cross-side pairs (these are the important ones)
crossPairs = [];
for c1 = 1:nCams
    for c2 = c1+1:nCams
        if cfg.walls(c1).normal(1) ~= cfg.walls(c2).normal(1)
            crossPairs = [crossPairs; c1, c2]; %#ok<AGROW>
        end
    end
end
samePairs = [];
for c1 = 1:nCams
    for c2 = c1+1:nCams
        if cfg.walls(c1).normal(1) == cfg.walls(c2).normal(1)
            samePairs = [samePairs; c1, c2]; %#ok<AGROW>
        end
    end
end

fprintf('Cross-side pairs: ');
for i = 1:size(crossPairs,1), fprintf('%d-%d ', crossPairs(i,1), crossPairs(i,2)); end
fprintf('\nSame-side pairs: ');
for i = 1:size(samePairs,1), fprintf('%d-%d ', samePairs(i,1), samePairs(i,2)); end
fprintf('\n\n');

% Pre-compute normalised points for all pairs (both orders) to avoid
% recomputing in the loop
fprintf('Pre-computing correspondences for all pairs...\n');
pairData = struct();
allPairs = [crossPairs; samePairs];
for pi = 1:size(allPairs,1)
    c1 = allPairs(pi,1); c2 = allPairs(pi,2);
    
    pts1 = []; pts2A = []; pts2B = [];
    sharedCount = 0;
    for f = 1:nFrames
        if ~any(isnan(obs_base{f}(c1,:))) && ~any(isnan(obs_base{f}(c2,:)))
            sharedCount = sharedCount + 1;
            pts1  = [pts1;  obs_base{f}(c1,1:2); obs_base{f}(c1,3:4)]; %#ok<AGROW>
            pts2A = [pts2A; obs_base{f}(c2,1:2); obs_base{f}(c2,3:4)]; %#ok<AGROW>  % same order
            pts2B = [pts2B; obs_base{f}(c2,3:4); obs_base{f}(c2,1:2)]; %#ok<AGROW>  % swapped
        end
    end
    
    pairData(pi).c1 = c1;
    pairData(pi).c2 = c2;
    pairData(pi).nShared = sharedCount;
    pairData(pi).pts1_raw = pts1;
    pairData(pi).pts2A_raw = pts2A;  % c2 same order
    pairData(pi).pts2B_raw = pts2B;  % c2 swapped
end

% Test all 16 swap combinations
fprintf('\nTesting all 16 swap combinations...\n\n');

meanF = mean([Ks{1}(1,1), Ks{2}(1,1), Ks{3}(1,1), Ks{4}(1,1)]);
thr = (5.0 / meanF)^2;

bestTotalInliers = 0;
bestCombo = [0 0 0 0];
allResults = zeros(16, 6);  % [swap1 swap2 swap3 swap4 totalCrossInliers totalAllInliers]

fprintf('%-20s  ', 'Swap combo');
for pi = 1:size(allPairs,1)
    fprintf('C%d-C%d  ', allPairs(pi,1), allPairs(pi,2));
end
fprintf('  Cross   All\n');
fprintf('%s\n', repmat('-', 1, 120));

comboIdx = 0;
for s1 = 0:1
    for s2 = 0:1
        for s3 = 0:1
            for s4 = 0:1
                comboIdx = comboIdx + 1;
                swaps = [s1, s2, s3, s4];
                
                totalCrossInliers = 0;
                totalAllInliers = 0;
                
                label = sprintf('[%d%d%d%d]', s1, s2, s3, s4);
                swappedCams = find(swaps);
                if isempty(swappedCams)
                    label = [label, ' (none)     '];
                else
                    label = [label, sprintf(' swap C%s', num2str(swappedCams))];
                    label = [label, repmat(' ', 1, max(0, 20-length(label)))];
                end
                fprintf('%-20s  ', label);
                
                for pi = 1:size(allPairs,1)
                    c1 = allPairs(pi,1);
                    c2 = allPairs(pi,2);
                    
                    if pairData(pi).nShared < 20
                        fprintf(' skip  ');
                        continue;
                    end
                    
                    % Determine which point set to use based on swap state
                    % pts1 comes from c1. If c1 is swapped, we need to swap pts1.
                    % pts2 comes from c2. pts2A = same order, pts2B = swapped.
                    
                    pts1_raw = pairData(pi).pts1_raw;
                    
                    % If c1 is swapped, swap the A/B in pts1
                    if swaps(c1)
                        % pts1 is [A1; B1; A1; B1; ...], swap to [B1; A1; B1; A1; ...]
                        pts1_use = pts1_raw;
                        for k = 1:2:size(pts1_raw,1)-1
                            pts1_use(k,:)   = pts1_raw(k+1,:);
                            pts1_use(k+1,:) = pts1_raw(k,:);
                        end
                    else
                        pts1_use = pts1_raw;
                    end
                    
                    % If c2 is swapped, use pts2B (swapped), else pts2A (same)
                    if swaps(c2)
                        pts2_use = pairData(pi).pts2B_raw;
                    else
                        pts2_use = pairData(pi).pts2A_raw;
                    end
                    
                    % Normalise
                    p1n = [(pts1_use(:,1)-Ks{c1}(1,3))/Ks{c1}(1,1), ...
                           (pts1_use(:,2)-Ks{c1}(2,3))/Ks{c1}(2,2)];
                    p2n = [(pts2_use(:,1)-Ks{c2}(1,3))/Ks{c2}(1,1), ...
                           (pts2_use(:,2)-Ks{c2}(2,3))/Ks{c2}(2,2)];
                    
                    % Estimate F with RANSAC
                    [~, inlMask] = estimateFundRANSAC(p1n, p2n, 500, thr);
                    nInl = sum(inlMask);
                    
                    isCross = cfg.walls(c1).normal(1) ~= cfg.walls(c2).normal(1);
                    if isCross
                        totalCrossInliers = totalCrossInliers + nInl;
                    end
                    totalAllInliers = totalAllInliers + nInl;
                    
                    fprintf('%4d    ', nInl);
                end
                
                fprintf('  %5d  %5d', totalCrossInliers, totalAllInliers);
                
                if totalCrossInliers > bestTotalInliers
                    bestTotalInliers = totalCrossInliers;
                    bestCombo = swaps;
                    fprintf('  <-- BEST');
                end
                fprintf('\n');
                
                allResults(comboIdx,:) = [s1 s2 s3 s4 totalCrossInliers totalAllInliers];
            end
        end
    end
end

fprintf('%s\n', repmat('-', 1, 120));
fprintf('\nBest combination: [%d %d %d %d]', bestCombo);
swappedCams = find(bestCombo);
if isempty(swappedCams)
    fprintf(' — no swaps needed\n');
else
    fprintf(' — swap cameras: %s\n', num2str(swappedCams));
end
fprintf('Cross-side inliers: %d\n', bestTotalInliers);

% Also report top 3
[~, sortIdx] = sort(allResults(:,5), 'descend');
fprintf('\nTop 3 combinations:\n');
for i = 1:min(3, size(allResults,1))
    r = allResults(sortIdx(i),:);
    swCams = find(r(1:4));
    if isempty(swCams), swStr = 'none'; else, swStr = num2str(swCams); end
    fprintf('  [%d%d%d%d] swap [%s]: cross=%d, all=%d\n', ...
        r(1), r(2), r(3), r(4), swStr, r(5), r(6));
end

fprintf('\n=== Done ===\n');
fprintf('\nTo use in Phase4, set:  swapCams = [%s];\n', num2str(find(bestCombo)));


% =========================================================================
function [F_best, inlierMask] = estimateFundRANSAC(p1, p2, nIter, thresh)
    nPts = size(p1,1);
    bestN = 0; F_best = eye(3); inlierMask = false(nPts,1);
    for it = 1:nIter
        idx = randperm(nPts, min(8, nPts));
        F = eightPointF(p1(idx,:), p2(idx,:));
        if isempty(F), continue; end
        errs = fundError(F, p1, p2);
        inl = errs < thresh;
        nI = sum(inl);
        if nI > bestN
            bestN = nI; inlierMask = inl; F_best = F;
        end
    end
    if sum(inlierMask) >= 8
        F_best = eightPointF(p1(inlierMask,:), p2(inlierMask,:));
        errs = fundError(F_best, p1, p2);
        inlierMask = errs < thresh;
    end
end

function F = eightPointF(p1, p2)
    nn = size(p1,1);
    A = zeros(nn, 9);
    for i = 1:nn
        A(i,:) = [p2(i,1)*p1(i,1), p2(i,1)*p1(i,2), p2(i,1), ...
                  p2(i,2)*p1(i,1), p2(i,2)*p1(i,2), p2(i,2), ...
                  p1(i,1), p1(i,2), 1];
    end
    [~,~,V] = svd(A);
    F = reshape(V(:,end), 3, 3)';
    [U,S,V] = svd(F);
    S(3,3) = 0;
    F = U * S * V';
end

function errs = fundError(F, p1, p2)
    nn = size(p1,1);
    errs = zeros(nn,1);
    for i = 1:nn
        pp1 = [p1(i,:), 1]';
        pp2 = [p2(i,:), 1]';
        Fp1 = F * pp1;
        Ftp2 = F' * pp2;
        num = (pp2' * F * pp1)^2;
        den = Fp1(1)^2 + Fp1(2)^2 + Ftp2(1)^2 + Ftp2(2)^2;
        errs(i) = num / max(den, 1e-12);
    end
end