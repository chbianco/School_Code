% =========================================================================
% DiagnoseEssentialMatrices.m
% =========================================================================
% PURPOSE:
%   Compute essential matrix decomposition for ALL camera pairs and report:
%     - Relative rotation angle
%     - Translation direction
%     - Number of inliers at various thresholds
%     - Whether the geometry is physically consistent
%
%   Run this BEFORE Phase4 to understand what the image data tells us
%   about the camera geometry.
% =========================================================================

clear; clc;
LPT_Config;

fprintf('=== Essential Matrix Diagnostics ===\n\n');

resultsDir = cfg.resultsDir;

% Load data
phase1 = load(fullfile(resultsDir, 'phase1_report.mat'));
phase2 = load(fullfile(resultsDir, 'wand_detections.mat'));

intrinsics = phase1.intrinsics;
detections = phase2.detections;
nCams      = cfg.nCams;
nFrames    = phase2.nFrames;

Ks = cell(nCams,1); distCoeffs = cell(nCams,1);
for c = 1:nCams
    Ks{c} = intrinsics{c}.K;
    dc = [phase1.calib{c}.params.RadialDistortion, ...
          phase1.calib{c}.params.TangentialDistortion];
    if numel(dc)<5, dc(end+1:5)=0; end
    distCoeffs{c} = dc;
end

fprintf('Focal lengths: [%.0f, %.0f, %.0f, %.0f] px\n\n', ...
    Ks{1}(1,1), Ks{2}(1,1), Ks{3}(1,1), Ks{4}(1,1));

% Build observation matrix
obs = cell(nFrames,1);
for f = 1:nFrames
    obs{f} = NaN(nCams,4);
    for c = 1:nCams
        det = detections{c}(f);
        if det.valid, obs{f}(c,:) = [det.pts(1,:), det.pts(2,:)]; end
    end
end

% For each camera pair
fprintf('=====================================================================\n');
fprintf('%-12s %-8s %-10s %-12s %-8s %s\n', ...
    'Pair', 'Shared', 'Inliers', 'Rot (deg)', 'Side', 'Translation dir');
fprintf('=====================================================================\n');

results = [];
nResults = 0;

for c1 = 1:nCams
    for c2 = c1+1:nCams
        % Collect shared frames
        pts1_px = []; pts2_px = [];
        for f = 1:nFrames
            if ~any(isnan(obs{f}(c1,:))) && ~any(isnan(obs{f}(c2,:)))
                pts1_px = [pts1_px; obs{f}(c1,1:2); obs{f}(c1,3:4)]; %#ok<AGROW>
                pts2_px = [pts2_px; obs{f}(c2,1:2); obs{f}(c2,3:4)]; %#ok<AGROW>
            end
        end

        nShared = size(pts1_px,1) / 2;  % number of frames
        if nShared < 20
            fprintf('Cam%d-Cam%d    %4d     skip (too few)\n', c1, c2, nShared);
            continue;
        end

        % Undistort and normalise
        pts1_ud = undistortPoints(pts1_px, Ks{c1}, distCoeffs{c1});
        pts2_ud = undistortPoints(pts2_px, Ks{c2}, distCoeffs{c2});
        pts1_n = [(pts1_ud(:,1)-Ks{c1}(1,3))/Ks{c1}(1,1), ...
                  (pts1_ud(:,2)-Ks{c1}(2,3))/Ks{c1}(2,2)];
        pts2_n = [(pts2_ud(:,1)-Ks{c2}(1,3))/Ks{c2}(1,1), ...
                  (pts2_ud(:,2)-Ks{c2}(2,3))/Ks{c2}(2,2)];

        % Try multiple thresholds
        meanF = mean([Ks{c1}(1,1), Ks{c2}(1,1)]);
        thresholds_px = [2, 5, 10, 20];

        foundResult = false;
        for ti = 1:numel(thresholds_px)
            thr = (thresholds_px(ti) / meanF)^2;  % Sampson error is squared
            [E, inlMask] = estimateEssentialRANSAC(pts1_n, pts2_n, 2000, thr);
            nInl = sum(inlMask);

            if nInl >= 20
                [R, t] = recoverPoseFromE(E, pts1_n(inlMask,:), pts2_n(inlMask,:));
                rotAngle = rad2deg(acos(max(-1,min(1,(trace(R)-1)/2))));

                isCross = cfg.walls(c1).normal(1) ~= cfg.walls(c2).normal(1);
                sideStr = ternary(isCross, 'cross', 'same');

                fprintf('Cam%d-Cam%d    %4d     %4d/%4d   %7.2f     %-6s [%.3f %.3f %.3f]  (@%dpx)\n', ...
                    c1, c2, nShared, nInl, size(pts1_n,1), rotAngle, sideStr, t(1), t(2), t(3), thresholds_px(ti));

                % Store result
                nResults = nResults + 1;
                results(nResults).c1 = c1;
                results(nResults).c2 = c2;
                results(nResults).thresh_px = thresholds_px(ti);
                results(nResults).nShared = nShared;
                results(nResults).nInliers = nInl;
                results(nResults).nTotal = size(pts1_n,1);
                results(nResults).rotAngle = rotAngle;
                results(nResults).tDir = t';
                results(nResults).R = R;
                results(nResults).isCross = isCross;
                foundResult = true;
                break;  % use tightest threshold that gives >=20 inliers
            end
        end

        if ~foundResult
            fprintf('Cam%d-Cam%d    %4d     <20 inliers at all thresholds\n', c1, c2, nShared);
        end
    end
end

fprintf('=====================================================================\n');

% -------------------------------------------------------------------------
% Physical consistency checks
% -------------------------------------------------------------------------
fprintf('\n=== Physical Consistency Checks ===\n\n');

fprintf('Expected geometry:\n');
fprintf('  Cross-side pairs (Cam1-Cam3, Cam1-Cam4, Cam2-Cam3, Cam2-Cam4):\n');
fprintf('    Rotation should be ~150-180 deg (cameras face each other)\n');
fprintf('  Same-side pairs (Cam1-Cam2, Cam3-Cam4):\n');
fprintf('    Rotation should be small, ~5-30 deg (cameras point roughly same direction)\n\n');

% Check each result
for i = 1:nResults
    r = results(i);
    
    ok = true;
    if r.isCross
        if r.rotAngle < 90 || r.rotAngle > 180
            fprintf('  [WARN] Cam%d-Cam%d: cross-side rotation %.1f deg (expected 90-180)\n', ...
                r.c1, r.c2, r.rotAngle);
            ok = false;
        end
    else
        if r.rotAngle > 60
            fprintf('  [WARN] Cam%d-Cam%d: same-side rotation %.1f deg (expected <60)\n', ...
                r.c1, r.c2, r.rotAngle);
            ok = false;
        end
    end
    
    if ok
        fprintf('  [OK]   Cam%d-Cam%d: %.1f deg (%s-side) — consistent\n', ...
            r.c1, r.c2, r.rotAngle, ternary(r.isCross,'cross','same'));
    end
end

% -------------------------------------------------------------------------
% Rotation consistency: R_13 should be close to R_12 * R_23
% -------------------------------------------------------------------------
fprintf('\n=== Rotation Chain Consistency ===\n\n');

% Find the R matrices for specific pairs
R_pairs = cell(nCams);
for i = 1:nResults
    r = results(i);
    R_pairs{r.c1, r.c2} = r.R;
end

% Check: R_13 vs R_12 * R_23 (if all three are available)
if ~isempty(R_pairs{1,3}) && ~isempty(R_pairs{1,2}) && ~isempty(R_pairs{2,3})
    R_13_direct = R_pairs{1,3};
    R_13_chain  = R_pairs{2,3} * R_pairs{1,2};  % R_23 * R_12
    dR = R_13_direct * R_13_chain';
    chainErr = rad2deg(acos(max(-1,min(1,(trace(dR)-1)/2))));
    fprintf('  R_13 vs R_12*R_23: %.2f deg discrepancy\n', chainErr);
    if chainErr < 5
        fprintf('  [OK] Rotation chain is consistent\n');
    else
        fprintf('  [WARN] Rotation chain inconsistency — possible issue with one pair\n');
    end
end

if ~isempty(R_pairs{1,3}) && ~isempty(R_pairs{1,4}) && ~isempty(R_pairs{3,4})
    R_14_direct = R_pairs{1,4};
    R_14_chain  = R_pairs{3,4} * R_pairs{1,3};
    dR = R_14_direct * R_14_chain';
    chainErr = rad2deg(acos(max(-1,min(1,(trace(dR)-1)/2))));
    fprintf('  R_14 vs R_13*R_34: %.2f deg discrepancy\n', chainErr);
    if chainErr < 5
        fprintf('  [OK] Rotation chain is consistent\n');
    else
        fprintf('  [WARN] Rotation chain inconsistency\n');
    end
end

% -------------------------------------------------------------------------
% What initial reprojection error would the E-matrix rotations give?
% -------------------------------------------------------------------------
fprintf('\n=== Predicted Initial Reprojection Error ===\n\n');

% Strategy: use R_13 as anchor, derive other rotations via chains
% R_12 = R_13 * inv(R_23)  (avoid the bad direct R_12)
% R_14 = direct R_14 (good)

% Method A: Direct R from each pair to Cam1
fprintf('--- Method A: Direct E-matrix rotations (R_1c) ---\n');
if ~isempty(R_pairs{1,3}) && ~isempty(R_pairs{1,4})
    testReproject(obs, nFrames, Ks, distCoeffs, cfg, nCams, ...
        R_pairs{1,2}, R_pairs{1,3}, R_pairs{1,4}, 'Direct R_1c');
end

% Method B: Chain-derived rotations
fprintf('\n--- Method B: Chain-derived rotations ---\n');
if ~isempty(R_pairs{1,3}) && ~isempty(R_pairs{2,3})
    % R_12_chain = R_13 * inv(R_23)
    % R_1c means "rotation from cam1 frame to cam_c frame"
    % If R_pairs{a,b} = R_b * R_a', then:
    %   R_12 = R_2 * R_1' = R_23' * R_13
    R_12_chain = R_pairs{2,3}' * R_pairs{1,3};
    R_12_angle = rad2deg(acos(max(-1,min(1,(trace(R_12_chain)-1)/2))));
    fprintf('  Chain R_12 (from R_13*inv(R_23)): %.2f deg\n', R_12_angle);
    
    R_14_direct = R_pairs{1,4};
    if isempty(R_14_direct) && ~isempty(R_pairs{1,3}) && ~isempty(R_pairs{3,4})
        R_14_direct = R_pairs{3,4} * R_pairs{1,3};
        fprintf('  Chain R_14 (from R_34*R_13): used\n');
    end
    
    testReproject(obs, nFrames, Ks, distCoeffs, cfg, nCams, ...
        R_12_chain, R_pairs{1,3}, R_14_direct, 'Chain-derived');
end

% Method C: All rotations from R_13 + R_34 + chain
fprintf('\n--- Method C: All from R_13 anchor ---\n');
if ~isempty(R_pairs{1,3}) && ~isempty(R_pairs{2,3}) && ~isempty(R_pairs{3,4})
    R_12_c = R_pairs{2,3}' * R_pairs{1,3};
    R_13_c = R_pairs{1,3};
    R_14_c = R_pairs{3,4} * R_pairs{1,3};
    
    testReproject(obs, nFrames, Ks, distCoeffs, cfg, nCams, ...
        R_12_c, R_13_c, R_14_c, 'R_13 anchor');
end

fprintf('\n=== Diagnostics Complete ===\n');


% =========================================================================
%  HELPER FUNCTIONS
% =========================================================================
function testReproject(obs, nFrames, Ks, distCoeffs, cfg, nCams, R_12, R_13, R_14, label)
% Test reprojection error with given relative rotations + C_approx positions.
    R_test = cell(nCams,1); t_test = cell(nCams,1);
    
    % Cam 1 at identity rotation, positioned at C_approx
    R_test{1} = eye(3);
    t_test{1} = -eye(3) * cfg.initPoses(1).C_approx(:);
    
    % Cam 2, 3, 4 use relative rotations from Cam 1
    R_rels = {[], R_12, R_13, R_14};  % R_rels{c} = R_1c
    for c = 2:nCams
        if ~isempty(R_rels{c})
            R_test{c} = R_rels{c};  % This IS the rotation of cam c (since cam1 = I)
            C_c = cfg.initPoses(c).C_approx(:);
            t_test{c} = -R_test{c} * C_c;
        end
    end
    
    % Check all cameras assigned
    for c = 1:nCams
        if isempty(R_test{c})
            fprintf('  [SKIP] Cam %d has no rotation estimate\n', c);
            return;
        end
    end
    
    % Compute reprojection error
    sampleFrames = round(linspace(1, nFrames, min(200, nFrames)));
    camErrors = zeros(nCams,1); camCounts = zeros(nCams,1);
    
    for f = sampleFrames
        seen = find(~any(isnan(obs{f}(:,1:2)),2));
        if numel(seen) < 2, continue; end
        
        % Triangulate from best cross-side pair
        done = false;
        for pi = 1:numel(seen)
            for pj = pi+1:numel(seen)
                ca = seen(pi); cb = seen(pj);
                isCross = cfg.walls(ca).normal(1) ~= cfg.walls(cb).normal(1);
                if ~isCross, continue; end
                
                P1 = Ks{ca} * [R_test{ca}, t_test{ca}];
                P2 = Ks{cb} * [R_test{cb}, t_test{cb}];
                
                for pt = 1:2
                    if pt==1, uv1=obs{f}(ca,1:2)'; uv2=obs{f}(cb,1:2)';
                    else,     uv1=obs{f}(ca,3:4)'; uv2=obs{f}(cb,3:4)'; end
                    
                    A=[uv1(1)*P1(3,:)-P1(1,:);uv1(2)*P1(3,:)-P1(2,:);
                       uv2(1)*P2(3,:)-P2(1,:);uv2(2)*P2(3,:)-P2(2,:)];
                    [~,~,V]=svd(A); X=V(1:3,end)/V(4,end);
                    
                    for ck = seen'
                        Xc = R_test{ck}*X + t_test{ck};
                        if Xc(3) > 0.01
                            u_p = Ks{ck}(1,1)*Xc(1)/Xc(3)+Ks{ck}(1,3);
                            v_p = Ks{ck}(2,2)*Xc(2)/Xc(3)+Ks{ck}(2,3);
                            if pt==1, e=sqrt((u_p-obs{f}(ck,1))^2+(v_p-obs{f}(ck,2))^2);
                            else,     e=sqrt((u_p-obs{f}(ck,3))^2+(v_p-obs{f}(ck,4))^2); end
                            camErrors(ck) = camErrors(ck)+e;
                            camCounts(ck) = camCounts(ck)+1;
                        end
                    end
                end
                done = true; break;
            end
            if done, break; end
        end
    end
    
    fprintf('  [%s]\n', label);
    for c = 1:nCams
        if camCounts(c)>0
            fprintf('    Cam %d: mean=%.1f px  (%d obs)\n', c, camErrors(c)/camCounts(c), camCounts(c));
        else
            fprintf('    Cam %d: no observations\n', c);
        end
    end
    overallMean = sum(camErrors)/max(sum(camCounts),1);
    fprintf('    Overall: %.1f px\n', overallMean);
    if overallMean<50, fprintf('    [GOOD] Should converge with joint BA\n');
    elseif overallMean<200, fprintf('    [OK] Rough but may converge\n');
    else, fprintf('    [POOR] Too rough for BA\n'); end
end
function [E_best, inlierMask] = estimateEssentialRANSAC(p1, p2, nIter, thresh)
    nPts=size(p1,1); bestN=0; E_best=eye(3); inlierMask=false(nPts,1);
    for it=1:nIter
        idx=randperm(nPts,min(8,nPts));
        E_c=eightPointE(p1(idx,:),p2(idx,:));
        if isempty(E_c), continue; end
        errs=sampsonError(E_c,p1,p2);
        inl=errs<thresh; nI=sum(inl);
        if nI>bestN, bestN=nI; inlierMask=inl; E_best=E_c; end
    end
    if sum(inlierMask)>=8
        E_best=eightPointE(p1(inlierMask,:),p2(inlierMask,:));
        errs=sampsonError(E_best,p1,p2);
        inlierMask=errs<thresh;
    end
end

function E=eightPointE(p1,p2)
    nn=size(p1,1); A=zeros(nn,9);
    for i=1:nn
        A(i,:)=[p2(i,1)*p1(i,1),p2(i,1)*p1(i,2),p2(i,1),...
                p2(i,2)*p1(i,1),p2(i,2)*p1(i,2),p2(i,2),...
                p1(i,1),p1(i,2),1];
    end
    [~,~,V]=svd(A); E=reshape(V(:,end),3,3)';
    [U,S,V]=svd(E); S(3,3)=0;
    S(1,1)=(S(1,1)+S(2,2))/2; S(2,2)=S(1,1);
    E=U*S*V';
end

function errs=sampsonError(E,p1,p2)
    nn=size(p1,1); errs=zeros(nn,1);
    for i=1:nn
        pp1=[p1(i,:),1]'; pp2=[p2(i,:),1]';
        Ep1=E*pp1; Etp2=E'*pp2;
        num=(pp2'*E*pp1)^2;
        den=Ep1(1)^2+Ep1(2)^2+Etp2(1)^2+Etp2(2)^2;
        errs(i)=num/max(den,1e-12);
    end
end

function [R,t]=recoverPoseFromE(E,p1,p2)
    [U,~,V]=svd(E); W=[0 -1 0;1 0 0;0 0 1];
    Rs={U*W*V',U*W*V',U*W'*V',U*W'*V'};
    ts={U(:,3),-U(:,3),U(:,3),-U(:,3)};
    for i=1:4, if det(Rs{i})<0, Rs{i}=-Rs{i}; ts{i}=-ts{i}; end; end
    bestS=-1; R=Rs{1}; t=ts{1};
    for i=1:4
        P1_=eye(3,4); P2_=[Rs{i},ts{i}]; sc=0;
        for k=1:min(size(p1,1),100)
            A=[p1(k,1)*P1_(3,:)-P1_(1,:);p1(k,2)*P1_(3,:)-P1_(2,:);
               p2(k,1)*P2_(3,:)-P2_(1,:);p2(k,2)*P2_(3,:)-P2_(2,:)];
            [~,~,Vt]=svd(A); X=Vt(1:3,end)/Vt(4,end);
            X2=Rs{i}*X+ts{i};
            if X(3)>0 && X2(3)>0, sc=sc+1; end
        end
        if sc>bestS, bestS=sc; R=Rs{i}; t=ts{i}; end
    end
end

function pts_out=undistortPoints(pts_in,K,dc)
    if numel(dc)<5, dc(end+1:5)=0; end
    pts_out=pts_in;
    for i=1:size(pts_in,1)
        xd_t=(pts_in(i,1)-K(1,3))/K(1,1);
        yd_t=(pts_in(i,2)-K(2,3))/K(2,2);
        xn=xd_t; yn=yd_t;
        for it=1:20
            [xd,yd]=applyDistortion(xn,yn,dc);
            xn=xn+(xd_t-xd); yn=yn+(yd_t-yd);
        end
        pts_out(i,1)=xn*K(1,1)+K(1,3);
        pts_out(i,2)=yn*K(2,2)+K(2,3);
    end
end

function s=ternary(c,a,b)
    if c, s=a; else, s=b; end
end