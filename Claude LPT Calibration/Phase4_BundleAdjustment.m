% =========================================================================
% Phase4_BundleAdjustment.m
% =========================================================================
% Hybrid initialisation: essential matrix rotations (reliable cross-side
% pairs) + geometric derivation (same-side cameras) + C_approx positions.
%
%   STAGE 1 — Initialise poses
%     a) E-matrix for best cross-side pairs (Cam1-Cam3, Cam1-Cam4)
%     b) Geometric rotation for same-side cameras (Cam2 from Cam1)
%     c) C_approx for all positions
%
%   STAGE 2 — Joint pinhole BA (cameras + wand)
%
%   STAGE 3 — Full refractive BA with Huber loss
% =========================================================================

clear; clc;
LPT_Config;

fprintf('=== Phase 4: Bundle Adjustment ===\n\n');

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

n = struct('air',cfg.n_air,'glass',cfg.n_glass,'water',cfg.n_water);
fprintf('Focal lengths: [%.0f, %.0f, %.0f, %.0f] px\n\n', ...
    Ks{1}(1,1), Ks{2}(1,1), Ks{3}(1,1), Ks{4}(1,1));

% Build observations
obs = cell(nFrames,1);
for f = 1:nFrames
    obs{f} = NaN(nCams,4);
    for c = 1:nCams
        det = detections{c}(f);
        if det.valid, obs{f}(c,:) = [det.pts(1,:), det.pts(2,:)]; end
    end
end

% ---- LED A/B correspondence fix ----
% The brightness-based LED ordering from Phase 2 is inconsistent across
% cameras. Exhaustive search of all 16 swap combinations found that
% swapping cameras 1 and 4 maximises cross-side epipolar inliers (3819).
% Swap: [uA vA uB vB] -> [uB vB uA vA]
swapCams = [1, 4];  % cameras whose A/B assignment should be flipped
fprintf('--- LED correspondence fix: swapping A/B for cameras [%s] ---\n', num2str(swapCams));
for f = 1:nFrames
    for c = swapCams
        if ~any(isnan(obs{f}(c,:)))
            obs{f}(c,:) = [obs{f}(c,3:4), obs{f}(c,1:2)];
        end
    end
end

% ---- Camera selection ----
% Exclude cameras with poor detections from the BA. They can be added
% later via PnP against the calibrated 3D wand points.
excludeCams = [];  % set to [] to use all cameras
if ~isempty(excludeCams)
    fprintf('--- Excluding cameras [%s] from BA ---\n', num2str(excludeCams));
    for f = 1:nFrames
        for c = excludeCams
            obs{f}(c,:) = NaN(1,4);
        end
    end
end

% ---- Detection quality filter ----
% Reject individual camera detections where the wand pixel length is
% implausible. This catches noise blobs, mismatched detections, and
% extreme wand orientations (nearly end-on).
minWandPx = 100;   % minimum wand length in pixels
maxWandPx = 800;   % maximum wand length in pixels
fprintf('--- Detection quality filter (wand length: %d-%d px) ---\n', minWandPx, maxWandPx);

nRejected = zeros(nCams,1);
nTotal    = zeros(nCams,1);
for f = 1:nFrames
    for c = 1:nCams
        if ~any(isnan(obs{f}(c,:)))
            nTotal(c) = nTotal(c) + 1;
            wandPxLen = sqrt((obs{f}(c,1)-obs{f}(c,3))^2 + (obs{f}(c,2)-obs{f}(c,4))^2);
            if wandPxLen < minWandPx || wandPxLen > maxWandPx
                obs{f}(c,:) = NaN(1,4);  % invalidate this detection
                nRejected(c) = nRejected(c) + 1;
            end
        end
    end
end

for c = 1:nCams
    fprintf('  Cam %d: rejected %d / %d detections (%.1f%%)\n', ...
        c, nRejected(c), nTotal(c), 100*nRejected(c)/max(nTotal(c),1));
end
fprintf('  Total rejected: %d / %d\n\n', sum(nRejected), sum(nTotal));

% Report filtered wand pixel lengths
fprintf('  Filtered wand image lengths (px):\n');
for c = 1:nCams
    lens = [];
    for f = 1:nFrames
        if ~any(isnan(obs{f}(c,:)))
            lens(end+1) = sqrt((obs{f}(c,1)-obs{f}(c,3))^2+(obs{f}(c,2)-obs{f}(c,4))^2); %#ok<AGROW>
        end
    end
    if ~isempty(lens)
        fprintf('    Cam %d: %.0f +/- %.0f (range %.0f-%.0f, n=%d)\n', ...
            c, mean(lens), std(lens), min(lens), max(lens), numel(lens));
    end
end
fprintf('\n');

allValid = find(arrayfun(@(f) sum(~any(isnan(obs{f}),2))>=2, (1:nFrames)'));
maxFrames = 400;
if numel(allValid)>maxFrames
    validFrames = allValid(round(linspace(1,numel(allValid),maxFrames)));
else
    validFrames = allValid;
end
nValid = numel(validFrames);
fprintf('Frames for BA: %d / %d\n\n', nValid, nFrames);


% =========================================================================
%  STAGE 1: HYBRID INITIALISATION
% =========================================================================
fprintf('========================================\n');
fprintf('  STAGE 1: Hybrid initialisation\n');
fprintf('========================================\n\n');

% Identify sides
leftCams=[]; rightCams=[];
for c=1:nCams
    if cfg.walls(c).normal(1)<0, leftCams(end+1)=c; else, rightCams(end+1)=c; end %#ok<AGROW>
end
fprintf('Left: [%s], Right: [%s]\n', num2str(leftCams), num2str(rightCams));

% --- Compute essential matrices for ALL cross-side pairs ---
fprintf('\nComputing essential matrices for cross-side pairs...\n');

E_results = struct();
nEResults = 0;

for cL = leftCams
    for cR = rightCams
        % Collect correspondences
        pts1_px=[]; pts2_px=[];
        for f=1:nFrames
            if ~any(isnan(obs{f}(cL,:))) && ~any(isnan(obs{f}(cR,:)))
                pts1_px=[pts1_px;obs{f}(cL,1:2);obs{f}(cL,3:4)]; %#ok<AGROW>
                pts2_px=[pts2_px;obs{f}(cR,1:2);obs{f}(cR,3:4)]; %#ok<AGROW>
            end
        end
        nShared = size(pts1_px,1)/2;
        if nShared < 30, continue; end

        pts1_ud = undistortPoints(pts1_px,Ks{cL},distCoeffs{cL});
        pts2_ud = undistortPoints(pts2_px,Ks{cR},distCoeffs{cR});
        meanF = mean([Ks{cL}(1,1),Ks{cR}(1,1)]);
        pts1_n = [(pts1_ud(:,1)-Ks{cL}(1,3))/Ks{cL}(1,1), ...
                  (pts1_ud(:,2)-Ks{cL}(2,3))/Ks{cL}(2,2)];
        pts2_n = [(pts2_ud(:,1)-Ks{cR}(1,3))/Ks{cR}(1,1), ...
                  (pts2_ud(:,2)-Ks{cR}(2,3))/Ks{cR}(2,2)];

        thr = (2.0/meanF)^2;
        [E,inlMask] = estimateEssentialRANSAC(pts1_n,pts2_n,3000,thr);
        nInl = sum(inlMask);

        if nInl >= 20
            [R,t] = recoverPoseFromE(E,pts1_n(inlMask,:),pts2_n(inlMask,:));
            rotAngle = rad2deg(acos(max(-1,min(1,(trace(R)-1)/2))));

            % Only accept if rotation is physically reasonable (>90 deg for cross-side)
            if rotAngle > 90
                nEResults = nEResults+1;
                E_results(nEResults).cL = cL;
                E_results(nEResults).cR = cR;
                E_results(nEResults).R = R;
                E_results(nEResults).t = t;
                E_results(nEResults).nInliers = nInl;
                E_results(nEResults).rotAngle = rotAngle;
                fprintf('  Cam%d-Cam%d: %.1f deg, %d inliers @2px [OK]\n', cL, cR, rotAngle, nInl);
            else
                fprintf('  Cam%d-Cam%d: %.1f deg, %d inliers @2px [REJECTED: angle too small]\n', ...
                    cL, cR, rotAngle, nInl);
            end
        else
            fprintf('  Cam%d-Cam%d: only %d inliers @2px [SKIPPED]\n', cL, cR, nInl);
        end
    end
end

% --- Build camera rotations ---
% Start ALL cameras with geometric look-at rotations from C_approx.
% Then optionally upgrade cross-side cameras with E-matrix if available.

fprintf('\nBuilding camera rotations...\n');
R_world = cell(nCams,1);
t_world = cell(nCams,1);

% Step 1: All cameras get look-at rotations
for c = 1:nCams
    C_c = cfg.initPoses(c).C_approx(:);
    tgt = cfg.initPoses(c).target(:);
    R_world{c} = lookAtRotation(C_c, tgt);
    t_world{c} = -R_world{c} * C_c;
    fprintf('  Cam %d: look-at rotation\n', c);
end

% Step 2: For right-side cameras, upgrade with E-matrix if a good one exists
for c = rightCams
    bestIdx = 0; bestInl = 0;
    for cL_check = leftCams
        for i = 1:nEResults
            if E_results(i).cL == cL_check && E_results(i).cR == c
                if E_results(i).nInliers > bestInl
                    bestInl = E_results(i).nInliers;
                    bestIdx = i;
                end
            end
        end
    end

    if bestIdx > 0 && bestInl >= 50  % only use E-matrix if we have strong evidence
        % The E-matrix gives R_c relative to the left camera
        % We need to convert to world rotation
        cL_ref = E_results(bestIdx).cL;
        R_rel = E_results(bestIdx).R;  % R such that x_c = R * x_L + t
        R_world{c} = R_rel * R_world{cL_ref};
        C_c = cfg.initPoses(c).C_approx(:);
        t_world{c} = -R_world{c} * C_c;
        fprintf('  Cam %d: UPGRADED with E-matrix from Cam%d (%.1f deg, %d inliers)\n', ...
            c, cL_ref, E_results(bestIdx).rotAngle, bestInl);
    end
end

% Print all camera centres
fprintf('\n  Camera centres:\n');
for c = 1:nCams
    C = -R_world{c}'*t_world{c};
    fprintf('    Cam %d: [%.4f %.4f %.4f]\n', c, C(1), C(2), C(3));
end

% --- Reprojection check ---
fprintf('\n  Initial reprojection check (pinhole):\n');
sampleF = validFrames(round(linspace(1,nValid,min(100,nValid))));
camErr = zeros(nCams,1); camCnt = zeros(nCams,1);
for f = sampleF'
    seen = find(~any(isnan(obs{f}(:,1:2)),2));
    if numel(seen)<2, continue; end
    for pi=1:numel(seen), for pj=pi+1:numel(seen)
        ca=seen(pi); cb=seen(pj);
        if cfg.walls(ca).normal(1)==cfg.walls(cb).normal(1), continue; end
        P1=Ks{ca}*[R_world{ca},t_world{ca}]; P2=Ks{cb}*[R_world{cb},t_world{cb}];
        for pt=1:2
            if pt==1, uv1=obs{f}(ca,1:2)';uv2=obs{f}(cb,1:2)';
            else, uv1=obs{f}(ca,3:4)';uv2=obs{f}(cb,3:4)'; end
            A=[uv1(1)*P1(3,:)-P1(1,:);uv1(2)*P1(3,:)-P1(2,:);
               uv2(1)*P2(3,:)-P2(1,:);uv2(2)*P2(3,:)-P2(2,:)];
            [~,~,V]=svd(A); X=V(1:3,end)/V(4,end);
            for ck=seen'
                Xc=R_world{ck}*X+t_world{ck};
                if Xc(3)>0.01
                    u_p=Ks{ck}(1,1)*Xc(1)/Xc(3)+Ks{ck}(1,3);
                    v_p=Ks{ck}(2,2)*Xc(2)/Xc(3)+Ks{ck}(2,3);
                    if pt==1, e=sqrt((u_p-obs{f}(ck,1))^2+(v_p-obs{f}(ck,2))^2);
                    else, e=sqrt((u_p-obs{f}(ck,3))^2+(v_p-obs{f}(ck,4))^2); end
                    camErr(ck)=camErr(ck)+e; camCnt(ck)=camCnt(ck)+1;
                end
            end
        end
        break;
    end; end
end
for c=1:nCams
    if camCnt(c)>0
        fprintf('    Cam %d: mean=%.1f px (%d obs)\n', c, camErr(c)/camCnt(c), camCnt(c));
    end
end
overallInit = sum(camErr)/max(sum(camCnt),1);
fprintf('    Overall: %.1f px\n', overallInit);


% =========================================================================
%  SEED WAND POSITIONS
% =========================================================================
fprintf('\n--- Triangulating wand positions ---\n');
wandPts_init = zeros(nValid,6); goodFrame = true(nValid,1);
for fi=1:nValid
    f=validFrames(fi);
    seen=find(~any(isnan(obs{f}(:,1:2)),2));
    if numel(seen)<2, goodFrame(fi)=false; continue; end
    bestXA=[]; bestXB=[]; bestScore=inf;
    for pi=1:numel(seen), for pj=pi+1:numel(seen)
        ca=seen(pi);cb=seen(pj);
        pA=struct('R',R_world{ca},'t',t_world{ca});
        pB=struct('R',R_world{cb},'t',t_world{cb});
        XA_t=dltTriangulate(obs{f}(ca,1:2)',obs{f}(cb,1:2)',Ks{ca},Ks{cb},pA,pB);
        XB_t=dltTriangulate(obs{f}(ca,3:4)',obs{f}(cb,3:4)',Ks{ca},Ks{cb},pA,pB);
        d=norm(XA_t-XB_t);
        if d<0.01||d>1.0, continue; end
        lengthErr=abs(d-cfg.wandLength)/cfg.wandLength;
        isCross=(cfg.walls(ca).normal(1)~=cfg.walls(cb).normal(1));
        score=lengthErr+ternary(~isCross,0.5,0);
        if score<bestScore, bestScore=score; bestXA=XA_t; bestXB=XB_t; end
    end; end
    if isempty(bestXA), goodFrame(fi)=false; continue; end
    mid=(bestXA+bestXB)/2; wdir=bestXB-bestXA; wn=norm(wdir);
    if wn<1e-8, wdir=[1;0;0]; else, wdir=wdir/wn; end
    wandPts_init(fi,:) = [directionToAxisAngle(wdir), mid'];
end
validFrames=validFrames(goodFrame); wandPts_init=wandPts_init(goodFrame,:);
nValid=numel(validFrames);
fprintf('  Good triangulations: %d\n', nValid);

% Wand length check
triLens=zeros(nValid,1);
for fi=1:nValid
    f=validFrames(fi); seen=find(~any(isnan(obs{f}(:,1:2)),2));
    if numel(seen)>=2
        ca=seen(1);cb=seen(2);
        XA=dltTriangulate(obs{f}(ca,1:2)',obs{f}(cb,1:2)',Ks{ca},Ks{cb},...
            struct('R',R_world{ca},'t',t_world{ca}),struct('R',R_world{cb},'t',t_world{cb}));
        XB=dltTriangulate(obs{f}(ca,3:4)',obs{f}(cb,3:4)',Ks{ca},Ks{cb},...
            struct('R',R_world{ca},'t',t_world{ca}),struct('R',R_world{cb},'t',t_world{cb}));
        triLens(fi)=norm(XA-XB);
    end
end
vl=triLens(triLens>0.01);
if ~isempty(vl)
    fprintf('  Wand lengths: %.2f +/- %.2f mm (true: %.1f)\n', ...
        mean(vl)*1000, std(vl)*1000, cfg.wandLength*1000);
end


% =========================================================================
%  STAGE 2: JOINT PINHOLE BA (two-phase: small then full)
% =========================================================================
fprintf('\n========================================\n');
fprintf('  STAGE 2: Joint pinhole BA\n');
fprintf('========================================\n\n');

camParams0=zeros(nCams,6);
for c=1:nCams
    aa=rotm2axang(R_world{c});
    camParams0(c,:)=[aa(1:3)*aa(4), t_world{c}'];
end

% --- Phase 2a: Small subset to refine cameras ---
nSmall = min(100, nValid);
smallIdx = round(linspace(1, nValid, nSmall));
smallFrames = validFrames(smallIdx);
smallWand   = wandPts_init(smallIdx,:);

fprintf('  Phase 2a: %d frames, %d parameters\n', nSmall, nCams*6 + nSmall*6);

x0_small = packParams(camParams0, smallWand, nCams, nSmall);

% Use regularised cost that prevents cameras from crossing their glass walls
camSideSign = zeros(nCams,1);  % +1 if camera should be at X > wall, -1 if X < wall
camWallX    = zeros(nCams,1);
for c = 1:nCams
    if cfg.walls(c).normal(1) < 0  % left wall, camera on negative-X side
        camSideSign(c) = -1;
        camWallX(c) = cfg.walls(c).point(1) + cfg.walls(c).thickness * abs(cfg.walls(c).normal(1));
    else  % right wall, camera on positive-X side
        camSideSign(c) = +1;
        camWallX(c) = cfg.walls(c).point(1) + cfg.walls(c).thickness * abs(cfg.walls(c).normal(1));
    end
end

costFun_small = @(x) pinholeBA_cost_bounded(x,obs,smallFrames,Ks,distCoeffs,...
    cfg.wandLength,nCams,camSideSign,camWallX,10.0);
r0s = costFun_small(x0_small);
fprintf('  Initial pinhole RMS: %.2f px\n', rms(r0s)/sqrt(2));

opts2a=optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
    'MaxIterations',300,'FunctionTolerance',1e-10,'StepTolerance',1e-10,'Display','iter');
[x_small,~,r_small,ef_small]=lsqnonlin(costFun_small,x0_small,[],[],opts2a);
fprintf('\n  Phase 2a exit: %d, RMS: %.4f px\n', ef_small, rms(r_small)/sqrt(2));

% Extract refined camera params
[camParams_refined, ~] = unpackParams(x_small, nCams, nSmall);

% Update camera poses
for c=1:nCams
    R_world{c}=axang2rotm(camParams_refined(c,1:3));
    t_world{c}=camParams_refined(c,4:6)';
end

fprintf('\n  Camera centres after Phase 2a:\n');
for c=1:nCams
    C=-R_world{c}'*t_world{c};
    fprintf('    Cam %d: [%.4f %.4f %.4f]\n',c,C(1),C(2),C(3));
end

% --- Re-triangulate ALL wand positions with refined cameras ---
fprintf('\n  Re-triangulating wand positions with refined cameras...\n');
wandPts_retri2 = zeros(nValid,6); goodFrame2 = true(nValid,1);
for fi=1:nValid
    f=validFrames(fi);
    seen=find(~any(isnan(obs{f}(:,1:2)),2));
    if numel(seen)<2, goodFrame2(fi)=false; continue; end
    bestXA=[]; bestXB=[]; bestScore=inf;
    for pi=1:numel(seen), for pj=pi+1:numel(seen)
        ca=seen(pi);cb=seen(pj);
        pA=struct('R',R_world{ca},'t',t_world{ca});
        pB=struct('R',R_world{cb},'t',t_world{cb});
        XA_t=dltTriangulate(obs{f}(ca,1:2)',obs{f}(cb,1:2)',Ks{ca},Ks{cb},pA,pB);
        XB_t=dltTriangulate(obs{f}(ca,3:4)',obs{f}(cb,3:4)',Ks{ca},Ks{cb},pA,pB);
        d=norm(XA_t-XB_t);
        if d<0.01||d>1.0, continue; end
        lengthErr=abs(d-cfg.wandLength)/cfg.wandLength;
        isCross=(cfg.walls(ca).normal(1)~=cfg.walls(cb).normal(1));
        score=lengthErr+ternary(~isCross,0.5,0);
        if score<bestScore, bestScore=score; bestXA=XA_t; bestXB=XB_t; end
    end; end
    if isempty(bestXA), goodFrame2(fi)=false; continue; end
    mid=(bestXA+bestXB)/2; wdir=bestXB-bestXA; wn=norm(wdir);
    if wn<1e-8, wdir=[1;0;0]; else, wdir=wdir/wn; end
    wandPts_retri2(fi,:) = [directionToAxisAngle(wdir), mid'];
end
validFrames=validFrames(goodFrame2); wandPts_retri2=wandPts_retri2(goodFrame2,:);
nValid=numel(validFrames);
fprintf('  Good re-triangulations: %d\n', nValid);

vl2=zeros(nValid,1);
for fi=1:nValid
    f=validFrames(fi); seen=find(~any(isnan(obs{f}(:,1:2)),2));
    if numel(seen)>=2
        ca=seen(1);cb=seen(2);
        XA=dltTriangulate(obs{f}(ca,1:2)',obs{f}(cb,1:2)',Ks{ca},Ks{cb},...
            struct('R',R_world{ca},'t',t_world{ca}),struct('R',R_world{cb},'t',t_world{cb}));
        XB=dltTriangulate(obs{f}(ca,3:4)',obs{f}(cb,3:4)',Ks{ca},Ks{cb},...
            struct('R',R_world{ca},'t',t_world{ca}),struct('R',R_world{cb},'t',t_world{cb}));
        vl2(fi)=norm(XA-XB);
    end
end
vl2g=vl2(vl2>0.01);
if ~isempty(vl2g)
    fprintf('  Wand lengths: %.2f +/- %.2f mm (true: %.1f)\n', ...
        mean(vl2g)*1000, std(vl2g)*1000, cfg.wandLength*1000);
end

% --- Phase 2b: Full BA with all frames ---
fprintf('\n  Phase 2b: full BA with %d frames, %d parameters\n', nValid, nCams*6+nValid*6);

x0_full = packParams(camParams_refined, wandPts_retri2, nCams, nValid);
costFun_full = @(x) pinholeBA_cost_bounded(x,obs,validFrames,Ks,distCoeffs,...
    cfg.wandLength,nCams,camSideSign,camWallX,5.0);
r0f = costFun_full(x0_full);
fprintf('  Initial pinhole RMS: %.2f px\n', rms(r0f)/sqrt(2));

opts2b=optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
    'MaxIterations',300,'FunctionTolerance',1e-10,'StepTolerance',1e-10,'Display','iter');
[x_ph,~,r_ph,ef_ph]=lsqnonlin(costFun_full,x0_full,[],[],opts2b);
fprintf('\n  Phase 2b exit: %d, RMS: %.4f px\n', ef_ph, rms(r_ph)/sqrt(2));

[camParams2,wandPts2]=unpackParams(x_ph,nCams,nValid);
for c=1:nCams, R_world{c}=axang2rotm(camParams2(c,1:3)); t_world{c}=camParams2(c,4:6)'; end

fprintf('\n  Camera centres after pinhole BA:\n');
for c=1:nCams
    C=-R_world{c}'*t_world{c};
    fprintf('    Cam %d: [%.4f %.4f %.4f]\n',c,C(1),C(2),C(3));
end


% =========================================================================
%  STAGE 2.5: REFRACTIVE RE-TRIANGULATION
% =========================================================================
fprintf('\n--- Refractive re-triangulation ---\n');
wallGeoms=cell(nCams,1); camPoses=cell(nCams,1);
for c=1:nCams
    wallGeoms{c}=struct('normal',cfg.walls(c).normal(:),'point',cfg.walls(c).point(:),...
        'thickness',cfg.walls(c).thickness);
    camPoses{c}=struct('R',R_world{c},'t',t_world{c});
end
wandPts_retri=zeros(nValid,6); goodF2=true(nValid,1); reconLens=zeros(nValid,1);
for fi=1:nValid
    f=validFrames(fi);
    oA=NaN(nCams,2);oB=NaN(nCams,2);
    for c=1:nCams
        if ~any(isnan(obs{f}(c,:))), oA(c,:)=obs{f}(c,1:2); oB(c,:)=obs{f}(c,3:4); end
    end
    rA=triangulateRefractive(oA,camPoses,Ks,distCoeffs,wallGeoms,n,cfg.recon);
    rB=triangulateRefractive(oB,camPoses,Ks,distCoeffs,wallGeoms,n,cfg.recon);
    if any(isnan(rA.pos))||any(isnan(rB.pos)), goodF2(fi)=false; continue; end
    reconLens(fi)=norm(rA.pos-rB.pos);
    maxE=max([rA.reprojErrors(~isnan(rA.reprojErrors));rB.reprojErrors(~isnan(rB.reprojErrors))]);
    if isempty(maxE)||maxE>15, goodF2(fi)=false; continue; end
    mid=(rA.pos+rB.pos)/2; wd=rB.pos-rA.pos; wn=norm(wd);
    if wn<1e-8,wd=[1;0;0];else,wd=wd/wn;end
    wandPts_retri(fi,:)=[directionToAxisAngle(wd),mid'];
end
gl=reconLens(goodF2);
if ~isempty(gl)
    med=median(gl);mad_v=1.4826*median(abs(gl-med));thr=max(3*mad_v,0.005);
    for fi=1:nValid, if goodF2(fi)&&abs(reconLens(fi)-med)>thr, goodF2(fi)=false;end;end
end
nRetri=sum(goodF2);
fprintf('  Refractive: %d / %d frames\n',nRetri,nValid);
if nRetri>0
    fprintf('  Wand length: %.2f +/- %.2f mm\n',mean(reconLens(goodF2))*1000,std(reconLens(goodF2))*1000);
end
for fi=1:nValid, if ~goodF2(fi), wandPts_retri(fi,:)=wandPts2(fi,:); end; end
if nRetri<20, wandPts_retri=wandPts2; end


% =========================================================================
%  STAGE 3: REFRACTIVE BA
% =========================================================================
fprintf('\n========================================\n');
fprintf('  STAGE 3: Refractive BA\n');
fprintf('========================================\n\n');

x0=packParams(camParams2,wandPts_retri,nCams,nValid);
huberThresh=5.0;
costFun_ref=@(x) refractiveBA_cost(x,obs,validFrames,Ks,distCoeffs,cfg.walls,n,cfg.wandLength,nCams,huberThresh);
r0r=costFun_ref(x0);
fprintf('  Initial refractive RMS: %.4f px\n',rms(r0r)/sqrt(2));
[~,iPC]=computePerCameraError(x0,obs,validFrames,Ks,distCoeffs,cfg.walls,n,cfg.wandLength,nCams);
for c=1:nCams, fprintf('    Cam %d: %.2f px\n',c,iPC(c)); end

costHistory=[];
nP=numel(x0); lb=-inf(nP,1); ub=inf(nP,1);
for c=1:nCams, idx=(c-1)*6+(4:6); lb(idx)=-10; ub(idx)=10; end

for pass=1:2
    if pass==1
        fprintf('\n  --- Pass 1 (Huber=%.1f) ---\n',huberThresh);
        mI=200;fT=1e-8;xT=1e-8;
    else
        huberThresh=2.0;
        costFun_ref=@(x)refractiveBA_cost(x,obs,validFrames,Ks,distCoeffs,cfg.walls,n,cfg.wandLength,nCams,huberThresh);
        fprintf('\n  --- Pass 2 (Huber=%.1f) ---\n',huberThresh);
        mI=cfg.ba.maxIter;fT=cfg.ba.fTol;xT=cfg.ba.xTol;
    end
    clear recordCost
    opts3=optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
        'MaxIterations',mI,'FunctionTolerance',fT,'StepTolerance',xT,...
        'Display','iter-detailed','OutputFcn',@recordCost,'SpecifyObjectiveGradient',false);
    [x0,~,residuals,exitflag,~,~,J_sparse]=lsqnonlin(costFun_ref,x0,lb,ub,opts3);
    fprintf('  Pass %d exit: %d, RMS: %.4f px\n',pass,exitflag,rms(residuals)/sqrt(2));
end
x_opt=x0;
fprintf('\n  Final RMS: %.4f px\n',rms(residuals)/sqrt(2));


% =========================================================================
%  RESULTS
% =========================================================================
[poses_opt,wandPts_opt]=unpackParams(x_opt,nCams,nValid);
nResid=numel(residuals);nParam=numel(x_opt);
if nResid>nParam
    s2=sum(residuals.^2)/(nResid-nParam);JtJ=full(J_sparse'*J_sparse);
    covMatrix=s2*pinv(JtJ);paramStd=sqrt(max(0,diag(covMatrix)));
else, covMatrix=[];paramStd=inf(nParam,1);end

poseStd=reshape(paramStd(1:nCams*6),6,nCams)';
fprintf('\nPose uncertainties:\n');
for c=1:nCams
    fprintf('  Cam %d rot:[%.2e %.2e %.2e] trans:[%.2e %.2e %.2e]\n',c,...
        poseStd(c,1),poseStd(c,2),poseStd(c,3),poseStd(c,4),poseStd(c,5),poseStd(c,6));
end

[~,perCamErr]=computePerCameraError(x_opt,obs,validFrames,Ks,distCoeffs,cfg.walls,n,cfg.wandLength,nCams);
fprintf('\nPer-camera RMS:\n');
for c=1:nCams, fprintf('  Cam %d: %.4f px\n',c,perCamErr(c)); end

[~,frameErrors]=computePerFrameError(x_opt,obs,validFrames,Ks,distCoeffs,cfg.walls,n,cfg.wandLength,nCams);
outlierFrames=validFrames(frameErrors>cfg.ba.reprojThresh);
fprintf('\nOutliers (>%.1f px): %d / %d\n',cfg.ba.reprojThresh,numel(outlierFrames),nValid);

cameras=struct();
for c=1:nCams
    cameras(c).R=axang2rotm(poses_opt(c,1:3));cameras(c).t=poses_opt(c,4:6)';
    cameras(c).K=Ks{c};cameras(c).distCoeffs=distCoeffs{c};
    cameras(c).wall=cfg.walls(c);cameras(c).reprojRMS=perCamErr(c);cameras(c).poseStd=poseStd(c,:);
end
save(fullfile(resultsDir,'extrinsics.mat'),'cameras','cfg');
save(fullfile(resultsDir,'bundle_adjustment.mat'),'cameras','wandPts_opt','residuals','covMatrix','frameErrors','outlierFrames','costHistory','cfg','validFrames','obs');
fprintf('\n[Phase 4] Saved to %s\n',resultsDir);

% Figures
fig=figure('Visible','off');
if ~isempty(costHistory),plot(costHistory,'b-','LineWidth',1.5);end
xlabel('Iter');ylabel('Cost');title('BA Convergence');grid on;
saveas(fig,fullfile(resultsDir,'ba_convergence.png'));close(fig);
fig2=figure('Visible','off');
histogram(frameErrors,40,'FaceColor',[0.2 0.5 0.8]);
xline(cfg.ba.reprojThresh,'r--','LineWidth',2);
xlabel('RMS (px)');ylabel('Count');title('Per-frame Errors');grid on;
saveas(fig2,fullfile(resultsDir,'ba_frame_errors.png'));close(fig2);
fig3=figure('Visible','off','Position',[0 0 900 500]);
subplot(1,2,1);
for c=1:nCams
    C0=cfg.initPoses(c).C_approx(:);
    plot3(C0(1),C0(2),C0(3),'ro','MarkerSize',10,'LineWidth',2);hold on;
    C1=-cameras(c).R'*cameras(c).t;
    plot3(C1(1),C1(2),C1(3),'b^','MarkerSize',10,'LineWidth',2);
    plot3([C0(1) C1(1)],[C0(2) C1(2)],[C0(3) C1(3)],'k--');
    text(C1(1),C1(2),C1(3)+0.03,sprintf('C%d',c));
end
grid on;axis equal;xlabel('X');ylabel('Y');zlabel('Z');title('init(o) vs opt(^)');
subplot(1,2,2);
if ~isempty(frameErrors)
    histogram(frameErrors,40,'FaceColor',[0.2 0.7 0.3]);
    xlabel('RMS(px)');ylabel('Count');title(sprintf('Med=%.3f px',median(frameErrors)));grid on;
end
sgtitle('Diagnostics');saveas(fig3,fullfile(resultsDir,'ba_diagnostics.png'));close(fig3);
fprintf('[Phase 4] Figures saved.\n');


% =========================================================================
function stop=recordCost(~,oV,state)
    persistent h;if strcmp(state,'init'),h=[];elseif strcmp(state,'iter'),h(end+1)=oV.resnorm;
    elseif strcmp(state,'done'),assignin('base','costHistory',h);end;stop=false;end

% =========================================================================
function R = lookAtRotation(C, target)
% Build rotation matrix for a camera at C looking at target.
    zAx = target - C; zAx = zAx/norm(zAx);
    up = [0;0;1];
    if abs(dot(zAx,up)) > 0.9, up = [0;1;0]; end
    xAx = cross(zAx,up); xAx = xAx/norm(xAx);
    yAx = cross(zAx,xAx);
    R = [xAx, yAx, zAx]';
end

% =========================================================================
%  ESSENTIAL MATRIX
% =========================================================================
function [E_best,inlierMask]=estimateEssentialRANSAC(p1,p2,nIter,thresh)
    nPts=size(p1,1);bestN=0;E_best=eye(3);inlierMask=false(nPts,1);
    for it=1:nIter
        idx=randperm(nPts,min(8,nPts));
        E_c=eightPointE(p1(idx,:),p2(idx,:));
        if isempty(E_c),continue;end
        errs=sampsonError(E_c,p1,p2);inl=errs<thresh;nI=sum(inl);
        if nI>bestN,bestN=nI;inlierMask=inl;E_best=E_c;end
    end
    if sum(inlierMask)>=8
        E_best=eightPointE(p1(inlierMask,:),p2(inlierMask,:));
        errs=sampsonError(E_best,p1,p2);inlierMask=errs<thresh;
    end
end
function E=eightPointE(p1,p2)
    nn=size(p1,1);A=zeros(nn,9);
    for i=1:nn,A(i,:)=[p2(i,1)*p1(i,1),p2(i,1)*p1(i,2),p2(i,1),...
        p2(i,2)*p1(i,1),p2(i,2)*p1(i,2),p2(i,2),p1(i,1),p1(i,2),1];end
    [~,~,V]=svd(A);E=reshape(V(:,end),3,3)';
    [U,S,V]=svd(E);S(3,3)=0;S(1,1)=(S(1,1)+S(2,2))/2;S(2,2)=S(1,1);E=U*S*V';
end
function errs=sampsonError(E,p1,p2)
    nn=size(p1,1);errs=zeros(nn,1);
    for i=1:nn,pp1=[p1(i,:),1]';pp2=[p2(i,:),1]';Ep1=E*pp1;Etp2=E'*pp2;
        num=(pp2'*E*pp1)^2;den=Ep1(1)^2+Ep1(2)^2+Etp2(1)^2+Etp2(2)^2;
        errs(i)=num/max(den,1e-12);end
end
function [R,t]=recoverPoseFromE(E,p1,p2)
    [U,~,V]=svd(E);W=[0 -1 0;1 0 0;0 0 1];
    Rs={U*W*V',U*W*V',U*W'*V',U*W'*V'};ts={U(:,3),-U(:,3),U(:,3),-U(:,3)};
    for i=1:4,if det(Rs{i})<0,Rs{i}=-Rs{i};ts{i}=-ts{i};end;end
    bestS=-1;R=Rs{1};t=ts{1};
    for i=1:4,P1_=eye(3,4);P2_=[Rs{i},ts{i}];sc=0;
        for k=1:min(size(p1,1),100)
            A=[p1(k,1)*P1_(3,:)-P1_(1,:);p1(k,2)*P1_(3,:)-P1_(2,:);
               p2(k,1)*P2_(3,:)-P2_(1,:);p2(k,2)*P2_(3,:)-P2_(2,:)];
            [~,~,Vt]=svd(A);X=Vt(1:3,end)/Vt(4,end);X2=Rs{i}*X+ts{i};
            if X(3)>0&&X2(3)>0,sc=sc+1;end;end
        if sc>bestS,bestS=sc;R=Rs{i};t=ts{i};end;end
end

function pts_out=undistortPoints(pts_in,K,dc)
    if numel(dc)<5,dc(end+1:5)=0;end;pts_out=pts_in;
    for i=1:size(pts_in,1)
        xd_t=(pts_in(i,1)-K(1,3))/K(1,1);yd_t=(pts_in(i,2)-K(2,3))/K(2,2);
        xn=xd_t;yn=yd_t;
        for it=1:20,[xd,yd]=applyDistortion(xn,yn,dc);xn=xn+(xd_t-xd);yn=yn+(yd_t-yd);end
        pts_out(i,1)=xn*K(1,1)+K(1,3);pts_out(i,2)=yn*K(2,2)+K(2,3);end
end

% =========================================================================
%  BA COSTS
% =========================================================================
function r=pinholeBA_cost(x,obs,vF,Ks,dc,wL,nC)
    nV=numel(vF);[p,w]=unpackParams(x,nC,nV);
    cR=cell(nC,1);cT=cell(nC,1);
    for c=1:nC,cR{c}=axang2rotm(p(c,1:3));cT{c}=p(c,4:6)';end
    r=[];
    for fi=1:nV,f=vF(fi);[XA,XB]=decodeWand(w(fi,:),wL);
        for c=1:nC
            if any(isnan(obs{f}(c,:))),continue;end
            K=Ks{c};R=cR{c};t=cT{c};
            XAc=R*XA+t;XBc=R*XB+t;
            if XAc(3)<=0.01||XBc(3)<=0.01,r=[r;50;50;50;50];continue;end %#ok<AGROW>
            xnA=XAc(1)/XAc(3);ynA=XAc(2)/XAc(3);[xdA,ydA]=applyDistortion(xnA,ynA,dc{c});
            uAp=[K(1,1)*xdA+K(1,3);K(2,2)*ydA+K(2,3)];
            xnB=XBc(1)/XBc(3);ynB=XBc(2)/XBc(3);[xdB,ydB]=applyDistortion(xnB,ynB,dc{c});
            uBp=[K(1,1)*xdB+K(1,3);K(2,2)*ydB+K(2,3)];
            r=[r;uAp-obs{f}(c,1:2)';uBp-obs{f}(c,3:4)']; %#ok<AGROW>
        end;end
end

function r=pinholeBA_cost_bounded(x,obs,vF,Ks,dc,wL,nC,camSideSign,camWallX,regW)
% Pinhole BA cost with camera-side constraint.
% Penalises cameras that cross to the wrong side of their glass wall.
    nV=numel(vF);[p,w]=unpackParams(x,nC,nV);
    cR=cell(nC,1);cT=cell(nC,1);
    for c=1:nC,cR{c}=axang2rotm(p(c,1:3));cT{c}=p(c,4:6)';end
    r=[];
    % Reprojection residuals
    for fi=1:nV,f=vF(fi);[XA,XB]=decodeWand(w(fi,:),wL);
        for c=1:nC
            if any(isnan(obs{f}(c,:))),continue;end
            K=Ks{c};R=cR{c};t=cT{c};
            XAc=R*XA+t;XBc=R*XB+t;
            if XAc(3)<=0.01||XBc(3)<=0.01,r=[r;50;50;50;50];continue;end %#ok<AGROW>
            xnA=XAc(1)/XAc(3);ynA=XAc(2)/XAc(3);[xdA,ydA]=applyDistortion(xnA,ynA,dc{c});
            uAp=[K(1,1)*xdA+K(1,3);K(2,2)*ydA+K(2,3)];
            xnB=XBc(1)/XBc(3);ynB=XBc(2)/XBc(3);[xdB,ydB]=applyDistortion(xnB,ynB,dc{c});
            uBp=[K(1,1)*xdB+K(1,3);K(2,2)*ydB+K(2,3)];
            r=[r;uAp-obs{f}(c,1:2)';uBp-obs{f}(c,3:4)']; %#ok<AGROW>
        end;end
    % Camera-side penalty: camera must stay on correct side of wall
    if regW > 0
        for c=1:nC
            C_c = -cR{c}'*cT{c};
            % How far is the camera on the wrong side? (negative = wrong side)
            margin = camSideSign(c) * (C_c(1) - camWallX(c));
            if margin < 0.05  % within 5cm of wall or past it
                penalty = regW * (0.05 - margin);
                r=[r; penalty; penalty; penalty]; %#ok<AGROW>
            else
                r=[r; 0; 0; 0]; %#ok<AGROW>
            end
        end
    end
end

function r=refractiveBA_cost(x,obs,vF,Ks,dc,walls,n,wL,nC,hT)
    nV=numel(vF);[p,w]=unpackParams(x,nC,nV);
    cP=buildCamPoses(p,nC);wG=buildWallGeoms(walls,nC);r=[];
    for fi=1:nV,f=vF(fi);[XA,XB]=decodeWand(w(fi,:),wL);
        for c=1:nC
            if any(isnan(obs{f}(c,:))),continue;end
            uAo=obs{f}(c,1:2)';uBo=obs{f}(c,3:4)';
            uAp=projectPointRefractive(XA,cP{c},Ks{c},dc{c},wG{c},n);
            uBp=projectPointRefractive(XB,cP{c},Ks{c},dc{c},wG{c},n);
            if any(isnan(uAp))||any(isnan(uBp))
                Cc=-cP{c}.R'*cP{c}.t;pen=20+5*norm(Cc);
                r=[r;pen;pen;pen;pen];continue;end %#ok<AGROW>
            raw=[uAp-uAo;uBp-uBo];
            for ri=1:numel(raw),a=abs(raw(ri));
                if a>hT,raw(ri)=sign(raw(ri))*sqrt(2*hT*a-hT^2);end;end
            r=[r;raw]; %#ok<AGROW>
    end;end
end

% =========================================================================
%  ERROR COMPUTATION
% =========================================================================
function [rt,pcRMS]=computePerCameraError(x,obs,vF,Ks,dc,walls,n,wL,nC)
    nV=numel(vF);[p,w]=unpackParams(x,nC,nV);cP=buildCamPoses(p,nC);wG=buildWallGeoms(walls,nC);
    eS=zeros(nC,1);ct=zeros(nC,1);
    for fi=1:nV,f=vF(fi);[XA,XB]=decodeWand(w(fi,:),wL);
        for c=1:nC,if any(isnan(obs{f}(c,:))),continue;end
            uA=projectPointRefractive(XA,cP{c},Ks{c},dc{c},wG{c},n);
            uB=projectPointRefractive(XB,cP{c},Ks{c},dc{c},wG{c},n);
            if any(isnan(uA))||any(isnan(uB)),continue;end
            eS(c)=eS(c)+norm(uA-obs{f}(c,1:2)')^2+norm(uB-obs{f}(c,3:4)')^2;ct(c)=ct(c)+2;
    end;end;pcRMS=sqrt(eS./max(ct,1));rt=sqrt(sum(eS)/max(sum(ct),1));
end
function [rt,fRMS]=computePerFrameError(x,obs,vF,Ks,dc,walls,n,wL,nC)
    nV=numel(vF);[p,w]=unpackParams(x,nC,nV);cP=buildCamPoses(p,nC);wG=buildWallGeoms(walls,nC);
    fRMS=zeros(nV,1);
    for fi=1:nV,f=vF(fi);[XA,XB]=decodeWand(w(fi,:),wL);eS=0;ct=0;
        for c=1:nC,if any(isnan(obs{f}(c,:))),continue;end
            uA=projectPointRefractive(XA,cP{c},Ks{c},dc{c},wG{c},n);
            uB=projectPointRefractive(XB,cP{c},Ks{c},dc{c},wG{c},n);
            if any(isnan(uA))||any(isnan(uB)),continue;end
            eS=eS+norm(uA-obs{f}(c,1:2)')^2+norm(uB-obs{f}(c,3:4)')^2;ct=ct+2;end
        fRMS(fi)=sqrt(eS/max(ct,1));end;rt=rms(fRMS);
end

% =========================================================================
%  UTILITIES
% =========================================================================
function x=packParams(p,w,nC,nF),x=[reshape(p(1:nC,:)',[],1);reshape(w(1:nF,:)',[],1)];end
function [p,w]=unpackParams(x,nC,nF),p=reshape(x(1:nC*6),6,nC)';w=reshape(x(nC*6+1:end),6,nF)';end
function cP=buildCamPoses(p,nC),cP=cell(nC,1);for c=1:nC,cP{c}=struct('R',axang2rotm(p(c,1:3)),'t',p(c,4:6)');end;end
function wG=buildWallGeoms(w,nC),wG=cell(nC,1);for c=1:nC,wG{c}=struct('normal',w(c).normal(:),'point',w(c).point(:),'thickness',w(c).thickness);end;end
function [XA,XB]=decodeWand(wp,wL)
    aa=wp(1:3);ang=norm(aa);if ang<1e-10,wd=[1;0;0];else,wd=axang2rotm([aa/ang,ang])*[1;0;0];end
    wd=wd/norm(wd);mid=wp(4:6)';XA=mid-(wL/2)*wd;XB=mid+(wL/2)*wd;end
function av=directionToAxisAngle(w),w=w(:)/norm(w);ref=[1;0;0];ax=cross(ref,w);axn=norm(ax);
    if axn<1e-8,if dot(ref,w)>0,av=[0 0 0];else,av=[0 pi 0];end
    else,ang=atan2(axn,dot(ref,w));av=(ax/axn*ang)';end;end
function s=ternary(c,a,b),if c,s=a;else,s=b;end;end
function R=axang2rotm(aa)
    if numel(aa)==3,angle=norm(aa);if angle<1e-10,R=eye(3);return;end;ax=aa(:)/angle;
    else,ax=aa(1:3)';angle=aa(4);end
    c=cos(angle);s=sin(angle);t=1-c;x=ax(1);y=ax(2);z=ax(3);
    R=[t*x*x+c,t*x*y-s*z,t*x*z+s*y;t*x*y+s*z,t*y*y+c,t*y*z-s*x;t*x*z-s*y,t*y*z+s*x,t*z*z+c];end
function aa=rotm2axang(R),angle=acos(max(-1,min(1,(trace(R)-1)/2)));
    if abs(angle)<1e-10,aa=[1 0 0 0];return;end
    ax=[(R(3,2)-R(2,3));(R(1,3)-R(3,1));(R(2,1)-R(1,2))]/(2*sin(angle));aa=[ax(:)',angle];end
function X=dltTriangulate(uv1,uv2,K1,K2,p1,p2)
    P1=K1*[p1.R,p1.t];P2=K2*[p2.R,p2.t];
    A=[uv1(1)*P1(3,:)-P1(1,:);uv1(2)*P1(3,:)-P1(2,:);uv2(1)*P2(3,:)-P2(1,:);uv2(2)*P2(3,:)-P2(2,:)];
    [~,~,V]=svd(A);X=V(1:3,end)/V(4,end);end