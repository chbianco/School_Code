function [Xp, sigma, info] = ftleMeshfree(tracks, i0, m, opts)
% FTLEMESHFREE  Mesh-free finite-time Lyapunov exponent from measured LPT tracks.
%
%   Estimates the flow-map gradient F = dPhi_{t0}^{t0+T}/dX at each surviving
%   particle by a least-squares fit of final-to-initial neighbour separations
%   (the particles ARE the flow map -- nothing is integrated), then
%
%       sigma = (1/|T|) * log( sigma_max(F) )
%
%   with sigma_max the largest singular value of F (== sqrt(lambda_max(F'F))).
%
%   INPUTS
%     tracks : struct with fields x,y,z  (nT x nTracks, NaN-padded) and t (nT x 1) [SI]
%     i0     : start-frame row index into 1..nT
%     m      : window length in FRAMES. m>0 forward (repelling LCS);
%              m<0 backward (attracting LCS). i1 = i0 + m.
%     opts   : (optional) struct
%       .k           neighbours for the LS gradient (default 12; >=4 required)
%       .radius      if set, use fixed-radius neighbours instead of k-NN [m].
%                    Fixed radius = uniform filter scale; k-NN adapts to local
%                    seeding density (denser -> smaller effective scale).
%       .requireSpan require every row in [i0,i1] finite, not just endpoints
%                    (default false; endpoints are all the flow map needs)
%       .condTol     reject a neighbourhood when svd ratio s3/s1 < condTol,
%                    i.e. neighbours do not span 3D (default 1e-3)
%
%   OUTPUTS
%     Xp    : (P x 3) initial positions of accepted seeds            [m]
%     sigma : (P x 1) FTLE at those seeds                            [1/s]
%     info  : struct of diagnostics (T, nValid, nRejected, nSeeds, k)
%
%   Requires knnsearch/rangesearch (Statistics & Machine Learning Toolbox).

    if nargin < 4, opts = struct; end
    k        = getdef(opts,'k',12);
    radius   = getdef(opts,'radius',[]);
    reqSpan  = getdef(opts,'requireSpan',false);
    condTol  = getdef(opts,'condTol',1e-3);

    nT = numel(tracks.t);
    i1 = i0 + m;
    assert(i0>=1 && i0<=nT && i1>=1 && i1<=nT, 'window [i0,i1]=[%d,%d] out of range 1..%d', i0,i1,nT);
    T = abs(tracks.t(i1) - tracks.t(i0));
    assert(T > 0, 'zero-length window');
    lo = min(i0,i1); hi = max(i0,i1);

    % --- endpoint positions, nTracks x 3 ---
    X0 = [tracks.x(i0,:); tracks.y(i0,:); tracks.z(i0,:)].';
    X1 = [tracks.x(i1,:); tracks.y(i1,:); tracks.z(i1,:)].';

    % --- surviving set S: finite at both endpoints (or across full span) ---
    if reqSpan
        good = all(isfinite(tracks.x(lo:hi,:)) & isfinite(tracks.y(lo:hi,:)) ...
                 & isfinite(tracks.z(lo:hi,:)), 1).';
    else
        good = all(isfinite(X0),2) & all(isfinite(X1),2);
    end
    S  = find(good);
    Xi = X0(S,:);  Xf = X1(S,:);          % initial / final config of the flow map
    nS = numel(S);
    assert(nS > 4, 'too few surviving particles (%d) for a 3D gradient', nS);

    % --- neighbour lists in the INITIAL configuration (one query for all seeds) ---
    useRadius = ~isempty(radius);
    if useRadius
        nbrCell = rangesearch(Xi, Xi, radius);   % each cell: self + neighbours
    else
        assert(nS > k, 'too few surviving particles (%d) for k=%d', nS, k);
        nbrIdx = knnsearch(Xi, Xi, 'K', k+1);
        nbrIdx = nbrIdx(:,2:end);                % drop self column -> nS x k
    end

    % --- per-seed least-squares deformation gradient ---
    sigma = nan(nS,1);
    nRej  = 0;
    for p = 1:nS
        if useRadius
            q = nbrCell{p};  q(q==p) = [];       % neighbours, excluding self
            if numel(q) < 4, nRej = nRej + 1; continue; end
        else
            q = nbrIdx(p,:);
        end
        dX = (Xi(q,:) - Xi(p,:)).';              % 3 x k  initial separations
        dx = (Xf(q,:) - Xf(p,:)).';              % 3 x k  final   separations
        s  = svd(dX);
        if s(3) < condTol*s(1)                   % neighbourhood not full 3D -> reject
            nRej = nRej + 1;  continue;
        end
        F        = dx / dX;                      % 3x3 LS solve of  F*dX = dx
        sigma(p) = log(max(svd(F))) / T;         % SVD of F: avoids squaring cond(F)
    end

    keep  = isfinite(sigma);
    Xp    = Xi(keep,:);
    sigma = sigma(keep);
    info  = struct('T',T,'i0',i0,'i1',i1,'nValid',nS, ...
                   'nRejected',nRej,'nSeeds',numel(sigma),'k',k);
end

function v = getdef(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end