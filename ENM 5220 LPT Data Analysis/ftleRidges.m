function [Xr, Sr, out] = ftleRidges(V, xg, yg, zg, opts)
% FTLERIDGES  Second-derivative (height) ridges of an FTLE field -> LCS.
%
%   Ridge = codim-1 locus of FTLE maxima. With H = grad grad V and its
%   ascending eigenpairs (lam1<=lam2<=lam3, e1..e3), a ridge point satisfies
%     (i)   grad V _|_ e1        : |gradV.e1| / |gradV| < alignTol
%     (ii)  lam1 << 0            : strong transverse concavity (adaptive)
%     (iii) V above a floor      : ignore weak background straining
%   Only sign-invariant quantities are used (|gradV.e1|, lam1), so the
%   eigenvector sign ambiguity never enters. This LOCATES structures by
%   geometry rather than by an absolute sigma level -- robust to T, box, h.
%
%   Feed it your ensemble-mean field (persistent LCS) or a single-window
%   field. Forward-FTLE ridges = repelling LCS; backward = attracting.
%
%   opts (all optional)
%     .smoothVox   Gaussian pre-smoothing width in VOXELS before differencing
%                  (default 1.0). Hessian of a noisy field is very noisy; keep
%                  >=1. This is the ridge resolution knob.
%     .sigmaThrPct percentile of valid V used as floor (iii)      (default 50)
%     .alignTol    orthogonality tolerance for (i), in [0,1]      (default 0.15)
%     .concavPct   keep voxels whose -lam1 exceeds this percentile
%                  among candidates -- adaptive threshold for (ii)(default 50)
%
%   OUTPUTS
%     Xr  (P x 3) ridge-point coordinates
%     Sr  (P x 1) smoothed FTLE at those points (for colouring)
%     out struct: .ridgeMask (3D logical), .Vsmooth, .validE, .sigThr

    if nargin < 5, opts = struct; end
    smoothVox   = getdef(opts,'smoothVox',1.0);
    sigmaThrPct = getdef(opts,'sigmaThrPct',50);
    alignTol    = getdef(opts,'alignTol',0.15);
    concavPct   = getdef(opts,'concavPct',50);

    [Xg,Yg,Zg] = meshgrid(xg,yg,zg);
    hx = xg(2)-xg(1); hy = yg(2)-yg(1); hz = zg(2)-zg(1);

    % restrict reported ridges to the genuinely-covered interior
    validE = erode1(isfinite(V));

    % fill holes so derivatives/smoothing are defined, then Gaussian smooth
    Vf = nanfill3(V, 8);
    if smoothVox > 0
        sz = 2*ceil(2*smoothVox)+1;
        Vf = smooth3(Vf, 'gaussian', sz, smoothVox);
    end

    % gradient and symmetric Hessian. MATLAB order: 1st out is d/dx (dim2, hx)
    [Gx,Gy,Gz]    = gradient(Vf, hx,hy,hz);
    [Gxx,Gxy,Gxz] = gradient(Gx, hx,hy,hz);   % Vxx Vxy Vxz
    [~,  Gyy,Gyz] = gradient(Gy, hx,hy,hz);   %     Vyy Vyz
    [~,  ~,  Gzz] = gradient(Gz, hx,hy,hz);   %         Vzz
    gmag = sqrt(Gx.^2 + Gy.^2 + Gz.^2) + eps;

    % candidate voxels: covered interior above the sigma floor
    sigThr = prctile(Vf(validE), sigmaThrPct);
    cand   = validE & (Vf > sigThr);
    ci     = find(cand);

    lam1  = nan(numel(ci),1);
    align = nan(numel(ci),1);
    for q = 1:numel(ci)
        p = ci(q);
        H = [Gxx(p) Gxy(p) Gxz(p);
             Gxy(p) Gyy(p) Gyz(p);
             Gxz(p) Gyz(p) Gzz(p)];
        [Vv,Dd] = eig(H);
        [ev,ord] = sort(diag(Dd),'ascend');   % force ascending
        e1 = Vv(:,ord(1));
        g  = [Gx(p); Gy(p); Gz(p)];
        lam1(q)  = ev(1);
        align(q) = abs(g.'*e1) / gmag(p);
    end

    % adaptive concavity threshold among candidates
    concav = -lam1;
    isRidge = false(numel(ci),1);
    if any(concav > 0)
        lamCut  = prctile(concav(concav>0), concavPct);
        isRidge = (concav > lamCut) & (lam1 < 0) & (align < alignTol);
    end

    ridgeIdx = ci(isRidge);
    Xr = [Xg(ridgeIdx), Yg(ridgeIdx), Zg(ridgeIdx)];
    Sr = Vf(ridgeIdx);

    ridgeMask = false(size(Vf));  ridgeMask(ridgeIdx) = true;
    out = struct('ridgeMask',ridgeMask,'Vsmooth',Vf,'validE',validE,'sigThr',sigThr);
end

% ---------- local helpers ----------
function v = getdef(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function e = erode1(M)
% true only where the voxel and all 6 face-neighbours are true
    e = M;
    e(2:end,:,:)   = e(2:end,:,:)   & M(1:end-1,:,:);
    e(1:end-1,:,:) = e(1:end-1,:,:) & M(2:end,:,:);
    e(:,2:end,:)   = e(:,2:end,:)   & M(:,1:end-1,:);
    e(:,1:end-1,:) = e(:,1:end-1,:) & M(:,2:end,:);
    e(:,:,2:end)   = e(:,:,2:end)   & M(:,:,1:end-1);
    e(:,:,1:end-1) = e(:,:,1:end-1) & M(:,:,2:end);
    e([1 end],:,:) = false; e(:,[1 end],:) = false; e(:,:,[1 end]) = false;
end

function Vf = nanfill3(V, K)
% fill NaNs by iterated 6-neighbour averaging (wrap-free), K passes
    Vf = V;
    for it = 1:K
        nanm = isnan(Vf);
        if ~any(nanm(:)), break; end
        A = Vf; A(nanm) = 0;  F = double(~nanm);
        S = zeros(size(Vf)); C = zeros(size(Vf));
        S(1:end-1,:,:) = S(1:end-1,:,:) + A(2:end,:,:);   C(1:end-1,:,:) = C(1:end-1,:,:) + F(2:end,:,:);
        S(2:end,:,:)   = S(2:end,:,:)   + A(1:end-1,:,:); C(2:end,:,:)   = C(2:end,:,:)   + F(1:end-1,:,:);
        S(:,1:end-1,:) = S(:,1:end-1,:) + A(:,2:end,:);   C(:,1:end-1,:) = C(:,1:end-1,:) + F(:,2:end,:);
        S(:,2:end,:)   = S(:,2:end,:)   + A(:,1:end-1,:); C(:,2:end,:)   = C(:,2:end,:)   + F(:,1:end-1,:);
        S(:,:,1:end-1) = S(:,:,1:end-1) + A(:,:,2:end);   C(:,:,1:end-1) = C(:,:,1:end-1) + F(:,:,2:end);
        S(:,:,2:end)   = S(:,:,2:end)   + A(:,:,1:end-1); C(:,:,2:end)   = C(:,:,2:end)   + F(:,:,1:end-1);
        fillv = S ./ max(C,1);
        take  = nanm & (C>0);
        Vf(take) = fillv(take);
    end
    if any(isnan(Vf(:)))          % isolated leftovers -> global mean
        Vf(isnan(Vf)) = mean(V(isfinite(V)));
    end
end