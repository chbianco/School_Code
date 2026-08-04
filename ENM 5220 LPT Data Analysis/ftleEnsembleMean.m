function [Vmean, Cnt] = ftleEnsembleMean(tracks, i0List, m, xg, yg, zg, opts)
% FTLEENSEMBLEMEAN  Mean FTLE over start frames, on a fixed Eulerian grid.
%
%   For each i0 in i0List the mesh-free FTLE (window m frames, same
%   T=|m|*dt) is computed with ftleMeshfree, interpolated onto the fixed
%   grid meshgrid(xg,yg,zg), and averaged where covered. The result
%   <sigma>(x) marks locations of PERSISTENT stretching, and is a drop-in
%   for your existing isosurface block (same shape/convention as Vg).
%
%   This is the mean of instantaneous FTLE fields -- NOT the FTLE of the
%   mean flow. Accumulation on a fixed grid is required because each
%   window's seed support differs; scattered sigma cannot be averaged
%   elementwise.
%
%   INPUTS
%     tracks  : struct with x,y,z (nT x nTracks, NaN-padded), t (nT x 1)
%     i0List  : vector of start-frame indices. Forward (m>0) needs
%               i0 <= nT-m; backward (m<0) needs i0 >= 1-m. Space them by
%               >~ the decorrelation time (T_L/dt frames) for independence.
%     m       : window length in frames (sign selects repelling/attracting)
%     xg,yg,zg: grid vectors (same you already build for isosurface)
%     opts    : (optional)
%       .ftle     struct passed straight to ftleMeshfree (k / radius / condTol)
%       .minCount require >= this many contributing windows per node (default 1)
%
%   OUTPUTS
%     Vmean : mean FTLE on meshgrid(xg,yg,zg); nodes below minCount are NaN
%     Cnt   : number of windows contributing to each node

    if nargin < 7, opts = struct; end
    ftleOpts = getdef(opts,'ftle',struct());
    minCount = getdef(opts,'minCount',1);

    [Xg,Yg,Zg] = meshgrid(xg,yg,zg);
    Vsum = zeros(size(Xg));
    Cnt  = zeros(size(Xg));

    for jj = 1:numel(i0List)
        i0 = i0List(jj);
        [Xp, sigma] = ftleMeshfree(tracks, i0, m, ftleOpts);
        if isempty(sigma), continue; end
        Fs = scatteredInterpolant(Xp(:,1),Xp(:,2),Xp(:,3), sigma, 'linear','none');
        Vg = Fs(Xg,Yg,Zg);
        ok = ~isnan(Vg);
        Vsum(ok) = Vsum(ok) + Vg(ok);
        Cnt(ok)  = Cnt(ok)  + 1;
        fprintf('  i0=%d (%d/%d): seeds=%d  covered=%d\n', ...
                i0, jj, numel(i0List), numel(sigma), nnz(ok));
    end

    Vmean = Vsum ./ Cnt;
    Vmean(Cnt < minCount) = NaN;
end

function v = getdef(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end