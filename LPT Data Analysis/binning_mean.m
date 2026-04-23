%% Preamble
close all; clc;
clearvars -except tracks t x y z
set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);

%% Load Data
tracks = load('jhtdb_long_LPT.mat').LPT;
t = tracks.t;
x = tracks.x;
y = tracks.y;
z = tracks.z;
[nT, nTracks] = size(x);
fprintf('Loaded: %d samples x %d tracks (%.2f GB per array)\n', ...
        nT, nTracks, nT*nTracks*8/1e9);

%% Subset: keep only the N longest tracks
N = nTracks;  

trackLengths = sum(~isnan(x), 1);
[~, order] = sort(trackLengths, 'descend');
keep = order(1:min(N, nTracks));

x = x(:, keep);
y = y(:, keep);
z = z(:, keep);

% Trim rows to longest surviving track
maxLen = max(sum(~isnan(x), 1));
x = x(1:maxLen, :);
y = y(1:maxLen, :);
z = z(1:maxLen, :);
t = t(1:maxLen);

[nT, nTracks] = size(x);
fprintf('Kept %d tracks, max length %d samples\n', nTracks, nT);

%% ------------GET JHU DATA-------------
%% Validate velocity 
authkey = 'edu.jhu.pha.turbulence.testing-201406';
dataset = 'channel';
variable = 'velocity';
temporal_method = 'pchip';
spatial_method = 'lag4';
spatial_operator = 'field';

% y-grid: 128 Chebyshev points, split into two tiles of 64
ny_dns = 128;
ny_tile = 64;
n_ytiles = ceil(ny_dns / ny_tile);

theta_dns = linspace(0, pi, ny_dns);
gridY_dns = -cos(theta_dns);

% x, z: 8 x 8 = 64 points
nx_dns = 8;
nz_dns = 8;

gridX_dns = linspace(2, 8*pi - 2, nx_dns);
gridZ_dns = linspace(0.5, 3*pi - 0.5, nz_dns);

%Number of time samples
n_times = 100;
times = linspace(0.5, 25, n_times);

% Accumulators
sumU_dns  = zeros(ny_dns, 1);
sumV_dns  = zeros(ny_dns, 1);
sumW_dns  = zeros(ny_dns, 1);
sumUU_dns = zeros(ny_dns, 1);
sumVV_dns = zeros(ny_dns, 1);
sumWW_dns = zeros(ny_dns, 1);
sumUV_dns = zeros(ny_dns, 1);
n_samples = 0;

fprintf('Querying JHTDB for velocity: %d times x %d y-tiles = %d total queries...\n', ...
        n_times, n_ytiles, n_times * n_ytiles);

for it = 1:n_times
    t_query = times(it);
    fprintf('  time %d/%d (t = %.2f)', it, n_times, t_query);

    for yt = 1:n_ytiles
        % y-indices for this tile
        jy_start = (yt - 1) * ny_tile + 1;
        jy_end   = min(yt * ny_tile, ny_dns);
        ny_this  = jy_end - jy_start + 1;
        gridY_tile = gridY_dns(jy_start:jy_end);

        % Build query points for this tile
        n_pts = nx_dns * ny_this * nz_dns;
        pts = zeros(n_pts, 3);
        idx = 0;
        for i = 1:nx_dns
            for j = 1:ny_this
                for k = 1:nz_dns
                    idx = idx + 1;
                    pts(idx, 1) = gridX_dns(i);
                    pts(idx, 2) = gridY_tile(j);
                    pts(idx, 3) = gridZ_dns(k);
                end
            end
        end

        % Query
        vel = getData(authkey, dataset, variable, t_query, ...
                      temporal_method, spatial_method, spatial_operator, pts);

        % Reshape: z fastest, then y, then x
        u3d = reshape(vel(:,1), [nz_dns, ny_this, nx_dns]);
        v3d = reshape(vel(:,2), [nz_dns, ny_this, nx_dns]);
        w3d = reshape(vel(:,3), [nz_dns, ny_this, nx_dns]);

        % Accumulate per y-location
        for jj = 1:ny_this
            jy_global = jy_start + jj - 1;

            u_vals = reshape(u3d(:, jj, :), [], 1);
            v_vals = reshape(v3d(:, jj, :), [], 1);
            w_vals = reshape(w3d(:, jj, :), [], 1);

            sumU_dns(jy_global)  = sumU_dns(jy_global)  + sum(u_vals);
            sumV_dns(jy_global)  = sumV_dns(jy_global)  + sum(v_vals);
            sumW_dns(jy_global)  = sumW_dns(jy_global)  + sum(w_vals);
            sumUU_dns(jy_global) = sumUU_dns(jy_global) + sum(u_vals.^2);
            sumVV_dns(jy_global) = sumVV_dns(jy_global) + sum(v_vals.^2);
            sumWW_dns(jy_global) = sumWW_dns(jy_global) + sum(w_vals.^2);
            sumUV_dns(jy_global) = sumUV_dns(jy_global) + sum(u_vals .* v_vals);
        end

        fprintf('.');
    end

    % Update sample count (same for all y-locations)
    n_samples = n_samples + nx_dns * nz_dns;

    fprintf(' done\n');
end

%Compute DNS statistics
Umean_dns = sumU_dns / n_samples;
Vmean_dns = sumV_dns / n_samples;
Wmean_dns = sumW_dns / n_samples;

uu_dns = sumUU_dns / n_samples - Umean_dns.^2;
vv_dns = sumVV_dns / n_samples - Vmean_dns.^2;
ww_dns = sumWW_dns / n_samples - Wmean_dns.^2;
uv_dns = sumUV_dns / n_samples - Umean_dns .* Vmean_dns;
tke_dns = 0.5 * (uu_dns + vv_dns + ww_dns);

yc_dns = gridY_dns(:);

fprintf('\nDNS statistics computed.\n');
fprintf('y-points: %d, samples per y-location: %d\n', ny_dns, n_samples);
fprintf('Total queries: %d\n', n_times * n_ytiles);


%% Validate pressure 
authkey = 'edu.jhu.pha.turbulence.testing-201406';
dataset = 'channel';
variable = 'pressure';
temporal_method = 'pchip';
spatial_method = 'lag4';
spatial_operator = 'field';

% Accumulator
sumP_dns = zeros(ny_dns, 1);
n_samples_p = 0;

fprintf('Querying JHTDB for pressure: %d times x %d y-tiles = %d total queries...\n', ...
        n_times, n_ytiles, n_times * n_ytiles);

for it = 1:n_times
    t_query = times(it);
    fprintf('  time %d/%d (t = %.2f)', it, n_times, t_query);

    for yt = 1:n_ytiles
        % y-indices for this tile
        jy_start = (yt - 1) * ny_tile + 1;
        jy_end   = min(yt * ny_tile, ny_dns);
        ny_this  = jy_end - jy_start + 1;
        gridY_tile = gridY_dns(jy_start:jy_end);

        % Build query points for this tile
        n_pts = nx_dns * ny_this * nz_dns;
        pts = zeros(n_pts, 3);
        idx = 0;
        for i = 1:nx_dns
            for j = 1:ny_this
                for k = 1:nz_dns
                    idx = idx + 1;
                    pts(idx, 1) = gridX_dns(i);
                    pts(idx, 2) = gridY_tile(j);
                    pts(idx, 3) = gridZ_dns(k);
                end
            end
        end

        % Query — pressure is scalar, returns (n_pts x 1)
        pres = getData(authkey, dataset, variable, t_query, ...
                      temporal_method, spatial_method, spatial_operator, pts);

        % Reshape: z fastest, then y, then x
        p3d = reshape(pres, [nz_dns, ny_this, nx_dns]);

        % Accumulate per y-location
        for jj = 1:ny_this
            jy_global = jy_start + jj - 1;
            p_vals = reshape(p3d(:, jj, :), [], 1);
            sumP_dns(jy_global) = sumP_dns(jy_global) + sum(p_vals);
        end

        fprintf('.');
    end

    n_samples_p = n_samples_p + nx_dns * nz_dns;
    fprintf(' done\n');
end

%DNS mean pressure profile
Pmean_dns = sumP_dns / n_samples_p;

% Subtract centerline value so p=0 at center
[~, ic_dns] = min(abs(yc_dns));
Pmean_dns = Pmean_dns - Pmean_dns(ic_dns);

fprintf('DNS pressure: %d samples per y-location\n', n_samples_p);


%% -----------TRACK DATA------------------

%% Sweep across bin sizes
Xbin_vec = [128];
Ybin_vec = [32];
Zbin_vec = [64]; 

%Preallocate bins for largest Ny
ny_max = max(Ybin_vec);
err_U  = NaN(length(Ybin_vec), ny_max);
err_uu = NaN(length(Ybin_vec), ny_max);
err_vv = NaN(length(Ybin_vec), ny_max);
err_ww = NaN(length(Ybin_vec), ny_max);
err_uv = NaN(length(Ybin_vec), ny_max);
err_P  = NaN(length(Ybin_vec), ny_max);

for bn = 1:length(Xbin_vec)

%Define bin grid
gridX = linspace(0,    8*pi, Xbin_vec(bn) + 1);
gridZ = linspace(0,    3*pi, Zbin_vec(bn) + 1);

ny_bins = Ybin_vec(bn);
theta = linspace(0, pi, ny_bins + 1);
gridY = -cos(theta);   %Cluster bins at walls

nx = numel(gridX) - 1;
ny = numel(gridY) - 1;
nz = numel(gridZ) - 1;

xc = 0.5*(gridX(1:end-1) + gridX(2:end));
yc = 0.5*(gridY(1:end-1) + gridY(2:end));
zc = 0.5*(gridZ(1:end-1) + gridZ(2:end));

%Chunked binning: accumulate sums over track chunks
counts = zeros(nx, ny, nz);
sumU   = zeros(nx, ny, nz);
sumV   = zeros(nx, ny, nz);
sumW   = zeros(nx, ny, nz);

chunkSize = 1000;   % tracks per chunk; tune for memory
dtCol = diff(t(:));

nChunks = ceil(nTracks / chunkSize);
fprintf('Processing %d chunks of up to %d tracks...\n', nChunks, chunkSize);

for c = 1:nChunks
    i0 = (c-1)*chunkSize + 1;
    i1 = min(c*chunkSize, nTracks);

    xc_ = x(:, i0:i1);
    yc_ = y(:, i0:i1);
    zc_ = z(:, i0:i1);

    % Forward differences -> velocities at midpoints
    U = diff(xc_, 1, 1) ./ dtCol;
    V = diff(yc_, 1, 1) ./ dtCol;
    W = diff(zc_, 1, 1) ./ dtCol;

    Xm = 0.5 * (xc_(1:end-1,:) + xc_(2:end,:));
    Ym = 0.5 * (yc_(1:end-1,:) + yc_(2:end,:));
    Zm = 0.5 * (zc_(1:end-1,:) + zc_(2:end,:));

    mask = ~isnan(Xm) & ~isnan(Ym) & ~isnan(Zm) & ...
           ~isnan(U)  & ~isnan(V)  & ~isnan(W);

    xs = Xm(mask);  ys = Ym(mask);  zs = Zm(mask);
    us = U(mask);   vs = V(mask);   ws = W(mask);

    % Bin lookup
    ix = discretize(xs, gridX);
    iy = discretize(ys, gridY);
    iz = discretize(zs, gridZ);

    inDomain = ~isnan(ix) & ~isnan(iy) & ~isnan(iz);
    ix = ix(inDomain); iy = iy(inDomain); iz = iz(inDomain);
    us = us(inDomain); vs = vs(inDomain); ws = ws(inDomain);

    if isempty(ix), continue; end

    subs = [ix, iy, iz];
    counts = counts + accumarray(subs, 1,  [nx ny nz]);
    sumU   = sumU   + accumarray(subs, us, [nx ny nz]);
    sumV   = sumV   + accumarray(subs, vs, [nx ny nz]);
    sumW   = sumW   + accumarray(subs, ws, [nx ny nz]);

    fprintf('  chunk %d/%d done\n', c, nChunks);
end

Umean = sumU ./ counts;   % NaN where counts==0
Vmean = sumV ./ counts;
Wmean = sumW ./ counts;

fprintf('Total samples binned: %d\n', sum(counts(:)));


%Bin Reynolds Stresses 
sumUU = zeros(nx, ny, nz);
sumVV = zeros(nx, ny, nz);
sumWW = zeros(nx, ny, nz);
sumUV = zeros(nx, ny, nz);
sumUW = zeros(nx, ny, nz);
sumVW = zeros(nx, ny, nz);

fprintf('Reynolds stress pass: %d chunks...\n', nChunks);

for c = 1:nChunks
    i0 = (c-1)*chunkSize + 1;
    i1 = min(c*chunkSize, nTracks);

    xc_ = x(:, i0:i1);
    yc_ = y(:, i0:i1);
    zc_ = z(:, i0:i1);

    U = diff(xc_, 1, 1) ./ dtCol;
    V = diff(yc_, 1, 1) ./ dtCol;
    W = diff(zc_, 1, 1) ./ dtCol;

    Xm = 0.5 * (xc_(1:end-1,:) + xc_(2:end,:));
    Ym = 0.5 * (yc_(1:end-1,:) + yc_(2:end,:));
    Zm = 0.5 * (zc_(1:end-1,:) + zc_(2:end,:));

    mask = ~isnan(Xm) & ~isnan(Ym) & ~isnan(Zm) & ...
           ~isnan(U)  & ~isnan(V)  & ~isnan(W);

    xs = Xm(mask);  ys = Ym(mask);  zs = Zm(mask);
    us = U(mask);   vs = V(mask);   ws = W(mask);

    ix = discretize(xs, gridX);
    iy = discretize(ys, gridY);
    iz = discretize(zs, gridZ);

    inDomain = ~isnan(ix) & ~isnan(iy) & ~isnan(iz);
    ix = ix(inDomain); iy = iy(inDomain); iz = iz(inDomain);
    us = us(inDomain); vs = vs(inDomain); ws = ws(inDomain);

    if isempty(ix), continue; end

    % Look up mean at each sample's bin via linear indexing
    lin = sub2ind([nx ny nz], ix, iy, iz);
    up = us - Umean(lin);
    vp = vs - Vmean(lin);
    wp = ws - Wmean(lin);

    subs = [ix, iy, iz];
    sumUU = sumUU + accumarray(subs, up.*up, [nx ny nz]);
    sumVV = sumVV + accumarray(subs, vp.*vp, [nx ny nz]);
    sumWW = sumWW + accumarray(subs, wp.*wp, [nx ny nz]);
    sumUV = sumUV + accumarray(subs, up.*vp, [nx ny nz]);
    sumUW = sumUW + accumarray(subs, up.*wp, [nx ny nz]);
    sumVW = sumVW + accumarray(subs, vp.*wp, [nx ny nz]);

    fprintf('  chunk %d/%d done\n', c, nChunks);
end

% Normalize
uu = sumUU ./ counts;
vv = sumVV ./ counts;
ww = sumWW ./ counts;
uv = sumUV ./ counts;
uw = sumUW ./ counts;
vw = sumVW ./ counts;

% Turbulent kinetic energy
tke = 0.5 * (uu + vv + ww);

%Binned Dispersive fluxes 

% Step 1: plane-averaged means (1D profiles in y)
U_plane = squeeze(mean(Umean, [1 3], 'omitnan'));   % ny x 1
V_plane = squeeze(mean(Vmean, [1 3], 'omitnan'));
W_plane = squeeze(mean(Wmean, [1 3], 'omitnan'));

% Step 2: spatial deviations of the time-mean field
% Broadcast the 1D profile back to 3D shape and subtract.
U_tilde = Umean - reshape(U_plane, [1 ny 1]);
V_tilde = Vmean - reshape(V_plane, [1 ny 1]);
W_tilde = Wmean - reshape(W_plane, [1 ny 1]);

% Step 3: dispersive flux components at each y (plane average of products)
uu_disp = squeeze(mean(U_tilde .* U_tilde, [1 3], 'omitnan'));
vv_disp = squeeze(mean(V_tilde .* V_tilde, [1 3], 'omitnan'));
ww_disp = squeeze(mean(W_tilde .* W_tilde, [1 3], 'omitnan'));
uv_disp = squeeze(mean(U_tilde .* V_tilde, [1 3], 'omitnan'));
uw_disp = squeeze(mean(U_tilde .* W_tilde, [1 3], 'omitnan'));
vw_disp = squeeze(mean(V_tilde .* W_tilde, [1 3], 'omitnan'));


%Pressure Poisson

%Velocity correlation tensors
Tuu = Umean.*Umean + uu;
Tvv = Vmean.*Vmean + vv;
Tww = Wmean.*Wmean + ww;
Tuv = Umean.*Vmean + uv;
Tuw = Umean.*Wmean + uw;
Tvw = Vmean.*Wmean + vw;

%For any empty bins, use zero (near wall)
Tuu(isnan(Tuu)) = 0; Tvv(isnan(Tvv)) = 0; Tww(isnan(Tww)) = 0;
Tuv(isnan(Tuv)) = 0; Tuw(isnan(Tuw)) = 0; Tvw(isnan(Tvw)) = 0;

%Wavenumbers for (x,z) directions
Lx = gridX(end) - gridX(1);
Ly = gridY(end) - gridY(1);
Lz = gridZ(end) - gridZ(1);

kx_vec = 2*pi/Lx * [0:floor(nx/2), -floor(nx/2)+1:-1]';
ky_vec = 2*pi/Ly * [0:floor(ny/2), -floor(ny/2)+1:-1]';
kz_vec = 2*pi/Lz * [0:floor(nz/2), -floor(nz/2)+1:-1]';

% Source term
Tuu_h = fftn(Tuu);
Tvv_h = fftn(Tvv);
Tww_h = fftn(Tww);
Tuv_h = fftn(Tuv);
Tuw_h = fftn(Tuw);
Tvw_h = fftn(Tvw);

KX = reshape(kx_vec, [nx, 1,  1 ]);
KY = reshape(ky_vec, [1,  ny, 1 ]);
KZ = reshape(kz_vec, [1,  1,  nz]);

S_hat =   (KX.^2)       .* Tuu_h ...
        + (KY.^2)       .* Tvv_h ...
        + (KZ.^2)       .* Tww_h ...
        + 2*(KX.*KY)    .* Tuv_h ...
        + 2*(KX.*KZ)    .* Tuw_h ...
        + 2*(KY.*KZ)    .* Tvw_h;

K2 = KX.^2 + KY.^2 + KZ.^2;
K2(K2 < 1e-16) = Inf;

P_hat = -S_hat ./ K2;
P_hat_mean = 0;
P_hat(1,1) = P_hat_mean;

%Inv FFT
Pmean = real(ifftn(P_hat));

%Error analysis
%Interpolate DNS profiles onto the LPT y-grid for this bin resolution
Umean_dns_interp = interp1(yc_dns, Umean_dns, yc, 'pchip');
uu_dns_interp = interp1(yc_dns, uu_dns,    yc, 'pchip');
vv_dns_interp = interp1(yc_dns, vv_dns,    yc, 'pchip');
ww_dns_interp = interp1(yc_dns, ww_dns,    yc, 'pchip');
uv_dns_interp = interp1(yc_dns, uv_dns,    yc, 'pchip');
Pmean_dns_interp = interp1(yc_dns, Pmean_dns, yc, 'pchip');

% LPT profiles for this bin resolution
Uprofile_lpt = squeeze(mean(Umean, [1 3], 'omitnan'));
uu_lpt = squeeze(mean(uu, [1 3], 'omitnan'));
vv_lpt = squeeze(mean(vv, [1 3], 'omitnan'));
ww_lpt = squeeze(mean(ww, [1 3], 'omitnan'));
uv_lpt = squeeze(mean(uv, [1 3], 'omitnan'));
P_prof = squeeze(mean(Pmean, [1 3], 'omitnan'));

% Subtract centerline pressure
[~, ic_lpt] = min(abs(yc));
P_prof = P_prof - P_prof(ic_lpt);

% Force all profiles to column vectors
Uprofile_lpt = Uprofile_lpt(:);
uu_lpt = uu_lpt(:); vv_lpt = vv_lpt(:);
ww_lpt = ww_lpt(:); uv_lpt = uv_lpt(:);
P_prof = P_prof(:);
Umean_dns_interp = Umean_dns_interp(:);
uu_dns_interp = uu_dns_interp(:); vv_dns_interp = vv_dns_interp(:);
ww_dns_interp = ww_dns_interp(:); uv_dns_interp = uv_dns_interp(:);
Pmean_dns_interp = Pmean_dns_interp(:);

% Pointwise errors
err_U(bn, 1:ny)  = (Uprofile_lpt - Umean_dns_interp)';
err_uu(bn, 1:ny) = (uu_lpt - uu_dns_interp)';
err_vv(bn, 1:ny) = (vv_lpt - vv_dns_interp)';
err_ww(bn, 1:ny) = (ww_lpt - ww_dns_interp)';
err_uv(bn, 1:ny) = (uv_lpt - uv_dns_interp)';
err_P(bn, 1:ny)  = (P_prof  - Pmean_dns_interp)';

% L2 relative error
L2_rel = @(lpt, dns) norm(lpt - dns) / norm(dns);

L2_U(bn)  = L2_rel(Uprofile_lpt, Umean_dns_interp);
L2_uu(bn) = L2_rel(uu_lpt, uu_dns_interp);
L2_vv(bn) = L2_rel(vv_lpt, vv_dns_interp);
L2_ww(bn) = L2_rel(ww_lpt, ww_dns_interp);
L2_uv(bn) = L2_rel(uv_lpt, uv_dns_interp);
L2_P(bn)  = L2_rel(P_prof,  Pmean_dns_interp);

% Linf relative error
Linf_rel = @(lpt, dns) max(abs(lpt - dns)) / max(abs(dns));

Linf_U(bn)  = Linf_rel(Uprofile_lpt, Umean_dns_interp);
Linf_uu(bn) = Linf_rel(uu_lpt, uu_dns_interp);
Linf_vv(bn) = Linf_rel(vv_lpt, vv_dns_interp);
Linf_ww(bn) = Linf_rel(ww_lpt, ww_dns_interp);
Linf_uv(bn) = Linf_rel(uv_lpt, uv_dns_interp);
Linf_P(bn)  = Linf_rel(P_prof,  Pmean_dns_interp);

% Store bin info for later plotting
ny_list(bn) = ny;
dy_mean(bn) = mean(diff(gridY));   % mean bin width in y

fprintf('\n--- Bin resolution %d: nx=%d, ny=%d, nz=%d ---\n', bn, nx, ny, nz);
fprintf('  L2  errors: U=%.4f, uu=%.4f, vv=%.4f, ww=%.4f, uv=%.4f, P=%.4f\n', ...
        L2_U(bn), L2_uu(bn), L2_vv(bn), L2_ww(bn), L2_uv(bn), L2_P(bn));
fprintf('  Linf errors: U=%.4f, uu=%.4f, vv=%.4f, ww=%.4f, uv=%.4f, P=%.4f\n', ...
        Linf_U(bn), Linf_uu(bn), Linf_vv(bn), Linf_ww(bn), Linf_uv(bn), Linf_P(bn));

counts_y = squeeze(sum(counts, [1 3]));
min_samples_per_y(bn) = min(counts_y);
fprintf('  Min samples in any y-bin: %d\n', min_samples_per_y(bn));

end

%% -----------PLOTTING BELOW HERE------------------

%% Binning Visualization
xc = 0.5*(gridX(1:end-1) + gridX(2:end));
yc = 0.5*(gridY(1:end-1) + gridY(2:end));
zc = 0.5*(gridZ(1:end-1) + gridZ(2:end));

Uprofile = squeeze(mean(Umean, [1 3], 'omitnan'));
figure;
plot(yc, Uprofile);
xlabel('$y$'); ylabel('$\langle u \rangle$');
title('Mean streamwise velocity profile');
grid on;

[~, kmid] = min(abs(zc - mean(zc)));
figure;
imagesc(xc, yc, squeeze(Umean(:,:,kmid))');
axis xy equal tight; 
c = colorbar;
c.Label.String = 'u';
c.Label.FontName = 'Latex';
xlabel('$x$'); ylabel('$y$');
title(sprintf('$\\langle u \\rangle$ at $z = %.2f$', zc(kmid)));

figure;
imagesc(xc, yc, squeeze(sum(counts,3))');
axis xy equal tight; 
c = colorbar;
c.Label.String = 'Samples';
c.Label.FontName = 'Latex';
xlabel('$x$'); ylabel('$y$');
title('Samples per $(x,y)$ bin (summed over $z$)');

%% Binned Reynolds stress plots
uu_prof = squeeze(mean(uu, [1 3], 'omitnan'));
vv_prof = squeeze(mean(vv, [1 3], 'omitnan'));
ww_prof = squeeze(mean(ww, [1 3], 'omitnan'));
uv_prof = squeeze(mean(uv, [1 3], 'omitnan'));
tke_prof = squeeze(mean(tke, [1 3], 'omitnan'));

figure;
hold on;
plot(yc, uu_prof);
plot(yc, vv_prof);
plot(yc, ww_prof);
plot(yc, uv_prof);
xlabel('$y$'); ylabel('Reynolds stress');
legend({'$\overline{u''u''}$','$\overline{v''v''}$','$\overline{w''w''}$','$\overline{u''v''}$'}, ...
       'Location','best');
title('Reynolds stress profiles');
grid on;
hold off;

figure;
plot(yc, tke_prof);
xlabel('$y$'); ylabel('$k = \frac{1}{2}\overline{u_i'' u_i''}$');
title('Turbulent kinetic energy profile');
grid on;

%% Dispersive fluxes vs Reynolds stresses plots
figure;
tiledlayout(2, 2);

nexttile;
plot(yc, uu_prof, 'LineWidth', 2); hold on;
plot(yc, uu_disp, 'LineWidth', 2);
xlabel('$y$'); ylabel('$\overline{u''u''}$ vs $\langle \tilde u \tilde u \rangle$');
legend({'Reynolds', 'Dispersive'}, 'Location','best');
title('$uu$'); grid on;

nexttile;
plot(yc, vv_prof, 'LineWidth', 2); hold on;
plot(yc, vv_disp, 'LineWidth', 2);
xlabel('$y$'); ylabel('$\overline{v''v''}$ vs $\langle \tilde v \tilde v \rangle$');
title('$vv$'); grid on;

nexttile;
plot(yc, ww_prof, 'LineWidth', 2); hold on;
plot(yc, ww_disp, 'LineWidth', 2);
xlabel('$y$'); ylabel('$\overline{w''w''}$ vs $\langle \tilde w \tilde w \rangle$');
title('$ww$'); grid on;

nexttile;
plot(yc, uv_prof, 'LineWidth', 2); hold on;
plot(yc, uv_disp, 'LineWidth', 2);
xlabel('$y$'); ylabel('$\overline{u''v''}$ vs $\langle \tilde u \tilde v \rangle$');
title('$uv$'); grid on;

figure;
plot(yc, uu_disp, 'DisplayName', '$|\langle\tilde u \tilde u\rangle|$'); hold on;
plot(yc, vv_disp, 'DisplayName', '$|\langle\tilde v \tilde v\rangle|$');
plot(yc, ww_disp, 'DisplayName', '$|\langle\tilde w \tilde w\rangle|$');
plot(yc, uv_disp, 'DisplayName', '$|\langle\tilde u \tilde v\rangle|$');
xlabel('$y$'); ylabel('Dispersive flux (abs)');
title('Dispersive flux magnitudes (log scale)');
legend('Interpreter','latex','Location','best'); grid on;

%% Pressure Plots
%Plane averaged pressure
P_prof = squeeze(mean(Pmean, [1 3], 'omitnan'));

figure;
plot(yc, P_prof, 'LineWidth', 2); hold on;
plot(yc, -vv_prof, 'LineWidth', 2);
plot(yc, P_prof + vv_prof, 'k--', 'LineWidth', 2);
xlabel('$y$');
legend({'$\bar{p}(y)$', '$-\overline{v''v''}(y)$', ...
        '$\bar{p} + \overline{v''v''}$ (should be const)'}, ...
       'Location', 'best');
title('Mean pressure profile from Poisson solve');
grid on;

% Spurious (x,z) variation — should be near zero for channel flow
P_std_xz = squeeze(std(Pmean, 0, [1 3]));
figure;
semilogy(yc, P_std_xz, 'LineWidth', 2);
xlabel('$y$'); ylabel('std of $\bar{p}$ over $(x,z)$');
title('Spurious $(x,z)$ variation (should be near zero)');
grid on;

% Mid-z slice of pressure field
[~, kmid] = min(abs(zc - mean(zc)));
figure;
imagesc(xc, yc, squeeze(Pmean(:,:,kmid))');
axis xy equal tight;
c = colorbar;
c.Label.String = '$\bar{p}$';
c.Label.Interpreter = 'latex';
xlabel('$x$'); ylabel('$y$');
title(sprintf('$\\bar{p}(x,y)$ at $z \\approx %.2f$', zc(kmid)));

%% Velocity comparison plots
Uprofile_lpt = squeeze(mean(Umean, [1 3], 'omitnan'));
uu_lpt = squeeze(mean(uu, [1 3], 'omitnan'));
vv_lpt = squeeze(mean(vv, [1 3], 'omitnan'));
ww_lpt = squeeze(mean(ww, [1 3], 'omitnan'));
uv_lpt = squeeze(mean(uv, [1 3], 'omitnan'));
tke_lpt = 0.5 * (uu_lpt + vv_lpt + ww_lpt);

figure;
plot(yc, Uprofile_lpt, 'b-', 'LineWidth', 2); hold on;
plot(yc_dns, Umean_dns, 'r--', 'LineWidth', 2);
xlabel('$y$'); ylabel('$\langle u \rangle$');
legend({'LPT', 'DNS'}, 'Location', 'best');
title('Mean streamwise velocity');
grid on;

figure;
tiledlayout(2, 2);

nexttile;
plot(yc, uu_lpt, 'b-', 'LineWidth', 2); hold on;
plot(yc_dns, uu_dns, 'r--', 'LineWidth', 2);
xlabel('$y$'); title('$\overline{u''u''}$');
legend({'LPT', 'DNS'}, 'Location', 'best'); grid on;

nexttile;
plot(yc, vv_lpt, 'b-', 'LineWidth', 2); hold on;
plot(yc_dns, vv_dns, 'r--', 'LineWidth', 2);
xlabel('$y$'); title('$\overline{v''v''}$');
grid on;

nexttile;
plot(yc, ww_lpt, 'b-', 'LineWidth', 2); hold on;
plot(yc_dns, ww_dns, 'r--', 'LineWidth', 2);
xlabel('$y$'); title('$\overline{w''w''}$');
grid on;

nexttile;
plot(yc, uv_lpt, 'b-', 'LineWidth', 2); hold on;
plot(yc_dns, uv_dns, 'r--', 'LineWidth', 2);
xlabel('$y$'); title('$\overline{u''v''}$');
grid on;

figure;
plot(yc, tke_lpt, 'b-', 'LineWidth', 2); hold on;
plot(yc_dns, tke_dns, 'r--', 'LineWidth', 2);
xlabel('$y$'); ylabel('$k$');
legend({'LPT', 'DNS'}, 'Location', 'best');
title('Turbulent kinetic energy');
grid on;



%% Compare DNS pressure to LPT Poisson pressure
figure;
plot(yc, P_prof, 'b-', 'LineWidth', 2); hold on;
plot(yc_dns, Pmean_dns, 'r--', 'LineWidth', 2);
xlabel('$y$'); ylabel('$\bar{p}(y)$');
legend({'LPT (Poisson solve)', 'DNS'}, 'Location', 'best');
title('Mean pressure profile');
grid on;

figure;
plot(yc_dns, Pmean_dns + vv_dns, 'r--', 'LineWidth', 2); hold on;
plot(yc, P_prof + vv_lpt, 'b-', 'LineWidth', 2);
xlabel('$y$');
legend({'DNS: $\bar{p} + \overline{v''v''}$', 'LPT: $\bar{p} + \overline{v''v''}$'}, ...
       'Location', 'best');
title('$\bar{p} + \overline{v''v''}$');
grid on;

%% Error plotting

if length(Ybin_vec) > 1

    figure;
    semilogy(ny_list, L2_U,  'o-', 'LineWidth', 2, 'DisplayName', '$\langle u \rangle$'); hold on;
    semilogy(ny_list, L2_uu, 's-', 'LineWidth', 2, 'DisplayName', '$\overline{u''u''}$');
    semilogy(ny_list, L2_vv, 'd-', 'LineWidth', 2, 'DisplayName', '$\overline{v''v''}$');
    semilogy(ny_list, L2_ww, '^-', 'LineWidth', 2, 'DisplayName', '$\overline{w''w''}$');
    semilogy(ny_list, L2_uv, 'v-', 'LineWidth', 2, 'DisplayName', '$\overline{u''v''}$');
    semilogy(ny_list, L2_P,  'p-', 'LineWidth', 2, 'DisplayName', '$\bar{p}$');
    xlabel('Number of $y$-bins'); ylabel('$L_2$ relative error');
    legend('Location', 'best');
    title('Convergence: $L_2$ error vs $y$-resolution');
    grid on;

    figure;
    semilogy(ny_list, Linf_U,  'o-', 'LineWidth', 2, 'DisplayName', '$\langle u \rangle$'); hold on;
    semilogy(ny_list, Linf_uu, 's-', 'LineWidth', 2, 'DisplayName', '$\overline{u''u''}$');
    semilogy(ny_list, Linf_vv, 'd-', 'LineWidth', 2, 'DisplayName', '$\overline{v''v''}$');
    semilogy(ny_list, Linf_ww, '^-', 'LineWidth', 2, 'DisplayName', '$\overline{w''w''}$');
    semilogy(ny_list, Linf_uv, 'v-', 'LineWidth', 2, 'DisplayName', '$\overline{u''v''}$');
    semilogy(ny_list, Linf_P,  'p-', 'LineWidth', 2, 'DisplayName', '$\bar{p}$');
    xlabel('Number of $y$-bins'); ylabel('$L_\infty$ relative error');
    legend('Location', 'best');
    title('Convergence: $L_\infty$ error vs $y$-resolution');
    grid on;

    % Error vs mean bin width (more physical than bin count)
    figure;
    loglog(dy_mean, L2_U,  'o-', 'LineWidth', 2, 'DisplayName', '$\langle u \rangle$'); hold on;
    loglog(dy_mean, L2_uu, 's-', 'LineWidth', 2, 'DisplayName', '$\overline{u''u''}$');
    loglog(dy_mean, L2_vv, 'd-', 'LineWidth', 2, 'DisplayName', '$\overline{v''v''}$');
    loglog(dy_mean, L2_ww, '^-', 'LineWidth', 2, 'DisplayName', '$\overline{w''w''}$');
    loglog(dy_mean, L2_uv, 'v-', 'LineWidth', 2, 'DisplayName', '$\overline{u''v''}$');
    loglog(dy_mean, L2_P,  'p-', 'LineWidth', 2, 'DisplayName', '$\bar{p}$');
    xlabel('Mean $\Delta y$'); ylabel('$L_2$ relative error');
    legend('Location', 'best');
    title('Convergence: $L_2$ error vs bin width');
    grid on;
    set(gca, 'XDir', 'reverse');   % smaller bins to the right

end


