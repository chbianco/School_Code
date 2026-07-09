%% Preamble
close all; clc;
clearvars 
set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);

%% Loading file
file = '07_06_26_freestream_tracks'; %Tracks to analyze.
% This file stores tracks in "long" (row-per-observation) format: a single
% matrix variable `data`, [nObservations x 19], with columns:
%   1: track ID (arbitrary integer label, not contiguous)
%   2: frame number (integer, increases by 1 within a track, no gaps)
%   3-5: x, y, z position
%   6-8: u, v, w velocity
%   9-19: additional columns (likely accelerations + per-camera image
%         coordinates based on typical OpenLPT export layout)

% Time vector. The file only stores integer frame numbers, not a
    % physical time/frame-rate, so dt_frame defaults to 1 (t is in units
    % of "frames"). Set dt_frame to your camera's actual sample interval
    % (e.g. 1/fps, in seconds) if you want t in physical time units
    dt_frame = 1/400;

%% Load tracks

    raw = load(file);
    data = raw.data;          % [trackID, frame, x, y, z, u, v, w, ...]

    trackID = data(:,1);
    frame   = data(:,2);
    xPos    = data(:,3);
    yPos    = data(:,4);
    zPos    = data(:,5);

    uVel = data(:,6);
    vVel = data(:,7);
    wVel = data(:,8);

    % Map arbitrary track IDs -> contiguous column indices 1..nTracks
    [~, ~, colIdx] = unique(trackID);
    nTracks = max(colIdx);

    % Build a shared frame axis spanning every sample; each sample's row
    % position is its frame number offset from the global minimum frame.
    frameMin = min(frame);
    frameMax = max(frame);
    nT = frameMax - frameMin + 1;
    rowIdx = frame - frameMin + 1;

    t = (0:nT-1)' * dt_frame;

    % Scatter long-format samples into NaN-padded wide matrices so the
    % rest of the script (written for nT x nTracks position matrices)
    % works unchanged.
    x = NaN(nT, nTracks);
    y = NaN(nT, nTracks);
    z = NaN(nT, nTracks);

    u = NaN(nT, nTracks);
    v = NaN(nT, nTracks);
    w = NaN(nT, nTracks);


    linIdx = sub2ind([nT, nTracks], rowIdx, colIdx);
    x(linIdx) = xPos;
    y(linIdx) = yPos;
    z(linIdx) = zPos;

    u(linIdx) = uVel;
    v(linIdx) = vVel;
    w(linIdx) = wVel;

    tracks.t = t; tracks.x = x; tracks.y = y; tracks.z = z;
    tracks.u = u; tracks.v = v; tracks.w = w;

    fprintf('Loaded: %d samples x %d tracks (%.6f GB per array)\n', ...
        nT, nTracks, nT*nTracks*8/1e9);

%% User inputs
%Number of tracks to analyze. Set to nTracks to keep all tracks
N = nTracks;  

%Number of trajectories to plot (clamped to nTracks)
n = min(nTracks, 200);

%Bounds of LPT view area, as a 2x1 (ie 0, 20)
% Auto-computed from the data with 5% padding
% Override manually below if you want a fixed
% window instead (e.g. for comparing multiple datasets on the same axes).
xLo = min(x(:)); xHi = max(x(:));
yLo = min(y(:)); yHi = max(y(:));
zLo = min(z(:)); zHi = max(z(:));
padFrac = 0.05;
Xlim = [xLo - padFrac*(xHi-xLo), xHi + padFrac*(xHi-xLo)];
Ylim = [yLo - padFrac*(yHi-yLo), yHi + padFrac*(yHi-yLo)];
Zlim = [zLo - padFrac*(zHi-zLo), zHi + padFrac*(zHi-zLo)];

%Number of bins in x, y, and z. Evenly spaced
Xbin = 5;
Ybin = 5;
Zbin = 5;

%Tracks per chunk; tune for memory. 1000 works well. Only need for a bunch
%of tracks
chunkSize = 1000;  

%Max lags for Lagrangian statistics
%Clamped to (longest available track - 1)
max_lag = min(2500, max(sum(~isnan(x), 1)) - 1);

%% Sort positions and velocities by track length 
x_lng = sum(~isnan(x), 1);
[~, sorted_idx] = sort(x_lng, 'descend');

x = x(:, sorted_idx);
y = y(:, sorted_idx);
z = z(:, sorted_idx);

u = u(:, sorted_idx);
v = v(:, sorted_idx);
w = w(:, sorted_idx);


%% Plot n tracks colored by speed
maxTrackLen = max(sum(~isnan(x), 1)); %sized to the dataset
speed = NaN(maxTrackLen, n);

figure(1)
hold on

for k = 1:n
    xn = rmmissing(x(:,k));
    yn = rmmissing(y(:,k));
    zn = rmmissing(z(:,k));

    trk_lng = min([length(xn),length(yn), length(zn)]);
    un = zeros(1, trk_lng -1);
    vn = zeros(1, trk_lng -1);
    wn = zeros(1, trk_lng -1);
    t_trk = t(1:trk_lng)';
    
    %Velocities
    un = rmmissing(u(:,k));
    vn = rmmissing(v(:,k));
    wn = rmmissing(w(:,k));
    
    speed(1:trk_lng, k) = sqrt(un'.^2 + vn'.^2 + wn'.^2)';

    scatter3(xn(1:end-1),yn(1:end-1),zn(1:end-1), 20, speed(1:trk_lng-1,k), 'filled')


end 

c = colorbar;
c.Label.String = 'Speed';
c.Label.FontName = 'Latex';
xlabel('x')
ylabel('y')
zlabel('z')
xlim(Xlim)
ylim(Ylim)
zlim(Zlim)
axis equal

view(3)

hold off

%% Keep N longest tracks and organize
trackLengths = sum(~isnan(x), 1);
keep = 1:min(N, nTracks);

x = x(:, keep);
y = y(:, keep);
z = z(:, keep);

u = u(:, keep);
v = v(:, keep);
w = w(:, keep);


% Trim rows to longest surviving track
maxLen = max(sum(~isnan(x), 1));
x = x(1:maxLen, :);
y = y(1:maxLen, :);
z = z(1:maxLen, :);

u = u(1:maxLen, :);
v = v(1:maxLen, :);
w = w(1:maxLen, :);

t = t(1:maxLen);

[nT, nTracks] = size(x);
fprintf('Kept %d tracks, max length %d samples\n', nTracks, nT);

%% Eulerian Binning
%Define bin grid
gridX = linspace(Xlim(1), Xlim(2), Xbin + 1);
gridY = linspace(Ylim(1), Ylim(2), Ybin + 1);
gridZ = linspace(Zlim(1), Zlim(2), Zbin + 1);

%Number of bins
nx = numel(gridX) - 1;
ny = numel(gridY) - 1;
nz = numel(gridZ) - 1;

%Centers
xc = 0.5*(gridX(1:end-1) + gridX(2:end));
yc = 0.5*(gridY(1:end-1) + gridY(2:end));
zc = 0.5*(gridZ(1:end-1) + gridZ(2:end));

%Chunked binning: accumulate sums over track chunks
counts = zeros(nx, ny, nz);
sumU   = zeros(nx, ny, nz);
sumV   = zeros(nx, ny, nz);
sumW   = zeros(nx, ny, nz);

dtCol = diff(t(:));
nChunks = ceil(nTracks / chunkSize);
fprintf('Processing %d chunks of up to %d tracks...\n', nChunks, chunkSize);

for c = 1:nChunks
    i0 = (c-1)*chunkSize + 1;
    i1 = min(c*chunkSize, nTracks)-1;

    xc_ = x(:, i0:i1);
    yc_ = y(:, i0:i1);
    zc_ = z(:, i0:i1);

    U = u(1:end-1, i0:i1);
    V = v(1:end-1, i0:i1);
    W = w(1:end-1, i0:i1);

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

%% Reynolds stresses
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

    U = u(1:end-1, i0:i1);
    V = v(1:end-1, i0:i1);
    W = w(1:end-1, i0:i1);

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

%% Dispersive fluxes
% Step 1: plane-averaged means (1D profiles in y)
U_plane = squeeze(mean(Umean, [1 3], 'omitnan'));   % ny x 1
V_plane = squeeze(mean(Vmean, [1 3], 'omitnan'));
W_plane = squeeze(mean(Wmean, [1 3], 'omitnan'));

% Step 2: spatial deviations of the time-mean field
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

%% Lagrangian Statistics
% Accumulators
Ruu_sum = zeros(max_lag + 1, 1);
Rvv_sum = zeros(max_lag + 1, 1);
Rww_sum = zeros(max_lag + 1, 1);
R_count = zeros(max_lag + 1, 1);   % number of valid pairs per lag

dt_track = t(2) - t(1);   % assuming uniform time step
tau_vec = (0:max_lag)' * dt_track;

% U_plane/V_plane/W_plane can contain NaN at y-bins with zero samples
% (common with sparse data + a fine grid). interp1 with 'pchip' will
% silently propagate those NaNs across the whole interpolated profile,
% so the empty bins are dropped here before building the interpolant.
validProfile = ~isnan(U_plane) & ~isnan(V_plane) & ~isnan(W_plane);
if nnz(validProfile) < 2
    error(['Fewer than 2 non-empty y-bins in U_plane/V_plane/W_plane - ', ...
           'reduce Ybin or widen Ylim so the Eulerian grid actually ', ...
           'captures samples before running the Lagrangian section.'])
end
yc_valid = yc(validProfile);
U_plane_valid = U_plane(validProfile);
V_plane_valid = V_plane(validProfile);
W_plane_valid = W_plane(validProfile);

fprintf('Computing Lagrangian autocorrelation across %d tracks...\n', nTracks);

for k = 1:nTracks
    % Extract this track's positions, remove NaNs
    xk = x(:, k);
    yk = y(:, k);
    zk = z(:, k);

    valid = ~isnan(xk) & ~isnan(yk) & ~isnan(zk);
    idx_valid = find(valid);

    if length(idx_valid) < 10
        continue   % skip very short tracks
    end

    % Contiguous segment (tracks might have gaps — take longest contiguous block)
    % Assume tracks are contiguous from start until first NaN
    i_start = idx_valid(1);
    i_end = idx_valid(1);
    for ii = 2:length(idx_valid)
        if idx_valid(ii) == idx_valid(ii-1) + 1
            i_end = idx_valid(ii);
        else
            break
        end
    end

    seg_len = i_end - i_start + 1;
    if seg_len < 10
        continue
    end

    % Velocity via central differences. CHANGE FOR OPENLPT DATA
    xseg = xk(i_start:i_end);
    yseg = yk(i_start:i_end);
    zseg = zk(i_start:i_end);

   % Velocity via central differences
    u_trk = gradient(xseg, dt_track);
    v_trk = gradient(yseg, dt_track);
    w_trk = gradient(zseg, dt_track);

    % Subtract LOCAL Eulerian mean at the particle's y-position
    % Interpolate the plane-averaged mean profile to each track point's y
    % U_plane, V_plane, W_plane are the (ny x 1) plane-averaged profiles
    % computed in the binning section. yc is the bin-center vector.
    u_mean_local = interp1(yc_valid, U_plane_valid, yseg, 'pchip', 'extrap');
    v_mean_local = interp1(yc_valid, V_plane_valid, yseg, 'pchip', 'extrap');
    w_mean_local = interp1(yc_valid, W_plane_valid, yseg, 'pchip', 'extrap');

    u_trk = u_trk - u_mean_local;
    v_trk = v_trk - v_mean_local;
    w_trk = w_trk - w_mean_local;

    % Compute autocorrelation for this track at each lag
    n = length(u_trk);
    max_lag_this = min(max_lag, n - 1);

    for lag = 0:max_lag_this
        n_pairs = n - lag;

        uu_corr = sum(u_trk(1:n_pairs) .* u_trk(1+lag:n_pairs+lag));
        vv_corr = sum(v_trk(1:n_pairs) .* v_trk(1+lag:n_pairs+lag));
        ww_corr = sum(w_trk(1:n_pairs) .* w_trk(1+lag:n_pairs+lag));

        Ruu_sum(lag+1) = Ruu_sum(lag+1) + uu_corr;
        Rvv_sum(lag+1) = Rvv_sum(lag+1) + vv_corr;
        Rww_sum(lag+1) = Rww_sum(lag+1) + ww_corr;
        R_count(lag+1) = R_count(lag+1) + n_pairs;
    end

    if mod(k, 10000) == 0
        fprintf('  %d/%d tracks processed\n', k, nTracks);
    end
end

% Normalize
Ruu_raw = Ruu_sum ./ R_count;
Rvv_raw = Rvv_sum ./ R_count;
Rww_raw = Rww_sum ./ R_count;

% Normalize to correlation coefficient (divide by zero-lag value)
Ruu = Ruu_raw / Ruu_raw(1);
Rvv = Rvv_raw / Rvv_raw(1);
Rww = Rww_raw / Rww_raw(1);

% Lagrangian integral timescales (integrate R from 0 to first zero crossing)
% Using trapezoidal rule up to the first time the tracks decays to 1/e 
% (or max_lag if no crossing). Assuming exponential decay 
i_zero_u = find(Ruu < exp(-1), 1, 'first');
i_zero_v = find(Rvv < exp(-1), 1, 'first');
i_zero_w = find(Rww < exp(-1), 1, 'first');

if isempty(i_zero_u), i_zero_u = max_lag + 1; end
if isempty(i_zero_v), i_zero_v = max_lag + 1; end
if isempty(i_zero_w), i_zero_w = max_lag + 1; end

TL_u = trapz(tau_vec(1:i_zero_u), Ruu(1:i_zero_u));
TL_v = trapz(tau_vec(1:i_zero_v), Rvv(1:i_zero_v));
TL_w = trapz(tau_vec(1:i_zero_w), Rww(1:i_zero_w));

% Eddy diffusivities
D_u = Ruu_raw(1) * TL_u;
D_v = Rvv_raw(1) * TL_v;
D_w = Rww_raw(1) * TL_w;

%% ----------------PLOTTING-------------------------------------------
%% Mean velocities
xc = 0.5*(gridX(1:end-1) + gridX(2:end));
yc = 0.5*(gridY(1:end-1) + gridY(2:end));
zc = 0.5*(gridZ(1:end-1) + gridZ(2:end));

%Mean streamwise velocity in y
Uprofile_y = squeeze(mean(Umean, [1 3], 'omitnan'));
figure;
plot(yc, Uprofile_y);
xlabel('$y$'); ylabel('$\langle \overline{u} \rangle$');
title('Mean streamwise velocity profile');
grid on;

%Mean streamwise velocity in z
Uprofile_z = squeeze(mean(Umean, [1 2], 'omitnan'));
figure;
plot(zc, Uprofile_z);
xlabel('$z$'); ylabel('$\langle \overline{u} \rangle$');
title('Mean streamwise velocity profile');
grid on;

%% Binning visualization
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

figure;
imagesc(xc, zc, squeeze(sum(counts,2))');
axis xy tight; 
c = colorbar;
c.Label.String = 'Samples';
c.Label.FontName = 'Latex';
xlabel('$x$'); ylabel('$z$');
title('Samples per $(x,z)$ bin (summed over $y$)');

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

%% Plot autocorrelations
figure;
plot(tau_vec, Ruu, 'LineWidth', 2); hold on;
plot(tau_vec, Rvv, 'LineWidth', 2);
plot(tau_vec, Rww, 'LineWidth', 2);
yline(0, 'k--');
xlabel('$\tau$'); ylabel('$R_{ii}(\tau)$');
legend({'$R_{uu}$', '$R_{vv}$', '$R_{ww}$'}, 'Location', 'best');
title('Lagrangian velocity autocorrelation');
grid on;
