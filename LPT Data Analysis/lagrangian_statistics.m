%% Preamble
close all; clc;
clearvars -except tracks t x y z
set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);

%% Load tracks
if exist('tracks', 'var') 
    t = tracks.t;
    x = tracks.x;
    y = tracks.y;
    z = tracks.z;
    [nT, nTracks] = size(x);
    fprintf('Track data already loaded')
else
    tracks = load('jhtdb_long_LPT.mat').LPT;
    t = tracks.t;
    x = tracks.x;
    y = tracks.y;
    z = tracks.z;
    [nT, nTracks] = size(x);
    fprintf('Loaded: %d samples x %d tracks (%.2f GB per array)\n', ...
        nT, nTracks, nT*nTracks*8/1e9);
end

%% Keep N longest tracks 
N = nTracks; %Set to nTracks to keep all tracks 

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

%% Binning to get Eularian mean velocities 
Xbin_vec = [128];
Ybin_vec = [32];
Zbin_vec = [64]; 

bn = 1;

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

Uprofile_lpt = squeeze(mean(Umean, [1 3], 'omitnan'));
Vprofile_lpt = squeeze(mean(Vmean, [1 3], 'omitnan'));
Wprofile_lpt = squeeze(mean(Wmean, [1 3], 'omitnan'));

fprintf('Total samples binned: %d\n', sum(counts(:)));


%% -----------LAGRANGIAN STATISTICS---------------
%% Calculate lagrangian velocity autocorrelations averaged in space
% Compute R_uu(tau), R_vv(tau), R_ww(tau) from individual track time series.

%Max lag in time steps 
max_lag = 2500;

% Accumulators
Ruu_sum = zeros(max_lag + 1, 1);
Rvv_sum = zeros(max_lag + 1, 1);
Rww_sum = zeros(max_lag + 1, 1);
R_count = zeros(max_lag + 1, 1);   % number of valid pairs per lag

dt_track = t(2) - t(1);   % assuming uniform time step
tau_vec = (0:max_lag)' * dt_track;

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

    % Velocity via central differences
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
    u_mean_local = interp1(yc, Uprofile_lpt, yseg, 'pchip', 'extrap');
    v_mean_local = interp1(yc, Vprofile_lpt, yseg, 'pchip', 'extrap');
    w_mean_local = interp1(yc, Wprofile_lpt, yseg, 'pchip', 'extrap');

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
% Using trapezoidal rule up to the first zero crossing (or max_lag if no crossing)
i_zero_u = find(Ruu < 0, 1, 'first');
i_zero_v = find(Rvv < 0, 1, 'first');
i_zero_w = find(Rww < 0, 1, 'first');

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

fprintf('\nLagrangian integral timescales:\n');
fprintf('  T_L(u) = %.4f\n', TL_u);
fprintf('  T_L(v) = %.4f\n', TL_v);
fprintf('  T_L(w) = %.4f\n', TL_w);
fprintf('Eddy diffusivities:\n');
fprintf('  D_u = %.6f\n', D_u);
fprintf('  D_v = %.6f\n', D_v);
fprintf('  D_w = %.6f\n', D_w);

%% Calculate autocorrelations binned in space (based on initial y)
% Requires from above: x, y, z, t, nTracks, yc, Uprofile_lpt, 
%                      Vprofile_lpt, Wprofile_lpt

dt_track = t(2) - t(1);
max_lag = 3000;
tau_vec = (0:max_lag)' * dt_track;

% Define y-bands for conditioning
y_band_centers = [-0.95, -0.8, -0.5, 0.0];
y_band_halfwidth = 0.05;
n_bands = length(y_band_centers);

% Accumulators: one set per y-band, plus one for "all tracks"
Ruu_sum = zeros(max_lag + 1, n_bands + 1);
Rvv_sum = zeros(max_lag + 1, n_bands + 1);
Rww_sum = zeros(max_lag + 1, n_bands + 1);
R_count = zeros(max_lag + 1, n_bands + 1);

fprintf('\n=== Lagrangian autocorrelation ===\n');
fprintf('max_lag = %d (tau_max = %.2f)\n', max_lag, tau_vec(end));
fprintf('y-bands: ');
for b = 1:n_bands
    fprintf('[%.2f ± %.2f]  ', y_band_centers(b), y_band_halfwidth);
end
fprintf('+ [all]\n');

for k = 1:nTracks
    xk = x(:, k);
    yk = y(:, k);
    zk = z(:, k);

    valid = ~isnan(xk) & ~isnan(yk) & ~isnan(zk);
    idx_valid = find(valid);
    if length(idx_valid) < 10, continue; end

    % Find longest contiguous segment from start
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
    if seg_len < 10, continue; end

    % Extract segment
    xseg = xk(i_start:i_end);
    yseg = yk(i_start:i_end);
    zseg = zk(i_start:i_end);

    % Velocity via central differences
    u_trk = gradient(xseg, dt_track);
    v_trk = gradient(yseg, dt_track);
    w_trk = gradient(zseg, dt_track);

    % Subtract local Eulerian mean at particle's y-position
    u_trk = u_trk - interp1(yc, Uprofile_lpt, yseg, 'pchip', 'extrap');
    v_trk = v_trk - interp1(yc, Vprofile_lpt, yseg, 'pchip', 'extrap');
    w_trk = w_trk - interp1(yc, Wprofile_lpt, yseg, 'pchip', 'extrap');

    % Determine which y-band this track starts in
    y0 = yseg(1);
    band_idx = 0;
    for b = 1:n_bands
        if abs(y0 - y_band_centers(b)) <= y_band_halfwidth
            band_idx = b;
            break
        end
    end

    % Compute autocorrelation
    n = length(u_trk);
    max_lag_this = min(max_lag, n - 1);

    for lag = 0:max_lag_this
        n_pairs = n - lag;
        uu_corr = sum(u_trk(1:n_pairs) .* u_trk(1+lag:n_pairs+lag));
        vv_corr = sum(v_trk(1:n_pairs) .* v_trk(1+lag:n_pairs+lag));
        ww_corr = sum(w_trk(1:n_pairs) .* w_trk(1+lag:n_pairs+lag));

        % Accumulate into "all tracks" column (last column)
        Ruu_sum(lag+1, end) = Ruu_sum(lag+1, end) + uu_corr;
        Rvv_sum(lag+1, end) = Rvv_sum(lag+1, end) + vv_corr;
        Rww_sum(lag+1, end) = Rww_sum(lag+1, end) + ww_corr;
        R_count(lag+1, end) = R_count(lag+1, end) + n_pairs;

        % Also accumulate into specific y-band if applicable
        if band_idx > 0
            Ruu_sum(lag+1, band_idx) = Ruu_sum(lag+1, band_idx) + uu_corr;
            Rvv_sum(lag+1, band_idx) = Rvv_sum(lag+1, band_idx) + vv_corr;
            Rww_sum(lag+1, band_idx) = Rww_sum(lag+1, band_idx) + ww_corr;
            R_count(lag+1, band_idx) = R_count(lag+1, band_idx) + n_pairs;
        end
    end

    if mod(k, 10000) == 0
        fprintf('  %d/%d tracks processed\n', k, nTracks);
    end
end

%Normalize and compute integral timescales
n_cols = n_bands + 1;   % bands + "all"
col_labels = cell(1, n_cols);
for b = 1:n_bands
    col_labels{b} = sprintf('$y_0 = %.2f$', y_band_centers(b));
end
col_labels{end} = 'All tracks';

Ruu = zeros(max_lag + 1, n_cols);
Rvv = zeros(max_lag + 1, n_cols);
Rww = zeros(max_lag + 1, n_cols);
TL_u = zeros(n_cols, 1);
TL_v = zeros(n_cols, 1);
TL_w = zeros(n_cols, 1);
D_u = zeros(n_cols, 1);
D_v = zeros(n_cols, 1);
D_w = zeros(n_cols, 1);

fprintf('\n--- Results ---\n');

for b = 1:n_cols
    if R_count(1, b) == 0
        fprintf('%s: no tracks found, skipping\n', col_labels{b});
        continue
    end

    Ruu_raw = Ruu_sum(:, b) ./ R_count(:, b);
    Rvv_raw = Rvv_sum(:, b) ./ R_count(:, b);
    Rww_raw = Rww_sum(:, b) ./ R_count(:, b);

    Ruu(:, b) = Ruu_raw / Ruu_raw(1);
    Rvv(:, b) = Rvv_raw / Rvv_raw(1);
    Rww(:, b) = Rww_raw / Rww_raw(1);

    % Integral timescales (to first zero crossing or end)
    iz_u = find(Ruu(:,b) < 0, 1, 'first');
    iz_v = find(Rvv(:,b) < 0, 1, 'first');
    iz_w = find(Rww(:,b) < 0, 1, 'first');
    if isempty(iz_u), iz_u = max_lag + 1; end
    if isempty(iz_v), iz_v = max_lag + 1; end
    if isempty(iz_w), iz_w = max_lag + 1; end

    TL_u(b) = trapz(tau_vec(1:iz_u), Ruu(1:iz_u, b));
    TL_v(b) = trapz(tau_vec(1:iz_v), Rvv(1:iz_v, b));
    TL_w(b) = trapz(tau_vec(1:iz_w), Rww(1:iz_w, b));

    D_u(b) = Ruu_raw(1) * TL_u(b);
    D_v(b) = Rvv_raw(1) * TL_v(b);
    D_w(b) = Rww_raw(1) * TL_w(b);

    fprintf('%s: n=%d, T_L(u)=%.3f, T_L(v)=%.3f, T_L(w)=%.3f, D_u=%.5f, D_v=%.5f, D_w=%.5f\n', ...
            col_labels{b}, R_count(1, b), ...
            TL_u(b), TL_v(b), TL_w(b), D_u(b), D_v(b), D_w(b));
end

%% Calculate FTLE 




%% -----------PLOTTING------------------

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

% Log-linear to check for exponential decay
figure;
semilogy(tau_vec, abs(Ruu), 'LineWidth', 2); hold on;
semilogy(tau_vec, abs(Rvv), 'LineWidth', 2);
semilogy(tau_vec, abs(Rww), 'LineWidth', 2);
xlabel('$\tau$'); ylabel('$|R_{ii}(\tau)|$');
legend({'$R_{uu}$', '$R_{vv}$', '$R_{ww}$'}, 'Location', 'best');
title('Lagrangian autocorrelation (log scale)');
grid on;

% Exponential fit for comparison
tau_fit = tau_vec(tau_vec < 2*TL_u);
R_exp_fit = exp(-tau_fit / TL_u);
figure;
plot(tau_vec, Ruu, 'b-', 'LineWidth', 2); hold on;
plot(tau_fit, R_exp_fit, 'r--', 'LineWidth', 2);
yline(0, 'k--');
xlabel('$\tau$'); ylabel('$R_{uu}(\tau)$');
legend({'Data', sprintf('$e^{-\\tau/T_L}$, $T_L = %.3f$', TL_u)}, 'Location', 'best');
title('Lagrangian $R_{uu}$ with exponential fit');
grid on;

%% Spatial banded autocorrelations 
% --- R_uu, R_vv, R_ww per y-band ---
figure;
tiledlayout(1, 3);

nexttile;
hold on;
for b = 1:n_cols
    if R_count(1, b) == 0, continue; end
    if b == n_cols
        plot(tau_vec, Ruu(:, b), 'k--', 'LineWidth', 2, 'DisplayName', col_labels{b});
    else
        plot(tau_vec, Ruu(:, b), 'LineWidth', 2, 'DisplayName', col_labels{b});
    end
end
yline(0, 'k:', 'HandleVisibility', 'off');
xlabel('$\tau$'); ylabel('$R_{uu}(\tau)$');
title('$R_{uu}$'); legend('Location', 'best'); grid on;

nexttile;
hold on;
for b = 1:n_cols
    if R_count(1, b) == 0, continue; end
    if b == n_cols
        plot(tau_vec, Rvv(:, b), 'k--', 'LineWidth', 2, 'DisplayName', col_labels{b});
    else
        plot(tau_vec, Rvv(:, b), 'LineWidth', 2, 'DisplayName', col_labels{b});
    end
end
yline(0, 'k:', 'HandleVisibility', 'off');
xlabel('$\tau$'); ylabel('$R_{vv}(\tau)$');
title('$R_{vv}$'); grid on;

nexttile;
hold on;
for b = 1:n_cols
    if R_count(1, b) == 0, continue; end
    if b == n_cols
        plot(tau_vec, Rww(:, b), 'k--', 'LineWidth', 2, 'DisplayName', col_labels{b});
    else
        plot(tau_vec, Rww(:, b), 'LineWidth', 2, 'DisplayName', col_labels{b});
    end
end
yline(0, 'k:', 'HandleVisibility', 'off');
xlabel('$\tau$'); ylabel('$R_{ww}(\tau)$');
title('$R_{ww}$'); grid on;

% --- Integral timescale vs y ---
figure;
plot(y_band_centers, TL_u(1:n_bands), 'o-', 'LineWidth', 2, 'DisplayName', '$T_L(u)$'); hold on;
plot(y_band_centers, TL_v(1:n_bands), 's-', 'LineWidth', 2, 'DisplayName', '$T_L(v)$');
plot(y_band_centers, TL_w(1:n_bands), 'd-', 'LineWidth', 2, 'DisplayName', '$T_L(w)$');
% Add horizontal lines for the "all tracks" values
yline(TL_u(end), 'b--', 'HandleVisibility', 'off');
yline(TL_v(end), 'r--', 'HandleVisibility', 'off');
yline(TL_w(end), 'Color', [0.93 0.69 0.13], 'LineStyle', '--', 'HandleVisibility', 'off');
xlabel('$y_0$'); ylabel('$T_L$');
legend('Location', 'best');
title('Lagrangian integral timescale vs initial $y$-position');
grid on;

% --- Eddy diffusivity vs y ---
figure;
plot(y_band_centers, D_u(1:n_bands), 'o-', 'LineWidth', 2, 'DisplayName', '$D_u$'); hold on;
plot(y_band_centers, D_v(1:n_bands), 's-', 'LineWidth', 2, 'DisplayName', '$D_v$');
plot(y_band_centers, D_w(1:n_bands), 'd-', 'LineWidth', 2, 'DisplayName', '$D_w$');
xlabel('$y_0$'); ylabel('$D_t$');
legend('Location', 'best');
title('Eddy diffusivity vs initial $y$-position');
grid on;
