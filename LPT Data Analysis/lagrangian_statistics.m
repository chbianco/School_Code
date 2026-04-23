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

%% -----------LAGRANGIAN STATISTICS---------------
%% Calculate lagrangian velocity autocorrelations 
% Compute R_uu(tau), R_vv(tau), R_ww(tau) from individual track time series.

%Max lag in time steps 
max_lag = 1500;

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

    % Central differences for interior, forward/backward at ends
    u_trk = gradient(xseg, dt_track);
    v_trk = gradient(yseg, dt_track);
    w_trk = gradient(zseg, dt_track);

    % Remove track mean to get fluctuations
    % Note: this is the Lagrangian mean along the track, not the Eulerian bin mean.
    % For homogeneous turbulence they're equivalent; for channel flow there's a
    % subtlety because the particle samples different y-locations over time.
    % MAY NEED TO CHANGE THIS LATER 
    u_trk = u_trk - mean(u_trk);
    v_trk = v_trk - mean(v_trk);
    w_trk = w_trk - mean(w_trk);

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

