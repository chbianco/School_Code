%% Preamble
close all; clc;
clearvars -except tracks t x y z
set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);

%% Loading file
file = 'jhtb_long_LPT.mat'; %Tracks to analyze. 
%Should be a mat file containing struct LPT with column vectors t, x, y, and z

%% Load tracks
if exist('tracks', 'var') 
    t = tracks.t;
    x = tracks.x;
    y = tracks.y;
    z = tracks.z;
    [nT, nTracks] = size(x);
    fprintf('Track data already loaded')
else
    tracks = load(file).LPT;
    t = tracks.t;
    x = tracks.x;
    y = tracks.y;
    z = tracks.z;
    [nT, nTracks] = size(x);
    fprintf('Loaded: %d samples x %d tracks (%.2f GB per array)\n', ...
        nT, nTracks, nT*nTracks*8/1e9);
end
%% User inputs
%Number of tracks to analyze. Set to nTracks to keep all tracks
N = nTracks;  

%Number of trajectories to plot
n = 100;

%Bounds of LPT view area, as a 2x1 (ie 0, 20)
Xlim = [0, 8*pi];
Ylim = [-1, 1];
Zlim = [0, 3*pi];

%Number of bins in x, y, and z. Evenly spaced
Xbin = 128;
Ybin = 128;
Zbin = 128;

%Tracks per chunk; tune for memory. 1000 works well
chunkSize = 1000;  

%Max lags for Lagrangian statistics
max_lag = 2500;


%% Plot n tracks colored by speed
speed = NaN(4000, n);

u_vals = NaN(4000, n); 
v_vals = NaN(4000, n);
w_vals = NaN(4000, n);

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
    
    %Velocities. CHANGE FOR OPENLPT DATA
    for j = 1: trk_lng - 1
        dt = t_trk(j+1) - t_trk(j);
        un(j) = (xn(j+1) - xn(j))/dt;
        vn(j) = (yn(j+1) - yn(j))/dt;
        wn(j) = (zn(j+1) - zn(j))/dt;
    end
    u_vals(1:trk_lng-1, k) = un';
    v_vals(1:trk_lng-1, k) = vn';
    w_vals(1:trk_lng-1, k) = wn';
    
    speed(1:trk_lng-1, k) = sqrt(un.^2 + vn.^2 + wn.^2)';

    scatter3(xn(1:end-1),yn(1:end-1),zn(1:end-1), 20, speed(1:trk_lng-1,k), 'filled')


end 

c = colorbar;
c.Label.String = 'Speed';
c.Label.FontName = 'Latex';
xlabel('x')
ylabel('y')
zlabel('z')
xlim([0 8*pi])
ylim([-1 1])
zlim([0 3*pi])
axis equal

view(3)

hold off

%% Keep N longest tracks and organize
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
    i1 = min(c*chunkSize, nTracks);

    xc_ = x(:, i0:i1);
    yc_ = y(:, i0:i1);
    zc_ = z(:, i0:i1);

    % Forward differences -> velocities at midpoints. CHANGE FOR REAL DATA
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
