%% Preamble
close all; clear variables; clc;
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

%% Define bin grid
gridX = linspace(0,    8*pi, 64);   % 64 bins in x
gridZ = linspace(0,    3*pi, 33);   % 32 bins in z

ny_bins = 128; %128 bins in y
theta = linspace(0, pi, ny_bins + 1);
gridY = -cos(theta);   %Cluster bins at walls

nx = numel(gridX) - 1;
ny = numel(gridY) - 1;
nz = numel(gridZ) - 1;

%% Chunked binning: accumulate sums over track chunks
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

%% Visualization
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

%% Reynolds Stresses 
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

%% Reynolds stress plots
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

%% Dispersive fluxes 

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

%% Dispersive fluxes vs Reynolds stresses
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
semilogy(yc, abs(uu_disp), 'DisplayName', '$|\langle\tilde u \tilde u\rangle|$'); hold on;
semilogy(yc, abs(vv_disp), 'DisplayName', '$|\langle\tilde v \tilde v\rangle|$');
semilogy(yc, abs(ww_disp), 'DisplayName', '$|\langle\tilde w \tilde w\rangle|$');
semilogy(yc, abs(uv_disp), 'DisplayName', '$|\langle\tilde u \tilde v\rangle|$');
xlabel('$y$'); ylabel('Dispersive flux (abs)');
title('Dispersive flux magnitudes (log scale)');
legend('Interpreter','latex','Location','best'); grid on;

%% Testing/comparison

[~, iy_mid] = min(abs(yc - 0));    % centerline
[~, iy_near] = min(abs(yc - 0.9)); % near wall

figure;
subplot(1,2,1);
imagesc(xc, zc, squeeze(U_tilde(:, iy_mid, :))');
axis xy equal tight; colorbar;
xlabel('$x$'); ylabel('$z$');
title(sprintf('$\\tilde u$ at $y = %.2f$ (centerline)', yc(iy_mid)));

subplot(1,2,2);
imagesc(xc, zc, squeeze(U_tilde(:, iy_near, :))');
axis xy equal tight; colorbar;
xlabel('$x$'); ylabel('$z$');
title(sprintf('$\\tilde u$ at $y = %.2f$ (near wall)', yc(iy_near)));