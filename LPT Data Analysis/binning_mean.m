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
gridX = linspace(0,    8*pi, 129);   % 64 bins in x
gridZ = linspace(0,    3*pi, 65);   % 32 bins in z

ny_bins = 256; %128 bins in y
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
plot(yc, uu_disp, 'DisplayName', '$|\langle\tilde u \tilde u\rangle|$'); hold on;
plot(yc, vv_disp, 'DisplayName', '$|\langle\tilde v \tilde v\rangle|$');
plot(yc, ww_disp, 'DisplayName', '$|\langle\tilde w \tilde w\rangle|$');
plot(yc, uv_disp, 'DisplayName', '$|\langle\tilde u \tilde v\rangle|$');
xlabel('$y$'); ylabel('Dispersive flux (abs)');
title('Dispersive flux magnitudes (log scale)');
legend('Interpreter','latex','Location','best'); grid on;

%% Pressure Poisson

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

fprintf('Wavenumber vectors: kx(%d), ky(%d), kz(%d)\n', ...
        length(kx_vec), length(ky_vec), length(kz_vec));

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


%Inv FFT
Pmean = real(ifftn(P_hat));

% Check that imaginary part is negligible
imag_residual = max(abs(imag(ifftn(P_hat))), [], 'all');
fprintf('Max imaginary residual: %.2e (should be ~1e-14)\n', imag_residual);

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

%% Validate velocity 
% token
authkey = 'edu.jhu.pha.turbulence.testing-201406';  
% dataset (channel flow re_tau ~ 1000)
dataset =  'channel';
% variable of interest
variable = 'velocity';
% temporal interpolation
temporal_method = 'pchip'; 
% spatial interpolation
spatial_method = 'lag4';
% spatial data capture
spatial_operator  = 'field';

%Points per query 
ny_dns = 64;
nz_dns = 8;
nx_dns = 8;


theta_dns = linspace(0, pi, ny_dns);
gridY_dns = -cos(theta_dns); 

gridX_dns = linspace(2, 8*pi - 2, nx_dns);   % avoid edges
gridZ_dns = linspace(0.5, 3*pi - 0.5, nz_dns);

n_points = nx_dns * ny_dns * nz_dns;
points = zeros(n_points, 3);
idx = 0;
for i = 1:nx_dns
    for j = 1:ny_dns
        for k = 1:nz_dns
            idx = idx + 1;
            points(idx, 1) = gridX_dns(i);
            points(idx, 2) = gridY_dns(j);
            points(idx, 3) = gridZ_dns(k);
        end
    end
end

n_times = 20; 
times = linspace(0.5, 25, n_times);

sumU_dns  = zeros(ny_dns, 1);
sumV_dns  = zeros(ny_dns, 1);
sumW_dns  = zeros(ny_dns, 1);
sumUU_dns = zeros(ny_dns, 1);
sumVV_dns = zeros(ny_dns, 1);
sumWW_dns = zeros(ny_dns, 1);
sumUV_dns = zeros(ny_dns, 1);
n_samples = zeros(ny_dns, 1);

fprintf('Querying JHTDB: %d times x %d points = %d total queries...\n', ...
        n_times, n_points, n_times);

for it = 1:n_times
    t_query = times(it);
    fprintf('  time %d/%d (t = %.2f)...', it, n_times, t_query);

    vel = getData(authkey, dataset, variable, t_query, ...
                  temporal_method, spatial_method, spatial_operator, points);

    % vel is (n_points x 3): columns are u, v, w
    % Reshape to (nx_dns x ny_dns x nz_dns x 3)
    u_snap = reshape(vel(:,1), [nz_dns, ny_dns, nx_dns]);
    v_snap = reshape(vel(:,2), [nz_dns, ny_dns, nx_dns]);
    w_snap = reshape(vel(:,3), [nz_dns, ny_dns, nx_dns]);

    % Average over x and z dimensions for this snapshot
    u_prof = squeeze(mean(u_snap, [1 3]));   % (ny_dns x 1)
    v_prof = squeeze(mean(v_snap, [1 3]));
    w_prof = squeeze(mean(w_snap, [1 3]));

    % Accumulate first moments
    sumU_dns = sumU_dns + u_prof;
    sumV_dns = sumV_dns + v_prof;
    sumW_dns = sumW_dns + w_prof;

    % Accumulate second moments (per-point, then average over x,z)
    uu_prof = squeeze(mean(u_snap.^2, [1 3]));
    vv_prof_snap = squeeze(mean(v_snap.^2, [1 3]));
    ww_prof_snap = squeeze(mean(w_snap.^2, [1 3]));
    uv_prof_snap = squeeze(mean(u_snap .* v_snap, [1 3]));

    sumUU_dns = sumUU_dns + uu_prof;
    sumVV_dns = sumVV_dns + vv_prof_snap;
    sumWW_dns = sumWW_dns + ww_prof_snap;
    sumUV_dns = sumUV_dns + uv_prof_snap;

    n_samples = n_samples + 1;

    fprintf(' done\n');
end

%Compute DNS statistics
Umean_dns = sumU_dns ./ n_samples;
Vmean_dns = sumV_dns ./ n_samples;
Wmean_dns = sumW_dns ./ n_samples;

% Reynolds stresses: <u'u'> = <uu> - <u>^2
uu_dns = sumUU_dns ./ n_samples - Umean_dns.^2;
vv_dns = sumVV_dns ./ n_samples - Vmean_dns.^2;
ww_dns = sumWW_dns ./ n_samples - Wmean_dns.^2;
uv_dns = sumUV_dns ./ n_samples - Umean_dns .* Vmean_dns;

yc_dns = gridY_dns(:);

%% Comparison plots
% ---- Mean velocity ----
figure;
plot(yc, Uprofile, 'b-', 'LineWidth', 2); hold on;
plot(yc_dns, Umean_dns, 'r--', 'LineWidth', 2);
xlabel('$y$'); ylabel('$\langle u \rangle$');
legend({'LPT', 'DNS'}, 'Location', 'best');
title('Mean streamwise velocity');
grid on;

% ---- Reynolds stresses ----
figure;
tiledlayout(2, 2);

nexttile;
plot(yc, uu_prof, 'b-', 'LineWidth', 2); hold on;
plot(yc_dns, uu_dns, 'r--', 'LineWidth', 2);
xlabel('$y$'); title('$\overline{u''u''}$');
legend({'LPT', 'DNS'}, 'Location', 'best'); grid on;

nexttile;
plot(yc, vv_prof, 'b-', 'LineWidth', 2); hold on;
plot(yc_dns, vv_dns, 'r--', 'LineWidth', 2);
xlabel('$y$'); title('$\overline{v''v''}$');
grid on;

nexttile;
plot(yc, ww_prof, 'b-', 'LineWidth', 2); hold on;
plot(yc_dns, ww_dns, 'r--', 'LineWidth', 2);
xlabel('$y$'); title('$\overline{w''w''}$');
grid on;

nexttile;
plot(yc, uv_prof, 'b-', 'LineWidth', 2); hold on;
plot(yc_dns, uv_dns, 'r--', 'LineWidth', 2);
xlabel('$y$'); title('$\overline{u''v''}$');
grid on;

% ---- TKE ----
tke_dns = 0.5 * (uu_dns + vv_dns + ww_dns);

figure;
plot(yc, tke_prof, 'b-', 'LineWidth', 2); hold on;
plot(yc_dns, tke_dns, 'r--', 'LineWidth', 2);
xlabel('$y$'); ylabel('$k$');
legend({'LPT', 'DNS'}, 'Location', 'best');
title('Turbulent kinetic energy');
grid on;

fprintf('\nDNS comparison complete.\n');
fprintf('DNS: u_tau = 0.0499, nu = 5e-5, Re_tau ~ 1000\n');
fprintf('DNS samples: %d times x %d (x,z) points = %d per y-location\n', ...
        n_times, nx_dns*nz_dns, n_times*nx_dns*nz_dns);
