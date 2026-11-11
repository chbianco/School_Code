%% MEAM 5380 HW7
%% Preamble
clc; clear all; close all;
set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextFontSize', 10);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);

%% Load Data
[u,v,w,p,x,y,z,nx,ny,nz] = load_channel_data('channel_ReTau180_uvwp.txt');
ny = numel(y);                     % lower half only (97 points)
 
SAVE_FIGS = true;                  % write each figure to PNG

set(groot, 'defaultTextInterpreter', 'latex', ...
           'defaultAxesTickLabelInterpreter', 'latex');
try  
    set(groot, 'defaultColorbarTickLabelInterpreter', 'latex', ...
               'defaultLegendInterpreter', 'latex');
end
Re_tau = 180;  nu = 1/Re_tau;  utau = 1;  dnu = nu/utau;   % viscous length
xp = x/dnu;  yp = y/dnu;  zp = z/dnu;                     % wall units
 
% Fluctuations
ix = 1:nx-1;  iz = 1:nz-1;
pmean = @(f) reshape(mean(mean(f(ix,:,iz),1),3), 1, []);   % 1 x ny
uf = u - reshape(pmean(u), 1, ny, 1);
vf = v - reshape(pmean(v), 1, ny, 1);
wf = w - reshape(pmean(w), 1, ny, 1);
pf = p - reshape(pmean(p), 1, ny, 1);
% u_tau = rho = 1  ->  u'+ = u', p'+ = p'
 
%% B) u'+(x+,z+) at five wall-normal locations
yp_target = [3, 15, 35, 50, 180];        % sublayer, buffer, log, log, centerline
lbl = {'viscous sublayer','buffer layer','log layer','log layer','centerline'};
figure('Name','B: u''+ planes','Position',[50 50 900 1200]);
for n = 1:5
    [~, j] = min(abs(yp - yp_target(n)));
    uplane = squeeze(uf(:,j,:));                       % nx x nz
    tmp = uplane(ix,iz);  s = std(tmp(:));
    subplot(5,1,n);
    pcolor(xp, zp, uplane'); shading flat; daspect([1 1 1]); colorbar;
    caxis(2.5*s*[-1 1]); colormap(jet);
    title(sprintf('$u^{\\prime +}$ at $y^+ = %.1f$ (%s), $u^{\\prime +}_{\\mathrm{rms}} = %.2f$', yp(j), lbl{n}, s));
    ylabel('$z^+$'); if n == 5, xlabel('$x^+$'); end
end
if SAVE_FIGS, print('-dpng','-r120','figB_uplanes.png'); end
 
%% C, D) Low-speed streaks in the buffer layer (y+ ~ 12)
[~, jb] = min(abs(yp - 12));
ub = squeeze(uf(:,jb,:));
figure('Name','C/D: streaks','Position',[50 50 1400 600]);
pcolor(xp, zp, ub'); shading flat; daspect([1 1 1]); colorbar; colormap(jet);
caxis([-4 4]); hold on;
contour(xp, zp, ub', [-2 -2], 'k', 'LineWidth', 1);   % outline low-speed streaks
set(gca,'XTick',0:100:1200,'YTick',0:100:600); grid on; set(gca,'Layer','top');
xlabel('$x^+$'); ylabel('$z^+$');
title(sprintf('$u^{\\prime +}$ at $y^+ = %.1f$; black: $u^{\\prime +} = -2$ (grid = 100 wall units)', yp(jb)));
if SAVE_FIGS, print('-dpng','-r120','figCD_streaks.png'); end
 
% Quantitative check of the visual estimates: two-point autocorrelations
% (periodic, via FFT) of u' on the buffer-layer plane.
U = ub(ix,iz);                                         % unique periodic points
Rz = mean(real(ifft(abs(fft(U,[],2)).^2,[],2)),1);  Rz = Rz/Rz(1);
Rx = mean(real(ifft(abs(fft(U,[],1)).^2,[],1)),2);  Rx = Rx/Rx(1);
dzp = (0:nz-2)*(zp(2)-zp(1));  dxp = (0:nx-2)*(xp(2)-xp(1));
half = 1:floor((nz-1)/2);
kmin = find(Rz(1:end-1) < 0 & diff(Rz) > 0, 1);       % first local minimum
fprintf('C) R_uu(dz+) = %s\n', mat2str(round(Rz(half)*100)/100));
fprintf('C) First R_uu(dz) minimum at dz+ = %.0f -> streak spacing ~ %.0f wall units\n', ...
        dzp(kmin), 2*dzp(kmin));
kx = find(Rx < 0.2, 1);
fprintf('D) R_uu(dx) drops below 0.2 at dx+ = %.0f\n', dxp(kx));
figure('Name','C/D: autocorrelation');
subplot(2,1,1); plot(dzp(half), Rz(half), 'o-'); grid on;
xlabel('$\Delta z^+$'); ylabel('$R_{uu}$'); title(sprintf('Spanwise, $y^+ = %.1f$', yp(jb)));
subplot(2,1,2); plot(dxp(1:floor((nx-1)/2)), Rx(1:floor((nx-1)/2)), 'o-'); grid on;
xlabel('$\Delta x^+$'); ylabel('$R_{uu}$'); title('Streamwise');
if SAVE_FIGS, print('-dpng','-r120','figCD_autocorr.png'); end
 
%% F) Velocity gradient tensor, Q, and vorticity
% x, z: periodic 2nd-order central differences.
% y   : 2nd-order central differences on the stretched grid (3-point
%       non-uniform stencil); 2nd-order one-sided at the wall and at the
%       top of the (truncated) half-channel domain.
dx = x(2)-x(1);  dz = z(2)-z(1);
ddx = @(f) cat(1, (circshift(f(ix,:,:),-1,1) - circshift(f(ix,:,:),1,1))/(2*dx), ...
                  (f(2,:,:) - f(nx-1,:,:))/(2*dx));
ddz = @(f) cat(3, (circshift(f(:,:,iz),-1,3) - circshift(f(:,:,iz),1,3))/(2*dz), ...
                  (f(:,:,2) - f(:,:,nz-1))/(2*dz));
 
h  = diff(y);                                  % (ny-1) x 1
h1 = reshape(h(1:end-1),1,[],1);  h2 = reshape(h(2:end),1,[],1);
cm = -h2./(h1.*(h1+h2));  c0 = (h2-h1)./(h1.*h2);  cp = h1./(h2.*(h1+h2));
a1 = h(1); a2 = h(2);                          % wall (forward)
b1 = h(end-1); b2 = h(end);                    % top (backward)
ddy = @(f) cat(2, ...
    -(2*a1+a2)/(a1*(a1+a2))*f(:,1,:) + (a1+a2)/(a1*a2)*f(:,2,:) - a1/(a2*(a1+a2))*f(:,3,:), ...
    cm.*f(:,1:end-2,:) + c0.*f(:,2:end-1,:) + cp.*f(:,3:end,:), ...
    b2/(b1*(b1+b2))*f(:,end-2,:) - (b1+b2)/(b1*b2)*f(:,end-1,:) + (2*b2+b1)/(b2*(b1+b2))*f(:,end,:));
 
% A{i,j} = du_i/dx_j (total velocity; mean shear contributes to Q too)
vel = {u, v, w};  A = cell(3,3);
for i = 1:3
    A{i,1} = ddx(vel{i});  A{i,2} = ddy(vel{i});  A{i,3} = ddz(vel{i});
end
clear vel
SS = zeros(size(u));  OO = zeros(size(u));
for i = 1:3
    for j = 1:3
        Sij = 0.5*(A{i,j} + A{j,i});  Oij = 0.5*(A{i,j} - A{j,i});
        SS = SS + Sij.^2;  OO = OO + Oij.^2;
    end
end
Q  = 0.5*(OO - SS);
wx = A{3,2} - A{2,3};      % dw/dy - dv/dz
wy = A{1,3} - A{3,1};      % du/dz - dw/dx
wz = A{2,1} - A{1,2};      % dv/dx - du/dy
divu = A{1,1} + A{2,2} + A{3,3};
clear A SS OO Sij Oij
 
% Wall units: Q+ = Q nu^2/utau^4, omega+ = omega nu/utau^2
Qp = Q*nu^2;  wxp = wx*nu;  wyp = wy*nu;  wzp = wz*nu;
fprintf('F) max|div u|/max|omega| = %.2e (sanity check)\n', max(abs(divu(:)))/max(abs(wx(:))));
 
%% G) yz plane at x = 5.592 (i = 90): Q+, in-plane streamlines, p'+
i0 = 90;  fprintf('G) x(%d) = %.3f\n', i0, x(i0));
[Zg, Yg] = meshgrid(zp, yp);                       % ny x nz
Qyz = squeeze(Qp(i0,:,:));  pyz = squeeze(pf(i0,:,:));
Vyz = squeeze(vf(i0,:,:));  Wyz = squeeze(wf(i0,:,:));
 
figure('Name','G: Q+ and streamlines','Position',[50 50 1400 600]);
pcolor(Zg, Yg, Qyz); shading interp; daspect([1 1 1]); colormap(hot); colorbar;
caxis([0 0.02]); hold on;
if exist('streamslice','file') == 2
    hs = streamslice(Zg, Yg, Wyz, Vyz, 4);
else                                               % Octave fallback
    [sz, sy] = meshgrid(zp(1:4:end), yp(3:4:end));
    hs = streamline(Zg, Yg, Wyz, Vyz, sz(:), sy(:));
end
set(hs, 'Color', [0.3 0.75 1], 'LineWidth', 0.8);
xlim([0 zp(end)]); ylim([0 yp(end)]);
xlabel('$z^+$'); ylabel('$y^+$'); title(sprintf('$Q^+$ and $(w^{\\prime}, v^{\\prime})$ streamlines, $x = %.3f$', x(i0)));
if SAVE_FIGS, print('-dpng','-r120','figG_Q_streamlines.png'); end
 
figure('Name','G: p''+','Position',[50 50 1400 600]);
pcolor(Zg, Yg, pyz); shading interp; daspect([1 1 1]); colormap(jet); colorbar;
caxis([-3 3]); hold on;
contour(Zg, Yg, Qyz, [0.005 0.005], 'k', 'LineWidth', 1);   % Q+ outline for reference
xlabel('$z^+$'); ylabel('$y^+$'); title(sprintf('$p^{\\prime +}$, $x = %.3f$ (black: $Q^+ = 0.005$)', x(i0)));
if SAVE_FIGS, print('-dpng','-r120','figG_pressure.png'); end
mask = Qyz > 0.005;
fprintf('G) mean p''+ where Q+>0.005: %.2f ; elsewhere: %.2f\n', mean(pyz(mask)), mean(pyz(~mask)));
 
%% H) 3D: Q isosurfaces colored by omega_x, plus u'+ = -2 streaks
Qrms   = sqrt(pmean(Qp.^2));                     % rms of Q+ on each y plane
Qthr   = max(Qrms);
fprintf('H) Q+_rms^max = %.4f (at y+ = %.1f)\n', Qthr, yp(Qrms == Qthr));
 
[X3, Y3, Z3] = meshgrid(x, y, z);                % ny x nx x nz
P = @(f) permute(f, [2 1 3]);
figure('Name','H: vortical structures');
fv = isosurface(X3, Y3, Z3, P(Qp), Qthr, P(wxp));
hQ = patch(fv, 'FaceColor', 'interp', 'EdgeColor', 'none');
colormap(jet); caxis([-0.2 0.2]); cb = colorbar; ylabel(cb, '$\omega_x^+$');
hold on;
hU = patch(isosurface(X3, Y3, Z3, P(uf), -2), 'FaceColor', [0.6 0.6 0.6], ...
           'EdgeColor', 'none', 'FaceAlpha', 0.25);
xlabel('$x$'); ylabel('$y$'); zlabel('$z$');
title(sprintf('$Q^+ = %.4f$ colored by $\\omega_x^+$; grey: $u^{\\prime +} = -2$', Qthr));
 
% camera / lighting (from template)
fig = gcf;
axes1 = get(fig,'CurrentAxes');
set(axes1,'CameraPosition',-1.35*[-12 18 -9],'CameraUpVector',[0 1 0], ...
    'CameraViewAngle',11,'DataAspectRatio',[1 1 1]);
camlight('headlight'); lighting gouraud;
ylim([0 1]); xlim([0 2*pi]); zlim([0 1*pi]);
set(gcf, 'Position', [100, 100, 1500, 1000]);
box on; try, set(gca,'BoxStyle','full'); end   % BoxStyle is MATLAB-only
if SAVE_FIGS, print('-dpng','-r100','figH_isosurfaces.png'); end
 
%% I) Supporting statistics for the structure characterization
core = Qp > Qthr;
nyc  = squeeze(sum(sum(core(ix,:,iz),1),3));     % points per y plane
fprintf('I) Fraction of Q+>thr points: y+<5: %.2f, 5-30: %.2f, 30-54: %.2f, >54: %.2f\n', ...
    sum(nyc(yp<5))/sum(nyc), sum(nyc(yp>=5 & yp<30))/sum(nyc), ...
    sum(nyc(yp>=30 & yp<54))/sum(nyc), sum(nyc(yp>=54))/sum(nyc));
[~, jpk] = max(nyc ./ max(nyc));
fprintf('I) Vortex-core volume fraction peaks at y+ = %.1f\n', yp(jpk));
fprintf('I) Within cores: P(omega_x>0) = %.2f\n', mean(wxp(core) > 0));
wzf = wzp - reshape(pmean(wzp), 1, ny, 1);       % remove mean-shear vorticity
fprintf('I) Within cores: <wx^2>:<wy^2>:<wz''^2> = %.3f : %.3f : %.3f\n', ...
    mean(wxp(core).^2), mean(wyp(core).^2), mean(wzf(core).^2));
figure('Name','I: vortex-core distribution');
plot(yp, nyc/((nx-1)*(nz-1)), 'o-'); grid on;
xlabel('$y^+$'); ylabel('volume fraction with $Q^+ > Q^+_{\mathrm{thr}}$');
if SAVE_FIGS, print('-dpng','-r120','figI_core_distribution.png'); end