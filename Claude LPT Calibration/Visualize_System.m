% =========================================================================
% Visualize_System.m
% =========================================================================
% PURPOSE:
%   Interactive 3-D visualisation of the calibrated camera system,
%   including:
%     • Camera positions and viewing frustums
%     • Glass wall planes with normal vectors
%     • Sample refracted ray paths
%     • Wand trajectory (scatter of 3-D wand midpoints)
%     • Reconstructed particle positions (if available)
%     • Per-camera reprojection error heat maps
% =========================================================================

clear; clc;
LPT_Config;

fprintf('=== Visualise Calibrated System ===\n\n');

resultsDir = cfg.resultsDir;

% Load data
cal  = load(fullfile(resultsDir, 'extrinsics.mat'));
cameras = cal.cameras;
nCams   = cfg.nCams;
n       = struct('air', cfg.n_air, 'glass', cfg.n_glass, 'water', cfg.n_water);

% Optional loads
hasBundleAdj  = exist(fullfile(resultsDir, 'bundle_adjustment.mat'), 'file');
hasRecon      = exist(fullfile(resultsDir, 'reconstruction.mat'),    'file');

if hasBundleAdj
    ba = load(fullfile(resultsDir, 'bundle_adjustment.mat'));
end
if hasRecon
    rc = load(fullfile(resultsDir, 'reconstruction.mat'));
end

% -------------------------------------------------------------------------
% Figure 1: Camera network and wall geometry
% -------------------------------------------------------------------------
fig1 = figure('Name','Camera Network & Wall Geometry', ...
    'Position',[50 50 1100 800], 'Color','w');
ax = axes('Parent',fig1);
hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
xlabel(ax,'X (m)'); ylabel(ax,'Y (m)'); zlabel(ax,'Z (m)');
title(ax,'Calibrated Camera Network');
view(ax,35,25);
lighting(ax,'gouraud');
camlight(ax,'headlight');

camColors = lines(nCams);

for c = 1:nCams
    R = cameras(c).R;
    t = cameras(c).t;
    C = -R' * t;   % Camera centre in world coords

    % Draw camera body (small box)
    drawCameraBody(ax, C, R, 0.03, camColors(c,:));

    % Draw frustum
    K     = cameras(c).K;
    imgSz = cfg.imageSize;
    drawFrustum(ax, C, R, K, imgSz, 0.15, camColors(c,:));

    % Label
    text(ax, C(1)+0.02, C(2)+0.02, C(3)+0.02, ...
        sprintf('Cam %d', c), 'Color', camColors(c,:), 'FontWeight','bold','FontSize',11);
end

% Draw glass walls
for c = 1:nCams
    wall = cameras(c).wall;
    drawWallPlane(ax, wall.normal(:), wall.point(:), wall.thickness, ...
                  [0.5 0.8 1.0], 0.2);
end

% Draw water volume (bounding box, schematic)
% Estimate from camera positions
allC = zeros(nCams,3);
for c = 1:nCams, allC(c,:) = (-cameras(c).R' * cameras(c).t)'; end
volMin = min(allC) - 0.1;
volMax = max(allC) + 0.1;
drawBoundingBox(ax, volMin, volMax, [0.2 0.4 0.8], 0.08);

% Draw sample refracted rays (from camera 1, corner of image)
drawSampleRays(ax, cameras, cfg, n, camColors);

saveas(fig1, fullfile(resultsDir, 'viz_camera_network.png'));
fprintf('Saved: viz_camera_network.png\n');

% -------------------------------------------------------------------------
% Figure 2: Wand coverage (Phase 4 wand midpoints)
% -------------------------------------------------------------------------
if hasBundleAdj
    fig2 = figure('Name','Wand Coverage', 'Position',[200 50 900 700], 'Color','w');
    ax2  = axes('Parent',fig2);
    hold(ax2,'on'); grid(ax2,'on'); axis(ax2,'equal');
    xlabel(ax2,'X (m)'); ylabel(ax2,'Y (m)'); zlabel(ax2,'Z (m)');
    title(ax2,'Wand Calibration Coverage');
    view(ax2,35,25);

    % Extract wand midpoints from bundle adjustment
    wandPts = ba.wandPts_opt;
    nWand   = size(wandPts,1);
    mids    = wandPts(:,4:6);

    scatter3(ax2, mids(:,1), mids(:,2), mids(:,3), 20, ...
        1:nWand, 'filled','MarkerFaceAlpha',0.6);
    colormap(ax2, parula);
    cb = colorbar(ax2); cb.Label.String = 'Frame index';
    title(ax2, sprintf('Wand midpoints — %d frames', nWand));

    % Camera positions
    for c = 1:nCams
        C = -cameras(c).R' * cameras(c).t;
        scatter3(ax2, C(1), C(2), C(3), 120, camColors(c,:), 'filled', ...
            'Marker','p','MarkerEdgeColor','k');
        text(ax2, C(1), C(2), C(3)+0.03, sprintf('Cam%d',c), ...
            'Color',camColors(c,:),'FontWeight','bold');
    end

    saveas(fig2, fullfile(resultsDir, 'viz_wand_coverage.png'));
    fprintf('Saved: viz_wand_coverage.png\n');
end

% -------------------------------------------------------------------------
% Figure 3: Reprojection error heat maps
% -------------------------------------------------------------------------
if hasBundleAdj
    fig3 = figure('Name','Reprojection Error Maps', ...
        'Position',[300 50 1400 800], 'Color','w');

    for c = 1:nCams
        ax3 = subplot(2, ceil(nCams/2), c, 'Parent', fig3);
        hold(ax3,'on'); box(ax3,'on');
        xlim(ax3,[0 cfg.imageSize(1)]); ylim(ax3,[0 cfg.imageSize(2)]);
        set(ax3,'YDir','reverse');
        axis(ax3,'equal');

        % Collect all projected points for this camera
        allU = []; allV = []; allE = [];
        validFrames = find(arrayfun(@(f) ~isempty(ba.frameErrors) && f<=numel(ba.frameErrors), 1:nFrames));

        % Re-project wand endpoints
        wandPts_ = ba.wandPts_opt;
        wL = cfg.wandLength;
        for fi = 1:size(wandPts_,1)
            wp = wandPts_(fi,:);
            aa = wp(1:3); ang = norm(aa);
            if ang < 1e-10, wdir=[1;0;0];
            else, wdir=axang2rotm_local([aa/ang,ang])*[1;0;0]; end
            wdir = wdir/norm(wdir);
            mid  = wp(4:6)';
            XA   = mid-(wL/2)*wdir; XB = mid+(wL/2)*wdir;

            pose = struct('R',cameras(c).R,'t',cameras(c).t);
            wg   = struct('normal',cameras(c).wall.normal(:),...
                          'point',cameras(c).wall.point(:),...
                          'thickness',cameras(c).wall.thickness);

            uvA = projectPointRefractive(XA, pose, cameras(c).K, cameras(c).distCoeffs, wg, n);
            uvB = projectPointRefractive(XB, pose, cameras(c).K, cameras(c).distCoeffs, wg, n);

            if ~any(isnan(uvA)) && ~any(isnan(uvB))
                allU = [allU; uvA(1); uvB(1)]; %#ok<AGROW>
                allV = [allV; uvA(2); uvB(2)]; %#ok<AGROW>
                allE = [allE; ba.frameErrors(fi); ba.frameErrors(fi)]; %#ok<AGROW>
            end
        end

        if ~isempty(allU)
            scatter(ax3, allU, allV, 15, allE, 'filled','MarkerFaceAlpha',0.7);
            colormap(ax3, hot); clim(ax3,[0, cfg.ba.reprojThresh]);
            cb3 = colorbar(ax3); cb3.Label.String = 'Reproj. err (px)';
        end
        title(ax3, sprintf('Camera %d  (RMS = %.3f px)', c, cameras(c).reprojRMS));
        xlabel(ax3,'u (px)'); ylabel(ax3,'v (px)');
    end
    sgtitle(fig3,'Reprojection Error Spatial Distribution');
    saveas(fig3, fullfile(resultsDir, 'viz_reproj_heatmaps.png'));
    fprintf('Saved: viz_reproj_heatmaps.png\n');
end

% -------------------------------------------------------------------------
% Figure 4: Reconstructed particles (if available)
% -------------------------------------------------------------------------
if hasRecon
    recon = rc.recon;
    allPts = [];
    for f = 1:numel(recon)
        if ~isempty(recon(f).positions)
            valid = ~any(isnan(recon(f).positions),2);
            allPts = [allPts; recon(f).positions(valid,:)]; %#ok<AGROW>
        end
    end

    if ~isempty(allPts)
        fig4 = figure('Name','Reconstructed Particles','Position',[400 50 900 700],'Color','w');
        ax4  = axes('Parent',fig4);
        scatter3(ax4, allPts(:,1), allPts(:,2), allPts(:,3), 5, ...
            allPts(:,3), 'filled', 'MarkerFaceAlpha',0.4);
        colormap(ax4, cool);
        colorbar(ax4);
        xlabel(ax4,'X (m)'); ylabel(ax4,'Y (m)'); zlabel(ax4,'Z (m)');
        title(ax4, sprintf('Reconstructed particles — %d points', size(allPts,1)));
        grid(ax4,'on'); axis(ax4,'equal'); view(ax4,35,25);
        saveas(fig4, fullfile(resultsDir, 'viz_particles.png'));
        fprintf('Saved: viz_particles.png\n');
    end
end

fprintf('\n[Visualise] Complete.\n');


% =========================================================================
%  DRAWING HELPERS
% =========================================================================
function drawCameraBody(ax, C, R, sz, col)
    corners = sz * (R' * [-1 1 1 -1 -1 1 1 -1; -1 -1 1 1 -1 -1 1 1; 0 0 0 0 1 1 1 1] - [0;0;0.5]);
    corners = corners + C;
    faces   = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    patch(ax,'Vertices',corners','Faces',faces,'FaceColor',col, ...
        'FaceAlpha',0.6,'EdgeColor',col*0.7,'LineWidth',0.5);
end

function drawFrustum(ax, C, R, K, imgSz, depth, col)
    % Project image corners into world
    corners_img = [0 imgSz(1) imgSz(1) 0; 0 0 imgSz(2) imgSz(2)];
    dirs = zeros(3,4);
    for i = 1:4
        xn = (corners_img(1,i) - K(1,3)) / K(1,1);
        yn = (corners_img(2,i) - K(2,3)) / K(2,2);
        dirs(:,i) = R' * [xn;yn;1];
        dirs(:,i) = dirs(:,i) / norm(dirs(:,i));
    end
    for i = 1:4
        P = C + depth * dirs(:,i);
        plot3(ax,[C(1),P(1)],[C(2),P(2)],[C(3),P(3)],'Color',col,'LineWidth',1,'LineStyle','--');
    end
    % Connect far corners
    Ps = C + depth * dirs;
    plot3(ax,[Ps(1,:),Ps(1,1)],[Ps(2,:),Ps(2,1)],[Ps(3,:),Ps(3,1)], ...
        'Color',col,'LineWidth',1);
end

function drawWallPlane(ax, normal, point, thickness, col, alpha)
    % Draw two parallel planes (outer and inner glass faces)
    normal  = normal / norm(normal);
    outerPt = point + thickness * normal;

    % Build orthogonal basis in the plane
    if abs(normal(1)) < 0.9
        u = cross(normal, [1;0;0]); else, u = cross(normal, [0;1;0]);
    end
    u = u/norm(u); v = cross(normal,u);

    sz = 0.25;
    corners_inner = point + sz*[-u+v, u+v, u-v, -u-v];
    corners_outer = outerPt + sz*[-u+v, u+v, u-v, -u-v];

    patch(ax,'Vertices',corners_inner','Faces',[1 2 3 4],'FaceColor',col, ...
        'FaceAlpha',alpha,'EdgeColor',col*0.6);
    patch(ax,'Vertices',corners_outer','Faces',[1 2 3 4],'FaceColor',col*0.8, ...
        'FaceAlpha',alpha/2,'EdgeColor',col*0.5,'LineStyle',':');

    % Normal arrow
    quiver3(ax,point(1),point(2),point(3), normal(1)*0.05,normal(2)*0.05,normal(3)*0.05, ...
        'Color',col*0.5,'LineWidth',2,'MaxHeadSize',0.5);
end

function drawBoundingBox(ax, mn, mx, col, alpha)
    X = [mn(1) mx(1)]; Y = [mn(2) mx(2)]; Z = [mn(3) mx(3)];
    verts = [X([1 2 2 1 1 2 2 1])', Y([1 1 2 2 1 1 2 2])', Z([1 1 1 1 2 2 2 2])'];
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 3 4 8 7; 1 4 8 5; 2 3 7 6];
    patch(ax,'Vertices',verts,'Faces',faces,'FaceColor',col,'FaceAlpha',alpha, ...
        'EdgeColor',col,'LineStyle','--','LineWidth',0.5);
end

function drawSampleRays(ax, cameras, cfg, n, camColors)
    % Draw 4 sample refracted rays, one from each camera centre
    for c = 1:numel(cameras)
        R = cameras(c).R;
        t = cameras(c).t;
        C = -R' * t;
        K = cameras(c).K;

        % Principal ray direction
        d_air = R' * [0;0;1];

        wall = cameras(c).wall;
        wg   = struct('normal', wall.normal(:), 'point', wall.point(:), ...
                       'thickness', wall.thickness);

        ray = refractRay(C, d_air, wg.normal, wg.point, wg.thickness, ...
                          n.air, n.glass, n.water);
        if ~ray.valid, continue; end

        % Air segment
        plot3(ax, [C(1), ray.hitAir(1)], [C(2), ray.hitAir(2)], [C(3), ray.hitAir(3)], ...
            '-','Color',camColors(c,:),'LineWidth',1.5);
        % Glass segment
        plot3(ax, [ray.hitAir(1), ray.hitGlass(1)], [ray.hitAir(2), ray.hitGlass(2)], ...
            [ray.hitAir(3), ray.hitGlass(3)], '-','Color',camColors(c,:)*0.7,'LineWidth',2);
        % Water segment (extend 0.1m into water)
        endPt = ray.origin + 0.10 * ray.dir;
        plot3(ax, [ray.origin(1), endPt(1)], [ray.origin(2), endPt(2)], ...
            [ray.origin(3), endPt(3)], '--','Color',camColors(c,:),'LineWidth',1.5);
    end
end

function R = axang2rotm_local(aa)
    if numel(aa)==3
        angle=norm(aa); if angle<1e-10, R=eye(3); return; end
        ax=aa(:)/angle;
    else
        ax=aa(1:3)'; angle=aa(4);
    end
    c=cos(angle);s=sin(angle);t=1-c;x=ax(1);y=ax(2);z=ax(3);
    R=[t*x*x+c,t*x*y-s*z,t*x*z+s*y;t*x*y+s*z,t*y*y+c,t*y*z-s*x;t*x*z-s*y,t*y*z+s*x,t*z*z+c];
end
