% =========================================================================
% LPT_Config.m  —  Central configuration for 4-camera LPT calibration
% =========================================================================
% Edit this file before running any calibration phase.
% All physical units are in METRES unless otherwise noted.
%
% Workflow:
%   1.  Set parameters in this file
%   2.  Run Phase1_AirCalibration.m   (intrinsics, in air, checkerboard)
%   3.  Run Phase2_WandDetection.m    (detect LED wand points in images)
%   4.  Run Phase4_BundleAdjustment.m (solve extrinsics + refine everything)
%   5.  Run Phase5_Reconstruction.m   (triangulate 3-D particle positions)
%   6.  Run Visualize_System.m        (optional diagnostics / 3-D plots)
% =========================================================================

cfg = struct();

% -------------------------------------------------------------------------
% 1.  CAMERA SYSTEM
% -------------------------------------------------------------------------
cfg.nCams           = 4;            % number of cameras
cfg.imageSize       = [1280 1024];  % [width height] pixels  — adjust to your model
% Edgertronic SC2+  →  1920×1080 or 1280×1024 depending on crop/binning

% -------------------------------------------------------------------------
% 2.  REFRACTIVE INDICES
% -------------------------------------------------------------------------
cfg.n_air           = 1.000;        % air
cfg.n_glass         = 1.473;        % borosilicate / soda-lime glass (typical flume)
cfg.n_water         = 1.333;        % fresh water at ~20 °C

% -------------------------------------------------------------------------
% 3.  FLUME / GLASS WALL GEOMETRY
% -------------------------------------------------------------------------
% Each camera looks through one planar glass interface.
% Specify for each camera:
%   wallNormal  : outward-pointing unit normal of the wall (in world coords)
%                 e.g. side wall facing +X  →  [1 0 0]
%                      floor              →  [0 0 -1]   (z points up)
%   wallPoint   : any point on the inner (water-side) face [m]
%   glassThick  : thickness of that glass pane [m]
%
% You may use the same wall entry for multiple cameras if they share a wall.
%
% Example for a rectangular flume, cameras on the +X side wall only:
%   wallNormal = [1 0 0],  wallPoint = [0.010  0  0],  glassThick = 0.010
%
% IMPORTANT: set these carefully — errors here propagate directly into
%            the refraction-corrected ray model.

% -------------------------------------------------------------------------
% WALL GEOMETRY
% Coordinate system: origin at measurement volume centre
%   X: across tank (left wall to right wall)
%   Y: along tank length (away from observer)
%   Z: vertical (upward)
% -------------------------------------------------------------------------
% Left wall (cameras 1 and 2 shoot through this)
cfg.walls(1).normal    = [-1, 0, 0];
cfg.walls(1).point     = [-0.216, 0, 0];
cfg.walls(1).thickness = 0.00635;

cfg.walls(2).normal    = [-1, 0, 0];
cfg.walls(2).point     = [-0.216, 0, 0];
cfg.walls(2).thickness = 0.00635;

% Right wall (cameras 3 and 4 shoot through this)
cfg.walls(3).normal    = [+1, 0, 0];
cfg.walls(3).point     = [+0.216, 0, 0];
cfg.walls(3).thickness = 0.00635;

cfg.walls(4).normal    = [+1, 0, 0];
cfg.walls(4).point     = [+0.216, 0, 0];
cfg.walls(4).thickness = 0.00635;

%IMPORTANT: Change triangulateRefractive to account for flume vs fish tank

% -------------------------------------------------------------------------
% APPROXIMATE CAMERA POSES (for bundle adjustment initialisation)
%   C_approx: approximate camera centre in world coords [m]
%   target:   a point the camera is roughly aimed at [m]
% -------------------------------------------------------------------------
% Camera 1 — left side, near end
cfg.initPoses(1).C_approx = [-1.270,  0.025,  0.000];
cfg.initPoses(1).target   = [ 0.000,  0.000,  0.000];

% Camera 2 — left side, far end
cfg.initPoses(2).C_approx = [-1.270,  0.508,  0.000];
cfg.initPoses(2).target   = [ 0.000,  0.000,  0.000];

% Camera 3 — right side, near end
cfg.initPoses(3).C_approx = [ 1.181, -0.051,  0.000];
cfg.initPoses(3).target   = [ 0.000,  0.000,  0.000];

% Camera 4 — right side, far end
cfg.initPoses(4).C_approx = [ 1.181,  0.305,  0.000];
cfg.initPoses(4).target   = [ 0.000,  0.000,  0.000];

% -------------------------------------------------------------------------
% 4.  WAND GEOMETRY
% -------------------------------------------------------------------------
cfg.wandLength      = 0.120;        % distance between the two LED tips [m]
                                    % Measure with calipers — accuracy matters!

% -------------------------------------------------------------------------
% 5.  DATA PATHS
% -------------------------------------------------------------------------
cfg.airCalibDir     = 'data\air_calibration';      % folder of checkerboard images per camera
                                              % sub-folders: cam1/ cam2/ cam3/ cam4/
cfg.wandImageDir    = 'data\wand_calibration';    % synchronised wand image sets
                                              % sub-folders: cam1/ cam2/ cam3/ cam4/
cfg.resultsDir      = 'results';             % all outputs written here

% -------------------------------------------------------------------------
% 6.  CHECKERBOARD (Phase 1)
% -------------------------------------------------------------------------
cfg.cbSquareSize    = 0.0192;        % physical size of one square [m]
cfg.cbBoardSize     = [9 12];        % inner corners [cols rows]

% -------------------------------------------------------------------------
% 7.  LED DETECTION (Phase 2)
% -------------------------------------------------------------------------
cfg.led.minArea         = 30;        % min blob area [px²]
cfg.led.maxArea         = 2200;      % max blob area [px²]
cfg.led.intensityPct    = 90;       % percentile threshold for bright-pixel mask
cfg.led.gaussSigma      = 2.0;      % Gaussian pre-blur sigma [px]
cfg.led.subpixel        = true;     % use intensity-weighted centroid refinement

% -------------------------------------------------------------------------
% 8.  BUNDLE ADJUSTMENT (Phase 4)
% -------------------------------------------------------------------------
cfg.ba.maxIter          = 200;
cfg.ba.fTol             = 1e-10;
cfg.ba.xTol             = 1e-10;
cfg.ba.reprojThresh     = 1.5;      % RANSAC-style outlier rejection [px]
cfg.ba.fixIntrinsics    = true;     % keep Phase-1 intrinsics fixed (recommended)
cfg.ba.refineWallNormals = false;   % allow small perturbations to wall normals
                                    % set true only if wall geometry is uncertain

% -------------------------------------------------------------------------
% 9.  RECONSTRUCTION (Phase 5)
% -------------------------------------------------------------------------
cfg.recon.minCams       = 2;        % minimum cameras that must see a particle
cfg.recon.maxReprojErr  = 2.0;      % discard triangulations above this [px]

%-----------------
% 9.5 Adding utils directory
%-----------------
addpath(fullfile(pwd, 'utils'));

% -------------------------------------------------------------------------
% Save config to results dir
% -------------------------------------------------------------------------
if ~exist(cfg.resultsDir, 'dir'), mkdir(cfg.resultsDir); end
save(fullfile(cfg.resultsDir, 'LPT_config.mat'), 'cfg');
fprintf('[Config] Saved to %s\n', fullfile(cfg.resultsDir, 'LPT_config.mat'));
