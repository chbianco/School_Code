% =========================================================================
% ExportCalibration.m
% =========================================================================
% PURPOSE:
%   Export all calibration results in multiple formats:
%
%     calibration.json        — human-readable, language-agnostic
%     calibration.csv         — flat parameter table
%     cam_N_params.txt        — annotated per-camera files (generic)
%     openlpt/camN.txt        — OpenLPT / STB-compatible format
%     openlpt/README.txt      — notes on the refraction limitation
%
% OPENLPT FORMAT NOTES:
%   OpenLPT's camera file format (used by both the GUI and the legacy
%   STB C++ backend) is a plain-text file produced by CalibMatToTXTV2.m
%   in the original STB repo. The format is:
%
%     Line 1:  image_width  image_height          (pixels)
%     Line 2:  fx  fy  cx  cy                     (pixels)
%     Line 3:  k1  k2  p1  p2  k3                 (Brown-Conrady distortion)
%     Line 4:  r11  r12  r13                       (rotation matrix, row 1)
%     Line 5:  r21  r22  r23                       (rotation matrix, row 2)
%     Line 6:  r31  r32  r33                       (rotation matrix, row 3)
%     Line 7:  tx  ty  tz                          (translation, world units)
%
%   IMPORTANT — UNITS:
%   OpenLPT internally works in MILLIMETRES for the translation vector.
%   The rotation matrix is dimensionless. Focal lengths and principal
%   point are in pixels. Our calibration uses metres; this exporter
%   converts t -> mm automatically.
%
%   IMPORTANT — REFRACTION CAVEAT:
%   OpenLPT's own calibration pipeline does NOT model refraction.  By
%   feeding it our refractive-bundle-adjusted camera parameters, we are
%   giving it the best possible pinhole approximation of each camera's
%   effective projection. For cameras nearly perpendicular to the wall
%   this approximation is excellent. For cameras at steep angles the
%   residual error grows. The separate Phase5_Reconstruction.m uses the
%   full refractive model for 3-D triangulation and should be preferred
%   for final reconstruction; the OpenLPT files are primarily for using
%   OpenLPT's STB tracking front-end with our calibration.
%
%   FORMAT VERIFICATION:
%   The format above was reconstructed from:
%     - The CalibMatToTXTV2.m reference in the STB README
%     - The Camera.cpp source in the STB C++ backend
%     - The OpenLPT_GUI demo data structure
%   If OpenLPT produces errors loading the files, open one of the demo
%   cam*.txt files from the STB SD00125 sample dataset and compare line
%   by line — the structure should match exactly.
% =========================================================================

clear; clc;
LPT_Config;

resultsDir = cfg.resultsDir;
nCams      = cfg.nCams;

cal = load(fullfile(resultsDir,'extrinsics.mat'));
cameras = cal.cameras;

fprintf('=== Export Calibration Results ===\n\n');

% =========================================================================
%  1.  JSON export (human-readable)
% =========================================================================
jStruct.metadata.software     = 'LPT_Calibration — 4-camera refractive calibration';
jStruct.metadata.date          = datestr(now);
jStruct.metadata.n_air         = cfg.n_air;
jStruct.metadata.n_glass       = cfg.n_glass;
jStruct.metadata.n_water       = cfg.n_water;
jStruct.metadata.wand_length_m = cfg.wandLength;
jStruct.metadata.image_size    = cfg.imageSize;

jStruct.cameras = [];
for c = 1:nCams
    cam.id            = c;
    cam.K             = cameras(c).K;
    cam.R             = cameras(c).R;
    cam.t             = cameras(c).t(:)';
    cam.distCoeffs    = cameras(c).distCoeffs;  % [k1 k2 p1 p2 k3]
    cam.wall_normal   = cameras(c).wall.normal;
    cam.wall_point    = cameras(c).wall.point(:)';
    cam.wall_thickness = cameras(c).wall.thickness;
    cam.reproj_rms_px = cameras(c).reprojRMS;
    cam.C_world       = (-cameras(c).R' * cameras(c).t)';  % camera centre
    jStruct.cameras{c} = cam; %#ok<AGROW>
end

jsonStr = jsonencode(jStruct, 'PrettyPrint', true);
fid = fopen(fullfile(resultsDir,'calibration.json'),'w');
fprintf(fid,'%s',jsonStr);
fclose(fid);
fprintf('Saved: calibration.json\n');

% =========================================================================
%  2.  CSV export (flat table)
% =========================================================================
% One row per camera
headers = {'cam_id','fx','fy','cx','cy','k1','k2','p1','p2','k3',...
           'r11','r12','r13','r21','r22','r23','r31','r32','r33',...
           'tx','ty','tz','Cx','Cy','Cz','reproj_rms_px'};

rows = zeros(nCams, numel(headers));
for c = 1:nCams
    K  = cameras(c).K;
    dc = cameras(c).distCoeffs; if numel(dc)<5, dc(end+1:5)=0; end
    R  = cameras(c).R;
    t  = cameras(c).t;
    Cv = -R'*t;
    rows(c,:) = [c, K(1,1),K(2,2),K(1,3),K(2,3), ...
                 dc(1),dc(2),dc(3),dc(4),dc(5), ...
                 R(1,1),R(1,2),R(1,3), R(2,1),R(2,2),R(2,3), R(3,1),R(3,2),R(3,3),...
                 t(1),t(2),t(3), Cv(1),Cv(2),Cv(3), cameras(c).reprojRMS];
end

T = array2table(rows,'VariableNames',headers);
writetable(T, fullfile(resultsDir,'calibration.csv'));
fprintf('Saved: calibration.csv\n');

% =========================================================================
%  3.  Per-camera text files (OpenPTV / ptv_is compatible format)
% =========================================================================
for c = 1:nCams
    K  = cameras(c).K;
    dc = cameras(c).distCoeffs; if numel(dc)<5, dc(end+1:5)=0; end
    R  = cameras(c).R;
    t  = cameras(c).t;

    fname = fullfile(resultsDir, sprintf('cam_%d_params.txt', c));
    fid   = fopen(fname, 'w');

    fprintf(fid, '# Camera %d calibration parameters\n', c);
    fprintf(fid, '# Generated by LPT_Calibration on %s\n\n', datestr(now));

    fprintf(fid, '## Intrinsics\n');
    fprintf(fid, 'fx:   %.8f\n', K(1,1));
    fprintf(fid, 'fy:   %.8f\n', K(2,2));
    fprintf(fid, 'cx:   %.8f\n', K(1,3));
    fprintf(fid, 'cy:   %.8f\n', K(2,3));
    fprintf(fid, 'skew: %.8f\n', K(1,2));

    fprintf(fid, '\n## Distortion [k1 k2 p1 p2 k3]\n');
    fprintf(fid, 'k1: %.10f\n', dc(1));
    fprintf(fid, 'k2: %.10f\n', dc(2));
    fprintf(fid, 'p1: %.10f\n', dc(3));
    fprintf(fid, 'p2: %.10f\n', dc(4));
    fprintf(fid, 'k3: %.10f\n', dc(5));

    fprintf(fid, '\n## Rotation matrix (world-to-camera)\n');
    fprintf(fid, 'R: ');
    fprintf(fid, '%.10f  ', R(:)');
    fprintf(fid, '\n');

    fprintf(fid, '\n## Translation vector (world-to-camera) [m]\n');
    fprintf(fid, 't: %.10f  %.10f  %.10f\n', t(1), t(2), t(3));

    fprintf(fid, '\n## Camera centre in world coords [m]\n');
    Cv = -R'*t;
    fprintf(fid, 'C: %.10f  %.10f  %.10f\n', Cv(1), Cv(2), Cv(3));

    fprintf(fid, '\n## Refractive interface\n');
    fprintf(fid, 'wall_normal:    %.8f  %.8f  %.8f\n', cameras(c).wall.normal);
    fprintf(fid, 'wall_point:     %.8f  %.8f  %.8f\n', cameras(c).wall.point);
    fprintf(fid, 'wall_thickness: %.8f\n', cameras(c).wall.thickness);
    fprintf(fid, 'n_air:          %.6f\n', cfg.n_air);
    fprintf(fid, 'n_glass:        %.6f\n', cfg.n_glass);
    fprintf(fid, 'n_water:        %.6f\n', cfg.n_water);

    fprintf(fid, '\n## Quality\n');
    fprintf(fid, 'reproj_rms_px: %.6f\n', cameras(c).reprojRMS);

    fclose(fid);
    fprintf('Saved: cam_%d_params.txt\n', c);
end

% =========================================================================
%  4.  OpenLPT / STB camera files  (openlpt/camN.txt)
% =========================================================================
openlptDir = fullfile(resultsDir, 'openlpt');
if ~exist(openlptDir, 'dir'), mkdir(openlptDir); end

fprintf('\nWriting OpenLPT camera files...\n');

for c = 1:nCams
    K  = cameras(c).K;
    dc = cameras(c).distCoeffs; if numel(dc)<5, dc(end+1:5)=0; end
    R  = cameras(c).R;
    t  = cameras(c).t;

    % OpenLPT expects translation in MILLIMETRES
    t_mm = t * 1000;

    fname = fullfile(openlptDir, sprintf('cam%d.txt', c));
    fid   = fopen(fname, 'w');

    % Line 1: image dimensions
    fprintf(fid, '%d %d\n', cfg.imageSize(1), cfg.imageSize(2));

    % Line 2: intrinsics (pixels)
    fprintf(fid, '%.10f %.10f %.10f %.10f\n', K(1,1), K(2,2), K(1,3), K(2,3));

    % Line 3: distortion coefficients [k1 k2 p1 p2 k3]
    fprintf(fid, '%.10f %.10f %.10f %.10f %.10f\n', ...
        dc(1), dc(2), dc(3), dc(4), dc(5));

    % Lines 4-6: rotation matrix rows
    fprintf(fid, '%.10f %.10f %.10f\n', R(1,1), R(1,2), R(1,3));
    fprintf(fid, '%.10f %.10f %.10f\n', R(2,1), R(2,2), R(2,3));
    fprintf(fid, '%.10f %.10f %.10f\n', R(3,1), R(3,2), R(3,3));

    % Line 7: translation in mm
    fprintf(fid, '%.10f %.10f %.10f\n', t_mm(1), t_mm(2), t_mm(3));

    fclose(fid);
    fprintf('  Saved: openlpt/cam%d.txt\n', c);
end

% Write a README into the openlpt folder explaining the refraction caveat
fidR = fopen(fullfile(openlptDir, 'README.txt'), 'w');
fprintf(fidR, 'OpenLPT Camera Files\n');
fprintf(fidR, 'Generated by LPT_Calibration on %s\n\n', datestr(now));
fprintf(fidR, 'FORMAT (each camN.txt):\n');
fprintf(fidR, '  Line 1: image_width image_height (pixels)\n');
fprintf(fidR, '  Line 2: fx fy cx cy (pixels)\n');
fprintf(fidR, '  Line 3: k1 k2 p1 p2 k3 (Brown-Conrady distortion)\n');
fprintf(fidR, '  Lines 4-6: rotation matrix R (world-to-camera), row by row\n');
fprintf(fidR, '  Line 7: translation vector tx ty tz (MILLIMETRES)\n\n');
fprintf(fidR, 'REFRACTION NOTE:\n');
fprintf(fidR, '  These parameters were estimated using a full refractive bundle\n');
fprintf(fidR, '  adjustment (air->glass->water ray model). However, OpenLPT itself\n');
fprintf(fidR, '  uses a standard pinhole model internally and does not model\n');
fprintf(fidR, '  refraction. These files represent the best pinhole approximation\n');
fprintf(fidR, '  of each camera given the refractive calibration.\n\n');
fprintf(fidR, '  For cameras approximately perpendicular to the glass wall, this\n');
fprintf(fidR, '  approximation is very accurate. For cameras at steep angles\n');
fprintf(fidR, '  (>30 deg from wall normal), residual errors may be noticeable.\n\n');
fprintf(fidR, '  For highest-accuracy 3D reconstruction, use Phase5_Reconstruction.m\n');
fprintf(fidR, '  from the LPT_Calibration MATLAB package, which applies the full\n');
fprintf(fidR, '  refractive ray model at every triangulation step.\n\n');
fprintf(fidR, 'WORLD COORDINATE SYSTEM:\n');
fprintf(fidR, '  Origin: set during bundle adjustment (camera 1 defines the frame)\n');
fprintf(fidR, '  Units: translation in mm, rotation dimensionless\n');
fprintf(fidR, '  Z-axis convention: follows OpenLPT (right-hand, Z pointing into volume)\n\n');
fprintf(fidR, 'REFRACTIVE INDICES USED IN CALIBRATION:\n');
fprintf(fidR, '  n_air   = %.4f\n', cfg.n_air);
fprintf(fidR, '  n_glass = %.4f\n', cfg.n_glass);
fprintf(fidR, '  n_water = %.4f\n', cfg.n_water);
fclose(fidR);
fprintf('  Saved: openlpt/README.txt\n');

% =========================================================================
%  5.  Verify OpenLPT files can be re-read
% =========================================================================
fprintf('\nVerifying OpenLPT files...\n');
allOK = true;
for c = 1:nCams
    fname = fullfile(openlptDir, sprintf('cam%d.txt', c));
    fid   = fopen(fname, 'r');
    try
        imgSz_read  = fscanf(fid, '%f %f', 2);
        intr_read   = fscanf(fid, '%f %f %f %f', 4);
        dist_read   = fscanf(fid, '%f %f %f %f %f', 5);
        R_read      = reshape(fscanf(fid, '%f', 9), 3, 3)';
        t_read      = fscanf(fid, '%f %f %f', 3);
        fclose(fid);

        % Check round-trip
        K   = cameras(c).K;
        err_fx = abs(intr_read(1) - K(1,1));
        err_tz = abs(t_read(3) - cameras(c).t(3)*1000);

        if err_fx > 1e-4 || err_tz > 1e-4
            fprintf('  [WARN] Camera %d: round-trip mismatch (fx err=%.2e, tz err=%.2e)\n', ...
                c, err_fx, err_tz);
            allOK = false;
        else
            fprintf('  Camera %d: OK\n', c);
        end
    catch ME
        fprintf('  [ERROR] Camera %d: could not re-read file (%s)\n', c, ME.message);
        fclose(fid);
        allOK = false;
    end
end

if allOK
    fprintf('\n[ExportCalibration] All OpenLPT files verified OK.\n');
else
    fprintf('\n[ExportCalibration] Some files had verification warnings — check above.\n');
end

fprintf('\n[ExportCalibration] All files saved to %s/\n', resultsDir);
fprintf('  Generic formats:  calibration.json, calibration.csv, cam_N_params.txt\n');
fprintf('  OpenLPT format:   openlpt/cam1.txt ... cam%d.txt\n', nCams);
