% =========================================================================
% PIVDataAnalysis.m
% Author: Christopher Bianco
% christopherhbianco@gmail.com
% =========================================================================
%
% DEPENDENCIES:
% - None
%
% PURPOSE:
% Analyze planar PIV data exported from PIVlab. Computes vertical profiles
% of mean velocity and turbulence intensity for water-channel PIV.
%
% Expects a .mat file formatted the same way as a PIVlab session export.
% PIVlab stores every calculated field as a cell array with one cell per
% frame; each cell holds a [ny x nx] grid. Relevant variables:
%   x, y                     : grid coordinates (same for every frame) [m]
%   u_filtered, v_filtered   : velocities incl. data validation        [m/s]
%   typevector_filtered      : 0 masked, 1 regular vector, 2 filtered vector
% Rows index the vertical coordinate y, columns index the streamwise
% coordinate x. Positions are in m and velocities in m/s (PIVlab
% calibration already applied), so every downstream quantity is SI:
% velocities [m/s], turbulence intensity [-].

%% ======================= Preamble =======================
close all; clc;
clearvars
set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);

%% ======================= Load File =======================
[file, path] = uigetfile('*.*', 'Select a File');

% Check if the user cancelled the dialog
if isequal(file, 0)
    disp('Aborting');
else
    % Create the full path to the file
    fullPath = fullfile(path, file);
    fprintf('Selected: %s\n', fullPath);
end
file = fullPath; %PIVlab export to analyze.

%% ======================= Load Frames =======================

raw = load(file);

% PIVlab stores each field as a [nFrames x 1] cell array. The grid is
% identical for every frame, so the coordinate vectors are taken from the
% first frame.
xGrid = raw.x{1};          % [ny x nx] streamwise coordinate [m]
yGrid = raw.y{1};          % [ny x nx] vertical coordinate   [m]

% Pick the most-processed velocity field that is actually populated:
% smoothed > filtered > original. A field is only accepted if its cells
% exist and hold finite data (PIVlab leaves cells empty / all-NaN for
% steps that were never run, e.g. u_smoothed in an export without
% smoothing). typevectorName tracks the matching mask field, if present.
fieldOrder = {'u_smoothed', 'v_smoothed', 'typevector_filtered'; ...
              'u_filtered', 'v_filtered', 'typevector_filtered'; ...
              'u_original', 'v_original', 'typevector_original'};

uName = ''; vName = ''; tvName = '';
for r = 1:size(fieldOrder, 1)
    uf = fieldOrder{r, 1};
    vf = fieldOrder{r, 2};
    if isfield(raw, uf) && isfield(raw, vf) && ~isempty(raw.(uf)) && ...
       any(cellfun(@(c) ~isempty(c) && any(isfinite(c(:))), raw.(uf)))
        uName = uf; vName = vf; tvName = fieldOrder{r, 3};
        break
    end
end
if isempty(uName)
    error(['No populated velocity field found in %s ', ...
           '(checked smoothed/filtered/original).'], file);
end
% Only use the mask field if it actually exists in this export.
if ~isfield(raw, tvName), tvName = ''; end
fprintf('Using velocity field: %s / %s\n', uName, vName);

uCells = raw.(uName);
vCells = raw.(vName);

nFrames = numel(uCells);
[ny, nx] = size(xGrid);

% Scatter the per-frame grids into [ny x nx x nFrames] stacks so the
% averaging below reads as plain array reductions over the frame
% dimension. Vectors flagged as masked (typevector == 0) are set to NaN
% so they drop out of every 'omitnan' average.
u = NaN(ny, nx, nFrames);
v = NaN(ny, nx, nFrames);

for f = 1:nFrames
    uf = uCells{f};
    vf = vCells{f};

    % Drop masked vectors (typevector == 0) when a mask field is present.
    if ~isempty(tvName)
        tv = raw.(tvName){f};
        uf(tv == 0) = NaN;
        vf(tv == 0) = NaN;
    end

    u(:,:,f) = uf;
    v(:,:,f) = vf;
end

% Vertical coordinate for the profiles: one value per row (horizontal
% plane). Columns share a common y within a row, so the first column is
% representative.
yc = yGrid(:,1);   % ny x 1 [m]

fprintf('Loaded: %d frames x (%d rows y x %d cols x) grid\n', ...
    nFrames, ny, nx);

%% ======================= User Parameters =======================
% Component(s) plotted for the mean-velocity profile: both the streamwise
% (u) and vertical (v) means are shown vs height y.

% Turbulence intensity is the local RMS fluctuation normalized by the
% local mean streamwise velocity, TI = u_rms / U. Both the streamwise and
% vertical RMS are normalized by U (standard single-reference convention).

%% ======================= Time (Ensemble) Mean =======================
% Mean over all frames at every grid point -> time-mean field [ny x nx].
Umean = mean(u, 3, 'omitnan');   % ny x nx [m/s]
Vmean = mean(v, 3, 'omitnan');

% Fluctuations relative to the local time mean at each grid point.
uprime = u - Umean;
vprime = v - Vmean;

% RMS of the fluctuations over the frame dimension -> [ny x nx] [m/s].
urms = sqrt(mean(uprime.^2, 3, 'omitnan'));
vrms = sqrt(mean(vprime.^2, 3, 'omitnan'));

%% ======================= Horizontal-Plane Averaging =======================
% Average along the streamwise direction (columns, dim 2) within each
% horizontal plane to collapse the field to 1D profiles in y.
U_profile = squeeze(mean(Umean, 2, 'omitnan'));   % ny x 1 [m/s]
V_profile = squeeze(mean(Vmean, 2, 'omitnan'));

urms_profile = squeeze(mean(urms, 2, 'omitnan')); % ny x 1 [m/s]
vrms_profile = squeeze(mean(vrms, 2, 'omitnan'));

% Turbulence intensity profiles, normalized by the local mean streamwise
% velocity U(y).
TI_u = urms_profile ./ U_profile;   % ny x 1 [-]
TI_v = vrms_profile ./ U_profile;

fprintf('Bulk mean streamwise velocity: %.4g m/s\n', mean(U_profile, 'omitnan'));
fprintf('Bulk streamwise TI           : %.4g %%\n', 100*mean(TI_u, 'omitnan'));

%% ======================= Plotting =======================
%% Mean velocity profiles
figure;
hold on;
plot(U_profile, yc);
plot(V_profile, yc);
xlabel('$\langle \overline{u_i} \rangle$ [m/s]'); ylabel('$y$ [m]');
legend({'$\langle \overline{u} \rangle$', '$\langle \overline{v} \rangle$'}, ...
       'Location', 'best');
title('Mean velocity profiles');
grid on;
hold off;

%% Turbulence intensity profiles
figure;
hold on;
plot(100*TI_u, yc);
plot(100*TI_v, yc);
xlabel('Turbulence intensity [\%]'); ylabel('$y$ [m]');
legend({'$u_{\mathrm{rms}}/\langle \overline{u} \rangle$', ...
        '$v_{\mathrm{rms}}/\langle \overline{u} \rangle$'}, ...
       'Location', 'best');
title('Turbulence intensity profiles');
grid on;
hold off;
