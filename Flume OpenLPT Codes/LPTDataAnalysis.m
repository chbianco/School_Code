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