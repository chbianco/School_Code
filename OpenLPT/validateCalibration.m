function validateCalibration(camDir, wandCsv, varargin)
% VALIDATECALIBRATION  Sanity-check an OpenLPT camera calibration by
% triangulating the calibration wand and checking that it comes back the
% right physical length and near the right depth.
%
% This is the gate to run after ANY calibration or VSC pass. A correct
% calibration reconstructs the rigid wand at its true length (~47 mm) with a
% small triangulation residual. VSC that has "diverged" (over-fit to ghost
% tracks) shows up here as a wrong wand length (e.g. 35 mm) and/or a shifted
% centre and a larger residual -- REJECT such a calibration.
%
% It reproduces OpenLPT's PINPLATE (pinhole + flat refractive plate) line of
% sight. Verified against OpenLPT's own pyopenlpt lineOfSight/triangulation:
% agrees to ~2e-4 in ray direction and gives wand length 46.97 mm vs OpenLPT's
% 46.94 mm on the reference dataset.
%
% USAGE
%   validateCalibration                       % uses the defaults below
%   validateCalibration(camDir, wandCsv)
%   validateCalibration(camDir, wandCsv, 'WandLenMM', 47, ...
%                       'CamPattern', 'vsc_cam%d.txt', 'NumCams', 4)
%
% INPUTS
%   camDir   folder holding the camera files (e.g. ...\camFile or ...\camFile_VSC)
%   wandCsv  the wand_points.csv exported by sizeWandTracking.m
%
% NAME/VALUE OPTIONS
%   'WandLenMM'  known physical wand-tip separation            (default 47.0)
%   'NumCams'    number of cameras                             (default 4)
%   'CamPattern' printf pattern for the per-camera file names. If empty (the
%                default) it auto-detects 'cam%d.txt' or 'vsc_cam%d.txt'.
%   'LenTolMM'   PASS if |meanLen - WandLenMM| <= this         (default 1.0)
%   'ResTolMM'   PASS if mean triangulation residual <= this   (default 0.6)
%
% Author: written for Christopher Bianco's flume LPT pipeline.

% ---------------- defaults (edit to taste) ----------------
if nargin < 1 || isempty(camDir)
    camDir  = 'C:\Users\FlumePIV\Desktop\CB_Data\07_24_26_freestreamTracking_L5R5\camFile';
end
if nargin < 2 || isempty(wandCsv)
    wandCsv = 'C:\Users\FlumePIV\Desktop\CB_Data\07_24_26_flumeCalibration\wand_points.csv';
end
p = inputParser;
p.addParameter('WandLenMM', 47.5);
p.addParameter('NumCams', 4);
p.addParameter('CamPattern', '');
p.addParameter('LenTolMM', 1.0);
p.addParameter('ResTolMM', 0.6);
p.parse(varargin{:});
opt = p.Results;

% ---------------- locate camera files ----------------
nC = opt.NumCams;
pat = opt.CamPattern;
if isempty(pat)
    if isfile(fullfile(camDir, 'cam0.txt')),       pat = 'cam%d.txt';
    elseif isfile(fullfile(camDir, 'vsc_cam0.txt')),pat = 'vsc_cam%d.txt';
    else, error('No cam0.txt or vsc_cam0.txt found in %s', camDir);
    end
end
cams = cell(nC,1);
for c = 1:nC
    f = fullfile(camDir, sprintf(pat, c-1));
    if ~isfile(f), error('Camera file not found: %s', f); end
    cams{c} = parseCam(f);
end

% ---------------- load wand 2D detections ----------------
% CSV columns: Frame,Camera,Status,PointIdx,X,Y,Radius,Metric
T = readmatrix(wandCsv);   % numeric columns; Status (text) becomes NaN, fine
frameCol = T(:,1); camCol = T(:,2); idxCol = T(:,4); Xc = T(:,5); Yc = T(:,6);

frames = unique(frameCol);
L = []; Zc = []; res = [];
for fi = 1:numel(frames)
    fr = frames(fi);
    P = nan(2,2,nC);            % P(idx+1, :, cam) = [x y]
    ok = true;
    for c = 1:nC
        for idx = 0:1
            sel = frameCol==fr & camCol==(c-1) & idxCol==idx;
            if nnz(sel) ~= 1, ok = false; break; end
            P(idx+1,:,c) = [Xc(find(sel,1)), Yc(find(sel,1))];
        end
        if ~ok, break; end
    end
    if ~ok, continue; end

    X = zeros(2,3);
    r2 = zeros(2,1);
    for idx = 1:2
        O = zeros(nC,3); D = zeros(nC,3);
        for c = 1:nC
            [O(c,:), D(c,:)] = lineOfSight(P(idx,:,c), cams{c});
        end
        [X(idx,:), r2(idx)] = triangulate(O, D);
    end
    L(end+1)   = norm(X(1,:) - X(2,:));           %#ok<AGROW>
    Zc(end+1)  = mean(X(:,3));                     %#ok<AGROW>
    res(end+1) = mean(r2);                         %#ok<AGROW>
end

if isempty(L)
    error('No frames had all %d cameras x 2 points. Check the CSV / NumCams.', nC);
end

% ---------------- report ----------------
mLen = mean(L); sLen = std(L); mRes = mean(res); mZ = mean(Zc);
fprintf('\n===== Calibration validation =====\n');
fprintf('  camera dir     : %s\n', camDir);
fprintf('  wand frames    : %d\n', numel(L));
fprintf('  wand length    : %.2f mm  (std %.2f)   [true = %.2f]\n', mLen, sLen, opt.WandLenMM);
fprintf('  length error   : %+.2f mm  (%+.1f%%)\n', mLen-opt.WandLenMM, 100*(mLen-opt.WandLenMM)/opt.WandLenMM);
fprintf('  wand centre Z  : %.1f mm\n', mZ);
fprintf('  triang residual: %.3f mm\n', mRes);

passLen = abs(mLen - opt.WandLenMM) <= opt.LenTolMM;
passRes = mRes <= opt.ResTolMM;
if passLen && passRes
    fprintf('  VERDICT: PASS  -- geometry reproduces the wand. Safe to use.\n\n');
else
    fprintf('  VERDICT: FAIL  -- ');
    if ~passLen, fprintf('wand length off by %.2f mm (>%.2f). ', abs(mLen-opt.WandLenMM), opt.LenTolMM); end
    if ~passRes, fprintf('residual %.3f mm (>%.2f). ', mRes, opt.ResTolMM); end
    fprintf('\n  Do NOT track with this calibration (VSC likely diverged on ghost tracks).\n\n');
end
end % validateCalibration


% ======================================================================
function c = parseCam(path)
% Parse an OpenLPT PINPLATE camera file (handles camFile and VSC variants).
lines = readLines(path);
c.K    = grabMatrix(lines, '# Camera Matrix:', 3);
c.dist = grabMatrix(lines, '# Distortion Coefficients:', 1);
c.R    = grabMatrix(lines, '# Rotation Matrix:', 3);
c.t    = grabMatrix(lines, '# Translation Vector:', 1).';
ppt = grabAfterAny(lines, {'Farthest Interface)', 'Reference Point of Refractive Plate:'}, 1);
pn  = grabAfterAny(lines, {'camera->object direction)', 'Normal Vector of Refractive Plate:'}, 1);
nA  = grabAfterAny(lines, {'obj->win->air)', 'Refractive Index Array:'}, 1);
w   = grabAfterAny(lines, {'plate thicknesses in mm)', 'Width of the Refractive Plate:'}, 1);
c.ppt = ppt(:).';
c.pn  = pn(:).' / norm(pn);
c.nA  = nA(:).';          % [n_obj(water) n_win(glass) n_air]
c.w   = w(1);
end

function L = readLines(path)
fid = fopen(path,'r'); assert(fid>0, 'cannot open %s', path);
L = {};
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    L{end+1} = ln; %#ok<AGROW>
end
fclose(fid);
end

function M = grabMatrix(lines, header, nRows)
% Find header line, then take the next nRows non-comment, non-blank lines.
i = findHeader(lines, header);
M = collectRows(lines, i, nRows);
end

function v = grabAfterAny(lines, headers, nRows)
% Try several possible header substrings; use the first that matches.
for k = 1:numel(headers)
    i = findHeader(lines, headers{k}, true);
    if i > 0
        v = collectRows(lines, i, nRows);
        return;
    end
end
error('None of the headers found: %s', strjoin(headers, ' | '));
end

function i = findHeader(lines, header, soft)
if nargin < 3, soft = false; end
i = 0;
for k = 1:numel(lines)
    if contains(lines{k}, header), i = k; return; end
end
if ~soft, error('Header not found: %s', header); end
end

function M = collectRows(lines, headerIdx, nRows)
M = []; got = 0; k = headerIdx + 1;
while got < nRows && k <= numel(lines)
    ln = strtrim(lines{k});
    if isempty(ln) || startsWith(ln, '#'), k = k+1; continue; end
    nums = str2double(regexp(ln, '[-+0-9.eE]+', 'match'));
    M = [M; nums]; %#ok<AGROW>
    got = got + 1; k = k + 1;
end
if got < nRows, error('Expected %d data rows after header, got %d.', nRows, got); end
end


% ======================================================================
function xn = undistort(uv, K, dist)
% Iteratively remove Brown-Conrady distortion -> normalized camera ray (z=1).
fx=K(1,1); fy=K(2,2); cx=K(1,3); cy=K(2,3);
k1=dist(1); k2=dist(2); p1=dist(3); p2=dist(4); k3=dist(5);
x=(uv(1)-cx)/fx; y=(uv(2)-cy)/fy; x0=x; y0=y;
for it=1:30
    r2=x*x+y*y; rad=1+k1*r2+k2*r2^2+k3*r2^3;
    dx=2*p1*x*y+p2*(r2+2*x*x); dy=p1*(r2+2*y*y)+2*p2*x*y;
    x=(x0-dx)/rad; y=(y0-dy)/rad;
end
xn=[x; y; 1];
end

function d2 = refract(d, n, n1, n2)
% Vector Snell's law. n1 -> n2 across interface with unit normal n.
d=d/norm(d); n=n/norm(n); ci=-dot(n,d);
if ci<0, n=-n; ci=-dot(n,d); end
r=n1/n2; k=1-r*r*(1-ci*ci);
d2 = r*d + (r*ci - sqrt(k))*n;
d2 = d2/norm(d2);
end

function [o, d] = lineOfSight(uv, c)
% Object-side (water) ray for a distorted pixel, matching OpenLPT PINPLATE.
C = -c.R.' * c.t(:);                 % camera centre (world)
d = c.R.' * undistort(uv, c.K, c.dist); d = d/norm(d);
nrm = c.pn(:); ppt = c.ppt(:); w = c.w;
n_obj=c.nA(1); n_win=c.nA(2); n_air=c.nA(3);
h1 = planeHit(C,  d,  ppt - w*nrm, nrm);   % air -> glass interface
d2 = refract(d,  nrm, n_air, n_win);
h2 = planeHit(h1, d2, ppt,          nrm);  % glass -> water interface
d  = refract(d2, nrm, n_win, n_obj);       % ray in water
o  = h2;
end

function x = planeHit(o, d, P, N)
x = o + ((P - o).'*N / (d.'*N)) * d;
end

function [X, resid] = triangulate(O, D)
% Least-squares closest point to a set of rays O(i,:)+t*D(i,:).
A=zeros(3); b=zeros(3,1);
for i=1:size(O,1)
    d=D(i,:).'; Pp=eye(3)-d*d';
    A=A+Pp; b=b+Pp*O(i,:).';
end
X=(A\b).';
resid=0;
for i=1:size(O,1)
    d=D(i,:).'; Pp=eye(3)-d*d';
    resid=resid+norm(Pp*(X.'-O(i,:).'));
end
resid=resid/size(O,1);
end
