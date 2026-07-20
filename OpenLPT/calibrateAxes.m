function calibrateAxes()
% CALIBRATE_AXES  Interactive coordinate-axis calibration from the black-sphere
% "axis wand" images (cam0..cam3_axes.tif) for particle tracking.
%
% The wand
% --------
% Three rods cross at a single point (the coordinate-system ORIGIN). Each rod has
% a black ball at both ends, so there are 3 lines (one per axis) defined by pairs
% of balls. Ball size encodes depth (closer to the camera => larger). The wand's
% orientation differs in every camera, which is why you pick the positive
% directions by hand for each image.
%
% Workflow (per camera)
% ---------------------
% For each axis you click the POSITIVE-direction ball, then the NEGATIVE-direction
% ball. Each click is snapped to the centroid of the dark ball you clicked on (so
% rough clicks are fine). The origin is computed as the least-squares intersection
% of the axis lines, and one row is written to the output CSV:
%
%   cam_id,center_x,center_y,plus_x_x,plus_x_y,plus_y_x,plus_y_y,plus_z_x,plus_z_y
%
% where center_* is the origin and plus_<axis>_* is the positive-direction ball.
%
% Controls (the figure must have focus)
% -------------------------------------
%   left click   select the ball for the current prompt (snaps to centroid)
%   n            mark the current NEGATIVE ball as "not visible" (skip it)
%   u            undo the last point
%   r            restart the current image from scratch
%   enter        accept the result (only once all required points are set)
%   esc          skip this camera (writes no row)
%
% Uses only base MATLAB (no toolboxes required).
%
% Run:  calibrate_axes

% ----------------------------- Configuration ------------------------------ %
folder = fileparts(mfilename('fullpath'));
if isempty(folder), folder = pwd; end
cameras = { 0, 'cam0_axes.tif';
            1, 'cam1_axes.tif';
            2, 'cam2_axes.tif';
            3, 'cam3_axes.tif' };
outputCsv   = fullfile(folder, 'axis_directions.csv');
searchHalf  = 75;      % half-window (px) searched around each click for the ball
saveOverlays = true;   % save an annotated PNG next to each source image

axisColors.x = [0.90 0.10 0.29];   % red   -> X
axisColors.y = [0.24 0.71 0.29];   % green -> Y
axisColors.z = [0.26 0.39 0.85];   % blue  -> Z

fprintf('%s\n', repmat('-', 1, 70));
fprintf(['Axis calibration.  For each camera click +X,-X, +Y,-Y, +Z,-Z balls.\n' ...
         'Keys: n=skip negative  u=undo  r=restart  enter=accept  esc=skip.\n']);
fprintf('%s\n', repmat('-', 1, 70));

rows = [];   % each row: [cam_id cx cy px_x px_y py_x py_y pz_x pz_y]

for k = 1:size(cameras, 1)
    camId = cameras{k, 1};
    fname = cameras{k, 2};
    imgPath = fullfile(folder, fname);
    if ~isfile(imgPath)
        fprintf('[skip] cam%d: %s not found\n', camId, fname);
        continue
    end
    fprintf('\n=== cam%d  (%s) ===\n', camId, fname);
    gray = readGray(imgPath);
    row  = pickAxes(gray, camId, searchHalf, axisColors, saveOverlays, folder);
    if isempty(row)
        fprintf('cam%d: skipped (no row written)\n', camId);
        continue
    end
    rows = [rows; row]; %#ok<AGROW>
    fprintf(['cam%d: origin=(%.1f, %.1f)  +X=(%.1f, %.1f)  ' ...
             '+Y=(%.1f, %.1f)  +Z=(%.1f, %.1f)\n'], camId, ...
             row(2), row(3), row(4), row(5), row(6), row(7), row(8), row(9));
end

if isempty(rows)
    fprintf('\nNo rows collected; nothing written.\n');
    return
end
writeCsv(outputCsv, rows);
fprintf('\nWrote %d row(s) to %s\n', size(rows, 1), outputCsv);
end


% ======================================================================== %
%                         Interactive picker                               %
% ======================================================================== %
function row = pickAxes(gray, camId, half, colors, saveOverlays, folder)
row = [];
stepAxis = {'x','x','y','y','z','z'};
stepSign = [ 1  -1   1  -1   1  -1 ];
stepReq  = [ 1   0   1   0   1   0 ];   % positive directions are required
nSteps   = numel(stepAxis);
[imgH, imgW] = size(gray);

% points(i,:) = [x y radius]; state: 0 unset, 1 set, 2 skipped
pts   = nan(nSteps, 3);
state = zeros(nSteps, 1);
step  = 1;

fig = figure('Name', sprintf('cam%d axes', camId), 'Color', 'w', ...
             'NumberTitle', 'off');
ax  = axes('Parent', fig);
set(fig, 'WindowButtonDownFcn', @onClick, 'WindowKeyPressFcn', @onKey);

redraw();
uiwait(fig);          % blocks until a callback deletes the figure
return                % 'row' is set (accept) or left empty (skip)

    % ------------------------------ click -------------------------------- %
    function onClick(~, ~)
        if ~isvalid(ax) || ~strcmp(get(fig, 'SelectionType'), 'normal')
            return    % only left-clicks select a ball
        end
        cp = get(ax, 'CurrentPoint');
        x = cp(1, 1); y = cp(1, 2);
        if x < 1 || x > imgW || y < 1 || y > imgH
            return    % click landed outside the image
        end
        if step > nSteps
            return
        end
        found = refineToBall(gray, x, y, half);
        if isempty(found)
            title(ax, 'No dark ball found there - click closer to a ball.');
            drawnow; return
        end
        pts(step, :) = found;
        state(step)  = 1;
        step = step + 1;
        redraw();
    end

    % ------------------------------ keys --------------------------------- %
    function onKey(~, ev)
        switch ev.Key
            case 'n'
                if step <= nSteps && ~stepReq(step)
                    state(step) = 2;   % negative ball not visible
                    step = step + 1;
                    redraw();
                end
            case 'u'
                [step, pts, state] = undo(step, pts, state);
                redraw();
            case 'r'
                pts(:) = NaN; state(:) = 0; step = 1;
                redraw();
            case 'escape'
                delete(fig);           % skip camera (row stays empty)
            case 'return'
                if step > nSteps
                    [ok, origin, plusPts] = solve(pts, state, stepAxis, stepSign);
                    if ok
                        row = [camId, origin(1), origin(2), ...
                               plusPts.x(1), plusPts.x(2), ...
                               plusPts.y(1), plusPts.y(2), ...
                               plusPts.z(1), plusPts.z(2)];
                        if saveOverlays
                            saveOverlay(fig, folder, camId);
                        end
                        delete(fig);   % accept -> resume uiwait
                    end
                end
        end
    end

    % --------------------------- nested drawing -------------------------- %
    function redraw()
        cla(ax);
        imagesc(ax, gray);
        colormap(ax, repmat((0:255)'/255, 1, 3));   % grayscale (var 'gray' shadows the fcn)
        set(ax, 'CLim', [0 255], 'YDir', 'reverse');
        axis(ax, 'image'); axis(ax, 'off'); hold(ax, 'on');

        for i = 1:nSteps
            if state(i) ~= 1, continue, end
            col = colors.(stepAxis{i});
            lw  = 2.5; if stepSign(i) < 0, lw = 1.25; end   % negatives drawn thinner
            cx = pts(i,1); cy = pts(i,2); r = pts(i,3);
            rectangle(ax, 'Position', [cx-r, cy-r, 2*r, 2*r], ...
                      'Curvature', [1 1], 'EdgeColor', col, 'LineWidth', lw);
            plot(ax, cx, cy, '+', 'Color', col, 'MarkerSize', 12, 'LineWidth', 2);
            text(ax, cx + r + 5, cy, sprintf('%c%s', signChar(stepSign(i)), ...
                 upper(stepAxis{i})), 'Color', col, 'FontWeight', 'bold', ...
                 'FontSize', 11, 'Clipping', 'on');
        end

        if step > nSteps
            [ok, o, pp] = solve(pts, state, stepAxis, stepSign);
            if ok
                axesList = {'x','y','z'};
                for j = 1:3
                    ptp = pp.(axesList{j});
                    if any(isnan(ptp)), continue, end
                    plot(ax, [o(1) ptp(1)], [o(2) ptp(2)], '-', ...
                         'Color', colors.(axesList{j}), 'LineWidth', 2.5);
                end
                plot(ax, o(1), o(2), 'o', 'MarkerFaceColor', 'y', ...
                     'MarkerEdgeColor', 'k', 'MarkerSize', 10);
                text(ax, o(1)+8, o(2)-8, 'origin', 'Color', [0.9 0.9 0], ...
                     'FontWeight', 'bold', 'FontSize', 11);
                title(ax, sprintf(['cam%d: ENTER = accept   ' ...
                      'u = undo   r = restart   esc = skip'], camId));
            else
                title(ax, ['Need >= 2 complete axis lines to locate the ' ...
                           'origin. Press r to restart.']);
            end
        else
            title(ax, promptText(camId, stepAxis{step}, stepSign(step)));
        end
        hold(ax, 'off');
        drawnow;
    end
end


function s = promptText(camId, axisChar, sgn)
if sgn > 0
    s = sprintf(['cam%d: click the  +%s  ball (positive %s direction)   ' ...
                 '[u undo | r restart | esc skip]'], camId, upper(axisChar), ...
                 upper(axisChar));
else
    s = sprintf(['cam%d: click the  -%s  ball, or press  n  if not visible   ' ...
                 '[u undo | r restart | esc skip]'], camId, upper(axisChar));
end
end


function c = signChar(sgn)
if sgn > 0, c = '+'; else, c = '-'; end
end


function [step, pts, state] = undo(step, pts, state)
if step > 1
    step = step - 1;
    pts(step, :) = NaN;
    state(step)  = 0;
end
end


function [ok, origin, plusPts] = solve(pts, state, stepAxis, stepSign)
% Least-squares origin from complete axis lines; positive-ball positions.
ok = false; origin = [NaN NaN];
plusPts = struct('x', [NaN NaN], 'y', [NaN NaN], 'z', [NaN NaN]);
axesList = {'x','y','z'};
pairs = {};
for j = 1:numel(axesList)
    ax = axesList{j};
    ip = find(strcmp(stepAxis, ax) & stepSign > 0);
    im = find(strcmp(stepAxis, ax) & stepSign < 0);
    if state(ip) == 1
        plusPts.(ax) = pts(ip, 1:2);
    else
        return   % every positive direction is required
    end
    if state(im) == 1
        pairs{end+1} = [pts(ip,1:2); pts(im,1:2)]; %#ok<AGROW>
    end
end
if numel(pairs) < 2
    return
end
origin = lsIntersection(pairs);
ok = true;
end


function o = lsIntersection(pairs)
% Best-fit intersection of lines. Each pair is [p_plus; p_minus] (2x2).
A = zeros(2); b = zeros(2, 1);
for i = 1:numel(pairs)
    p = pairs{i}(1, :)';
    d = pairs{i}(2, :)' - p;
    n = norm(d);
    if n == 0, continue, end
    d = d / n;
    P = eye(2) - d*d';        % projects onto the line's normal direction
    A = A + P;
    b = b + P*p;
end
o = (A \ b)';
end


% ======================================================================== %
%                         Image / blob helpers                             %
% ======================================================================== %
function gray = readGray(imgPath)
img = imread(imgPath);
if ndims(img) == 3
    img = mean(img, 3);       % avoids the rgb2gray toolbox dependency
end
gray = double(img);
end


function out = refineToBall(gray, x, y, half)
% Snap a rough click to the centroid of the dark ball under/near it.
% Returns [cx cy radius] in full-image pixel coords, or [] if none found.
out = [];
[H, W] = size(gray);
xc = round(x); yc = round(y);
x0 = max(1, xc - half); x1 = min(W, xc + half);
y0 = max(1, yc - half); y1 = min(H, yc + half);
if x0 >= x1 || y0 >= y1, return, end
win = gray(y0:y1, x0:x1);

T = otsuThreshold(win);
dark = win < T;
if ~any(dark(:)), return, end

% seed = clicked pixel in window coords; if not dark, jump to nearest dark px
sr = yc - y0 + 1; sc = xc - x0 + 1;
sr = min(max(sr, 1), size(win, 1));
sc = min(max(sc, 1), size(win, 2));
if ~dark(sr, sc)
    [dr, dc] = find(dark);
    [~, mi] = min((dr - sr).^2 + (dc - sc).^2);
    sr = dr(mi); sc = dc(mi);
end

[rr, cc] = floodFill(dark, sr, sc);
area = numel(rr);
cx = mean(cc) + x0 - 1;
cy = mean(rr) + y0 - 1;
radius = sqrt(area / pi);
out = [cx, cy, radius];
end


function T = otsuThreshold(win)
% Otsu threshold (0..255). The dark ball is the low-intensity class.
counts = histcounts(win(:), 0:256);   % 256 bins covering [0,255]
counts = double(counts(:));
total = sum(counts);
if total == 0, T = 128; return, end
p = counts / total;
levels = (0:255)';
omega = cumsum(p);
mu = cumsum(p .* levels);
muT = mu(end);
denom = omega .* (1 - omega);
sigmaB = (muT * omega - mu).^2 ./ denom;
sigmaB(denom == 0) = -Inf;
[~, idx] = max(sigmaB);
T = idx - 1;   % bin index -> intensity value
end


function [rows, cols] = floodFill(mask, r0, c0)
% 8-connected flood fill from (r0,c0). Returns row/col indices of the region.
[H, W] = size(mask);
visited = false(H, W);
cap = H * W;
stack = zeros(cap, 2);
region = zeros(cap, 2);
sp = 1; stack(1, :) = [r0, c0]; visited(r0, c0) = true;
rp = 0;
nbr = [-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1];
while sp > 0
    rc = stack(sp, :); sp = sp - 1;
    rp = rp + 1; region(rp, :) = rc;
    for k = 1:8
        nr = rc(1) + nbr(k, 1);
        nc = rc(2) + nbr(k, 2);
        if nr >= 1 && nr <= H && nc >= 1 && nc <= W ...
                && ~visited(nr, nc) && mask(nr, nc)
            visited(nr, nc) = true;
            sp = sp + 1; stack(sp, :) = [nr, nc];
        end
    end
end
region = region(1:rp, :);
rows = region(:, 1);
cols = region(:, 2);
end


% ======================================================================== %
%                              Output                                      %
% ======================================================================== %
function saveOverlay(fig, folder, camId)
out = fullfile(folder, sprintf('cam%d_axes_overlay.png', camId));
try
    exportgraphics(fig, out, 'Resolution', 120);   % R2020a+
catch
    saveas(fig, out);
end
end


function writeCsv(outputCsv, rows)
fid = fopen(outputCsv, 'w');
if fid == -1
    error('Could not open %s for writing.', outputCsv);
end
fprintf(fid, ['cam_id,center_x,center_y,plus_x_x,plus_x_y,' ...
              'plus_y_x,plus_y_y,plus_z_x,plus_z_y\n']);
for i = 1:size(rows, 1)
    r = rows(i, :);
    fprintf(fid, '%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n', ...
            r(1), r(2), r(3), r(4), r(5), r(6), r(7), r(8), r(9));
end
fclose(fid);
end
