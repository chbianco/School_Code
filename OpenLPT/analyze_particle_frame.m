function analyze_particle_frame(imagePath)
% ANALYZE_PARTICLE_FRAME  Frame analysis tool for Lagrangian particle tracking
%
% USAGE:
%   analyze_particle_frame()              % opens a file picker dialog
%   analyze_particle_frame('frame.tif')   % load a specific file

% ── Load image ───────────────────────────────────────────────────────────────
if nargin < 1 || isempty(imagePath)
    [fname, fpath] = uigetfile( ...
        {'*.tif;*.tiff;*.png;*.bmp;*.jpg;*.jpeg;*.pgm','Image Files';
         '*.*','All Files'}, ...
        'Select a particle tracking frame');
    if isequal(fname, 0), disp('No file selected.'); return; end
    imagePath = fullfile(fpath, fname);
end

info     = imfinfo(imagePath);
imgRaw   = imread(imagePath);
bitDepth = info(1).BitDepth;
maxVal   = 2^bitDepth - 1;

if size(imgRaw, 3) == 3
    imgGray = rgb2gray(imgRaw);
else
    imgGray = imgRaw;
end
[nRows, nCols] = size(imgGray);
imgVec = double(imgGray(:));
imgVec = imgVec(imgVec > 0);   % exclude zero-intensity pixels

% ── Statistics ───────────────────────────────────────────────────────────────
fprintf('\n--- Image: %s ---\n', imagePath);
fprintf('  Resolution : %d x %d px,  %d-bit\n', nCols, nRows, bitDepth);
fprintf('  Min/Max    : %.0f / %.0f\n',   min(imgVec), max(imgVec));
fprintf('  Mean +- Std : %.2f +- %.2f\n',  mean(imgVec), std(imgVec));
fprintf('  Saturated  : %.3f %%\n',        100*mean(imgVec >= maxVal));

% ════════════════════════════════════════════════════════════════════════════
%  A) INTENSITY HISTOGRAM
% ════════════════════════════════════════════════════════════════════════════
nBins = min(256, maxVal + 1);

figure('Name', 'Intensity Histogram');
histogram(imgVec, nBins, 'BinLimits', [0 maxVal], ...
          'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none');
xlabel('Pixel Intensity');
ylabel('Pixel Count');
title('Intensity Histogram');
xline(mean(imgVec),        'r-',  sprintf('Mean = %.1f',   mean(imgVec)),        'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(median(imgVec),      'k--', sprintf('Median = %.1f', median(imgVec)),      'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(prctile(imgVec, 95), 'g--', sprintf('p95 = %.0f',    prctile(imgVec, 95)), 'LineWidth', 1.2, 'LabelVerticalAlignment', 'bottom');
grid on;
xlim([0 maxVal]);

% ════════════════════════════════════════════════════════════════════════════
%  B) CUMULATIVE PIXEL COUNT
% ════════════════════════════════════════════════════════════════════════════
edges    = linspace(0, maxVal, nBins + 1);
counts   = histcounts(imgVec, edges);
bCenters = 0.5 * (edges(1:end-1) + edges(2:end));
cumPct   = 100 * cumsum(counts) / numel(imgVec);

figure('Name', 'Cumulative Pixel Count');
plot(bCenters, cumPct, 'b-', 'LineWidth', 1.8);
xlabel('Pixel Intensity');
ylabel('Cumulative Pixel Count  (%)');
title('Cumulative Intensity Distribution');
xline(mean(imgVec),        'r-',  sprintf('Mean = %.1f',   mean(imgVec)),        'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(prctile(imgVec, 50), 'k--', sprintf('Median = %.1f', prctile(imgVec, 50)), 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(prctile(imgVec, 95), 'g--', sprintf('p95 = %.0f',    prctile(imgVec, 95)), 'LineWidth', 1.2, 'LabelVerticalAlignment', 'bottom');
grid on;
xlim([0 maxVal]);
ylim([0 100]);

% ════════════════════════════════════════════════════════════════════════════
%  C) INTERACTIVE IMAGE VIEWER WITH PIXEL RULER
% ════════════════════════════════════════════════════════════════════════════
figImg = figure('Name', 'Particle Frame — Pixel Ruler');

axImg = axes('Parent', figImg);
if bitDepth > 8
    imshow(mat2gray(imgRaw), 'Parent', axImg);
else
    imshow(imgRaw, 'Parent', axImg);
end
axis(axImg, 'on');
set(axImg, 'YDir', 'reverse');
xlabel(axImg, 'x  (pixels)');
ylabel(axImg, 'y  (pixels)');
title(axImg, 'Click two points to measure  |  Press ''c'' to clear');

infoTxt = uicontrol('Parent', figImg, 'Style', 'text', ...
    'Units', 'normalized', 'Position', [0.05 0.01 0.90 0.04], ...
    'String', 'Click a first point to begin measuring.', ...
    'FontSize', 10, 'HorizontalAlignment', 'center', ...
    'BackgroundColor', get(figImg, 'Color'));

% ── Ruler state ──────────────────────────────────────────────────────────────
state.pt1     = [];
state.waiting = true;
state.handles = {};

set(figImg, 'WindowButtonDownFcn', @onMouseClick);
set(figImg, 'KeyPressFcn',         @onKeyPress);

% ── Callbacks ────────────────────────────────────────────────────────────────
    function onMouseClick(~, ~)
        if gca ~= axImg, return; end
        cp = get(axImg, 'CurrentPoint');
        x = cp(1,1);  y = cp(1,2);
        if x < 0.5 || x > nCols+0.5 || y < 0.5 || y > nRows+0.5, return; end

        if state.waiting
            state.pt1 = [x, y];
            h = plot(axImg, x, y, '+r', 'MarkerSize', 14, 'LineWidth', 2);
            state.handles{end+1} = h;
            set(infoTxt, 'String', sprintf('Pt 1: (%.1f, %.1f)  —  Now click the second point.', x, y));
            state.waiting = false;
        else
            x1 = state.pt1(1);  y1 = state.pt1(2);
            dist = sqrt((x-x1)^2 + (y-y1)^2);

            hl  = plot(axImg, [x1 x], [y1 y], 'r-', 'LineWidth', 1.8);
            hm1 = plot(axImg, x1, y1, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
            hm2 = plot(axImg, x,  y,  'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
            ht  = text(axImg, (x1+x)/2, (y1+y)/2, sprintf(' %.1f px', dist), ...
                       'Color', 'red', 'FontSize', 9, 'FontWeight', 'bold');

            state.handles(end+1:end+4) = {hl, hm1, hm2, ht};

            fprintf('  Measurement: %.2f px  (Dx=%.2f, Dy=%.2f)\n', dist, abs(x-x1), abs(y-y1));
            set(infoTxt, 'String', sprintf('Distance: %.2f px   (Dx=%.1f, Dy=%.1f)   |  Click a new first point or press c to clear.', ...
                dist, abs(x-x1), abs(y-y1)));

            state.pt1     = [];
            state.waiting = true;
        end
    end

    function onKeyPress(~, evt)
        if strcmpi(evt.Key, 'c')
            for i = 1:numel(state.handles)
                if ishandle(state.handles{i}), delete(state.handles{i}); end
            end
            state.handles = {};
            state.pt1     = [];
            state.waiting = true;
            set(infoTxt, 'String', 'Cleared.  Click a first point to begin measuring.');
        end
    end

end