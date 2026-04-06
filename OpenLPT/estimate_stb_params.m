function params = estimate_stb_params(imagePath)
% ESTIMATE_STB_PARAMS  Interactively estimate STB tracking parameters from a
%                      single sample frame.
%
% USAGE:
%   params = estimate_stb_params()
%   params = estimate_stb_params('frame.tif')
%
% WORKFLOW (follow the prompts in order):
%   Step 1 — Inspect the histogram, then enter an intensity threshold
%   Step 2 — Detected particles are overlaid; adjust threshold if needed
%   Step 3 — Measure a particle diameter with the ruler (two clicks)
%   Step 4 — Measure several nearest-neighbour spacings with the ruler
%   Step 5 — All parameters are computed and printed as a struct
%
% OUTPUT:
%   params  — struct with fields ready to pass to your STB pipeline:
%               .intensityThreshold
%               .tracerRadius
%               .searchRadius
%               .minFrameCount
%               .particleSpacing
%               .gridResolution
%               .shakeAmplitude
%               .ghostThreshold

% ── 0. Load image ─────────────────────────────────────────────────────────────
if nargin < 1 || isempty(imagePath)
    [fname, fpath] = uigetfile( ...
        {'*.tif;*.tiff;*.png;*.bmp;*.jpg;*.jpeg;*.pgm','Image Files';
         '*.*','All Files'}, ...
        'Select a sample particle tracking frame');
    if isequal(fname,0), disp('No file selected.'); params=[]; return; end
    imagePath = fullfile(fpath, fname);
end

info     = imfinfo(imagePath);
imgRaw   = imread(imagePath);
bitDepth = info(1).BitDepth;

if size(imgRaw,3) == 3
    % Greyscale cameras saved as 24-bit RGB have identical R=G=B channels.
    % Taking channel 1 directly avoids weighted-blend rounding artefacts.
    warning('Image loaded as RGB. Expected greyscale. Taking channel 1.');
    imgGray = imgRaw(:,:,1);
elseif size(imgRaw,3) == 1
    imgGray = imgRaw;
else
    error('Unexpected image dimensions: %s', mat2str(size(imgRaw)));
end

% Derive maxVal from the actual array class, not the raw file bit depth.
% This correctly handles 24-bit RGB (3x8-bit) files which yield uint8 after
% channel extraction, as well as true 16-bit greyscale images.
switch class(imgGray)
    case 'uint8',  maxVal = 255;
    case 'uint16', maxVal = 65535;
    otherwise,     maxVal = double(max(imgGray(:)));
end
[nRows, nCols] = size(imgGray);
imgDbl  = double(imgGray);
imgVec  = imgDbl(:);
imgVec  = imgVec(imgVec > 0);   % drop zero-intensity pixels

fprintf('\n=========================================\n');
fprintf(' STB Parameter Estimation\n');
fprintf('=========================================\n');
fprintf('  File       : %s\n', imagePath);
fprintf('  Resolution : %d x %d px,  %d-bit\n', nCols, nRows, bitDepth);
fprintf('  Mean / Std : %.1f / %.1f\n', mean(imgVec), std(imgVec));
fprintf('  p95 / Max  : %.0f / %.0f\n', prctile(imgVec,95), max(imgVec));

% ════════════════════════════════════════════════════════════════════════════
%  STEP 1 — HISTOGRAM  →  pick intensity threshold
% ════════════════════════════════════════════════════════════════════════════
nBins  = min(256, maxVal+1);
edges  = linspace(0, maxVal, nBins+1);
bCenters = 0.5*(edges(1:end-1)+edges(2:end));
counts = histcounts(imgVec, edges);
cumPct = 100 * cumsum(counts) / numel(imgVec);

figHist = figure('Name','Step 1 — Intensity Histogram  (non-zero pixels)');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
histogram(imgVec, nBins, 'BinLimits',[0 maxVal], ...
          'FaceColor',[0.2 0.5 0.8], 'EdgeColor','none');
xlabel('Pixel Intensity'); ylabel('Pixel Count');
title('Histogram'); grid on; xlim([0 maxVal]);
xline(prctile(imgVec,95),'g--',sprintf('p95=%.0f',prctile(imgVec,95)), ...
      'LineWidth',1.2,'LabelVerticalAlignment','bottom');
xline(prctile(imgVec,99),'r--',sprintf('p99=%.0f',prctile(imgVec,99)), ...
      'LineWidth',1.2,'LabelVerticalAlignment','bottom');

ax2 = nexttile;
plot(bCenters, cumPct, 'b-', 'LineWidth',1.8);
xlabel('Pixel Intensity'); ylabel('Cumulative (%)');
title('Cumulative Distribution'); grid on;
xlim([0 maxVal]); ylim([0 100]);
xline(prctile(imgVec,95),'g--','p95','LineWidth',1.2,'LabelVerticalAlignment','bottom');
xline(prctile(imgVec,99),'r--','p99','LineWidth',1.2,'LabelVerticalAlignment','bottom');

fprintf('\nSTEP 1 — Intensity Threshold\n');
fprintf('  Inspect the histogram. The threshold should sit in the valley\n');
fprintf('  between the background peak and the particle peak.\n');
fprintf('  Suggested starting points:  p95 = %.0f,  p99 = %.0f\n', ...
        prctile(imgVec,95), prctile(imgVec,99));

threshold = inputdlg( ...
    {sprintf('Enter intensity threshold  (range 0 – %d):', maxVal)}, ...
    'Intensity Threshold', 1, {sprintf('%.0f', prctile(imgVec,97))});
if isempty(threshold), disp('Cancelled.'); params=[]; return; end
threshold = str2double(threshold{1});

close(figHist);

% ════════════════════════════════════════════════════════════════════════════
%  STEP 2 — OVERLAY DETECTED PARTICLES  →  verify / adjust threshold
% ════════════════════════════════════════════════════════════════════════════
fprintf('\nSTEP 2 — Particle Detection Preview\n');

% Simple connected-component detection above threshold
bw      = imgDbl >= threshold;
bw      = bwareaopen(bw, 3);          % remove isolated 1-2 px specks
cc      = bwconncomp(bw);
props   = regionprops(cc, imgDbl, 'Centroid','EquivDiameter','MaxIntensity');
nFound  = cc.NumObjects;
fprintf('  Threshold = %.0f  →  %d particles detected\n', threshold, nFound);

figPrev = figure('Name','Step 2 — Detection Preview  (close when satisfied)');
axPrev  = axes('Parent', figPrev);
if maxVal > 255
    imshow(mat2gray(imgRaw), 'Parent', axPrev);
else
    imshow(imgRaw, 'Parent', axPrev);
end
axis(axPrev,'on'); hold(axPrev,'on');
set(axPrev,'YDir','reverse');
xlabel(axPrev,'x (pixels)'); ylabel(axPrev,'y (pixels)');
title(axPrev, sprintf('Threshold = %.0f  |  %d particles detected  |  Close window to continue', ...
      threshold, nFound));

if ~isempty(props)
    cents = vertcat(props.Centroid);
    radii = [props.EquivDiameter]/2;
    for k = 1:numel(props)
        viscircles(axPrev, cents(k,:), radii(k), ...
                   'Color','r', 'LineWidth',0.8, 'EnhanceVisibility',false);
    end
end

% Let user optionally re-enter threshold
uiwait(msgbox( ...
    sprintf('%d particles detected at threshold %.0f.\nClose this box — if the detection looks wrong, you can update the threshold next.', ...
    nFound, threshold), 'Step 2 — Check Detection','modal'));

answer = inputdlg( ...
    {'Update threshold? (leave unchanged to keep current value):'}, ...
    'Refine Threshold', 1, {sprintf('%.0f', threshold)});
if ~isempty(answer)
    newT = str2double(answer{1});
    if ~isnan(newT) && newT ~= threshold
        threshold = newT;
        bw    = imgDbl >= threshold;
        bw    = bwareaopen(bw, 3);
        cc    = bwconncomp(bw);
        props = regionprops(cc, imgDbl, 'Centroid','EquivDiameter','MaxIntensity');
        nFound = cc.NumObjects;
        fprintf('  Updated threshold = %.0f  →  %d particles detected\n', threshold, nFound);

        cla(axPrev); 
        if maxVal > 255
            imshow(mat2gray(imgRaw), 'Parent', axPrev);
        else
            imshow(imgRaw, 'Parent', axPrev);
        end
        hold(axPrev,'on');
        set(axPrev,'YDir','reverse');
        if ~isempty(props)
            cents = vertcat(props.Centroid);
            radii = [props.EquivDiameter]/2;
            for k = 1:numel(props)
                viscircles(axPrev, cents(k,:), radii(k), ...
                           'Color','r','LineWidth',0.8,'EnhanceVisibility',false);
            end
        end
        title(axPrev, sprintf('Threshold = %.0f  |  %d particles detected', threshold, nFound));
    end
end

% Auto-estimate radius from connected components
autoRadii = [props.EquivDiameter]/2;
autoRadius = median(autoRadii);
fprintf('  Auto-estimated tracer radius (median equiv. radius): %.2f px\n', autoRadius);

close(figPrev);

% ════════════════════════════════════════════════════════════════════════════
%  STEP 3 — MEASURE PARTICLE DIAMETER WITH RULER
% ════════════════════════════════════════════════════════════════════════════
fprintf('\nSTEP 3 — Measure Tracer Diameter\n');
fprintf('  Click on ONE edge of a well-focused particle, then the other edge.\n');
fprintf('  Press ''c'' to redo.  Close the window when done.\n');

diamMeasurements = [];

figRuler = figure('Name','Step 3 — Measure Particle Diameter  |  Click two edges of one particle');
axR = axes('Parent', figRuler);
if maxVal > 255
    imshow(mat2gray(imgRaw), 'Parent', axR);
else
    imshow(imgRaw, 'Parent', axR);
end
axis(axR,'on'); hold(axR,'on');
set(axR,'YDir','reverse');
xlabel(axR,'x (pixels)'); ylabel(axR,'y (pixels)');
title(axR,'Click both edges of a particle to measure diameter  |  ''c'' = clear  |  Close when done');

infoD = uicontrol('Parent',figRuler,'Style','text', ...
    'Units','normalized','Position',[0.05 0.01 0.90 0.04], ...
    'String','Click first edge of a particle.', ...
    'FontSize',10,'HorizontalAlignment','center', ...
    'BackgroundColor',get(figRuler,'Color'));

rulerState.pt1     = [];
rulerState.waiting = true;
rulerState.handles = {};

set(figRuler,'WindowButtonDownFcn', @(s,e) rulerClick(s,e,axR,nCols,nRows,infoD,'diam'));
set(figRuler,'KeyPressFcn',         @rulerKey);
set(figRuler,'CloseRequestFcn',     @(s,e) closeFig(s));

uiwait(figRuler);

if ~isempty(diamMeasurements)
    measuredDiam   = mean(diamMeasurements);
    measuredRadius = measuredDiam / 2;
    fprintf('  Measured diameter: %.2f px  →  radius = %.2f px\n', measuredDiam, measuredRadius);
else
    measuredRadius = autoRadius;
    fprintf('  No ruler measurement made — using auto-estimated radius: %.2f px\n', measuredRadius);
end

% ════════════════════════════════════════════════════════════════════════════
%  STEP 4 — MEASURE NEAREST-NEIGHBOUR SPACINGS
% ════════════════════════════════════════════════════════════════════════════
fprintf('\nSTEP 4 — Measure Particle Spacing\n');
fprintf('  Click centre-to-centre between several neighbouring particle pairs.\n');
fprintf('  Aim for ~5 measurements. Press ''c'' to clear. Close when done.\n');

spacingMeasurements = [];

figSpace = figure('Name','Step 4 — Measure Particle Spacing  |  Click centre-to-centre between neighbours');
axS = axes('Parent', figSpace);
if maxVal > 255
    imshow(mat2gray(imgRaw), 'Parent', axS);
else
    imshow(imgRaw, 'Parent', axS);
end
axis(axS,'on'); hold(axS,'on');
set(axS,'YDir','reverse');
xlabel(axS,'x (pixels)'); ylabel(axS,'y (pixels)');
title(axS,'Click centre-to-centre between neighbouring particles  |  ''c'' = clear  |  Close when done');

% Overlay detected centroids as guides
if ~isempty(props)
    cents = vertcat(props.Centroid);
    plot(axS, cents(:,1), cents(:,2), 'r+', 'MarkerSize',8, 'LineWidth',1);
end

infoS = uicontrol('Parent',figSpace,'Style','text', ...
    'Units','normalized','Position',[0.05 0.01 0.90 0.04], ...
    'String','Click centre of first particle.', ...
    'FontSize',10,'HorizontalAlignment','center', ...
    'BackgroundColor',get(figSpace,'Color'));

rulerState.pt1     = [];
rulerState.waiting = true;
rulerState.handles = {};

set(figSpace,'WindowButtonDownFcn', @(s,e) rulerClick(s,e,axS,nCols,nRows,infoS,'spacing'));
set(figSpace,'KeyPressFcn',         @rulerKey);
set(figSpace,'CloseRequestFcn',     @(s,e) closeFig(s));

uiwait(figSpace);

% Also compute auto nearest-neighbour from detected centroids
if ~isempty(props) && numel(props) > 1
    cents   = vertcat(props.Centroid);
    D       = pdist2(cents, cents);
    D(D==0) = Inf;
    autoSpacing = mean(min(D,[],2));
    fprintf('  Auto nearest-neighbour spacing (from %d detected particles): %.2f px\n', ...
            nFound, autoSpacing);
else
    autoSpacing = 10 * measuredRadius;
end

if ~isempty(spacingMeasurements)
    measuredSpacing = mean(spacingMeasurements);
    fprintf('  Ruler-measured spacing (mean of %d measurements): %.2f px\n', ...
            numel(spacingMeasurements), measuredSpacing);
    particleSpacing = measuredSpacing;
else
    particleSpacing = autoSpacing;
    fprintf('  No ruler measurements — using auto spacing: %.2f px\n', particleSpacing);
end

% ════════════════════════════════════════════════════════════════════════════
%  STEP 5 — COMPUTE ALL PARAMETERS
% ════════════════════════════════════════════════════════════════════════════
params.intensityThreshold = threshold;
params.tracerRadius       = measuredRadius;
params.searchRadius       = round(2.5 * measuredRadius * 10) / 10;   % 2.5× radius
params.minFrameCount      = 3;                                         % standard default
params.particleSpacing    = round(particleSpacing * 10) / 10;
params.gridResolution     = round(4 * particleSpacing * 10) / 10;     % ~4× spacing
params.shakeAmplitude     = 0.2;                                       % px, standard start
params.ghostThreshold     = round(0.4 * threshold);                   % 40% of int. threshold

fprintf('\n=========================================\n');
fprintf(' Estimated STB Parameters\n');
fprintf('=========================================\n');
fprintf('  intensityThreshold  : %.0f\n',  params.intensityThreshold);
fprintf('  tracerRadius        : %.2f px\n', params.tracerRadius);
fprintf('  --- Initial Phase ---\n');
fprintf('  searchRadius        : %.1f px\n', params.searchRadius);
fprintf('  minFrameCount       : %d frames\n', params.minFrameCount);
fprintf('  --- Convergence Phase ---\n');
fprintf('  particleSpacing     : %.1f px\n', params.particleSpacing);
fprintf('  --- Predict Field ---\n');
fprintf('  gridResolution      : %.1f px\n', params.gridResolution);
fprintf('  --- Shake ---\n');
fprintf('  shakeAmplitude      : %.2f px\n', params.shakeAmplitude);
fprintf('  ghostThreshold      : %.0f\n',  params.ghostThreshold);
fprintf('=========================================\n\n');

% ════════════════════════════════════════════════════════════════════════════
%  NESTED CALLBACKS  (shared by Steps 3 and 4)
% ════════════════════════════════════════════════════════════════════════════
    function rulerClick(~, ~, ax, nC, nR, txtH, mode)
        if gca ~= ax, return; end
        cp = get(ax,'CurrentPoint');
        x = cp(1,1);  y = cp(1,2);
        if x<0.5||x>nC+0.5||y<0.5||y>nR+0.5, return; end

        if rulerState.waiting
            rulerState.pt1 = [x,y];
            h = plot(ax, x, y, '+r', 'MarkerSize',14, 'LineWidth',2);
            rulerState.handles{end+1} = h;
            set(txtH,'String',sprintf('Pt 1: (%.1f, %.1f) — now click the second point.',x,y));
            rulerState.waiting = false;
        else
            x1=rulerState.pt1(1); y1=rulerState.pt1(2);
            dist = sqrt((x-x1)^2+(y-y1)^2);

            hl  = plot(ax,[x1 x],[y1 y],'r-','LineWidth',1.8);
            hm1 = plot(ax,x1,y1,'ro','MarkerFaceColor','r','MarkerSize',5);
            hm2 = plot(ax,x, y, 'ro','MarkerFaceColor','r','MarkerSize',5);
            ht  = text(ax,(x1+x)/2,(y1+y)/2,sprintf(' %.1f px',dist), ...
                       'Color','red','FontSize',9,'FontWeight','bold');
            rulerState.handles(end+1:end+4) = {hl,hm1,hm2,ht};

            if strcmp(mode,'diam')
                diamMeasurements(end+1) = dist; 
                set(txtH,'String',sprintf('Diameter = %.2f px.  Close window or click to measure again.',dist));
            else
                spacingMeasurements(end+1) = dist;
                set(txtH,'String',sprintf('Spacing = %.2f px  (n=%d).  Click next pair or close window.', ...
                    dist, numel(spacingMeasurements)));
            end

            rulerState.pt1    = [];
            rulerState.waiting = true;
        end
    end

    function rulerKey(~, evt)
        if strcmpi(evt.Key,'c')
            for i = 1:numel(rulerState.handles)
                if ishandle(rulerState.handles{i}), delete(rulerState.handles{i}); end
            end
            rulerState.handles = {};
            rulerState.pt1     = [];
            rulerState.waiting  = true;
        end
    end

    function closeFig(src)
        delete(src);   % releases uiwait
    end

end