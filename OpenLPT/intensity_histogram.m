% Load filepath 
[fname, fpath] = uigetfile( ...
    {'*.tif;*.tiff;*.png;*.bmp;*.jpg;*.jpeg;*.pgm','Image Files';
    '*.*','All Files'}, ...
    'Select a particle tracking frame');
if isequal(fname, 0), disp('No file selected.'); return; end
imagePath = fullfile(fpath, fname);

imgRaw   = imread(imagePath);
bitDepth = imfinfo(imagePath);
bitDepth = bitDepth(1).BitDepth;

if size(imgRaw, 3) == 3
    % Unexpected RGB load -- take channel 1 to avoid weighted-blend artefacts.
    warning('Image loaded as RGB. Expected greyscale. Taking channel 1.');
    imgGray = imgRaw(:,:,1);
elseif size(imgRaw, 3) == 1
    imgGray = imgRaw;
else
    error('Unexpected image dimensions: %s', mat2str(size(imgRaw)));
end

switch class(imgGray)
    case 'uint8',  maxVal = 255;
    case 'uint16', maxVal = 65535;
    otherwise,     maxVal = double(max(imgGray(:)));
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

