% plot_wand_points.m
% Plots wand detections from wand_points.csv, one subplot per camera.
% PointIdx 0 (Small) = green, PointIdx 1 (Large) = red

data = readtable('results/wand_points.csv');

nCams = 4;
figure('Position', [100 100 1400 900]);

for c = 0:nCams-1
    ax = subplot(2, 2, c+1);
    hold(ax, 'on');

    camMask = data.Camera == c;

    % PointIdx 0 — Small — green
    mask0 = camMask & data.PointIdx == 0;
    scatter(ax, data.X(mask0), data.Y(mask0), 6, 'g', 'filled', ...
        'DisplayName', 'Idx 0 (Small)');

    % PointIdx 1 — Large — red
    mask1 = camMask & data.PointIdx == 1;
    scatter(ax, data.X(mask1), data.Y(mask1), 6, 'r', 'filled', ...
        'DisplayName', 'Idx 1 (Large)');

    set(ax, 'YDir', 'reverse');   % image coords: Y increases downward
    axis(ax, 'tight');
    grid(ax, 'on');
    title(ax, sprintf('Camera %d', c));
    xlabel(ax, 'X (px)');
    ylabel(ax, 'Y (px)');
    legend(ax, 'Location', 'best');
    hold(ax, 'off');
end

sgtitle('Wand Detections by Camera');