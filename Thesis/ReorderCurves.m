% Open the .fig file
h = openfig('caseC_lift.fig');

% Get axes handle
ax = gca;

% Find all line objects in this axes
L = findobj(ax, 'Type', 'line');

% findobj returns objects in reverse stacking order
% Flip so L(1) = first plotted (bottom)
L = flipud(L);

% Move the first-plotted line to the top
uistack(L(1), 'top');

% -------- LEGEND CONTROL --------
% User-defined legend labels (edit as needed)
legends = {
    'Low Flex'
    'Rigid'
    'High Flex'
};

% Re-fetch lines in correct (visual) order
L = flipud(findobj(ax, 'Type', 'line'));

% Assign legend names
for k = 1:length(L)
    if k <= numel(legends)
        L(k).DisplayName = legends{k};
    else
        L(k).DisplayName = sprintf('Curve %d', k);
    end
end

% Create / update legend
lgd = legend(ax, 'Location', 'best');

% -------- FONT SIZE CONTROL --------
axisLabelFontSize = 36;   % <-- change as desired
legendFontSize    = 28;   % <-- change as desired

% Axis labels
ax.XLabel.FontSize = axisLabelFontSize;
ax.YLabel.FontSize = axisLabelFontSize;
ax.FontSize = 24;   % controls x- and y-tick label font size


% Legend text
lgd.FontSize = legendFontSize;

% Optional grid
grid on
