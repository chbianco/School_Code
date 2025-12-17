% Open the .fig file
h = openfig('static_force_vs_aoa_curves_CD.fig');

% Find all error bar objects
eb = findobj(h, 'Type', 'errorbar');

% Define a list of marker symbols to cycle through
markers = {'o','s','*','^','v','>','<','p','h','x','+'};
legends = {'High Flex', 'Low Flex', 'Rigid'};

% Loop through error bar objects
for k = 1:length(eb)

    % Assign a unique marker (cycle if more objects than markers)
    eb(k).Marker = markers{mod(k-1, numel(markers)) + 1};

    % Increase sizes
    eb(k).MarkerSize = 10;
    eb(k).LineWidth  = 2;
    eb(k).CapSize    = 12;

    % Ensure marker face is visible (optional but often helpful)
    eb(k).MarkerFaceColor = eb(k).Color;

    % Set legend label
    eb(k).DisplayName = legends{k};
end

% Turn grid on
grid on

% Add legend
legend('Location','best')
