% Open the figure without displaying it (optional)
hFig = openfig('caseA_pos.fig','invisible');

% Find all axes in the figure
ax = findall(hFig,'Type','axes');

% Axes order is not guaranteed; sort by horizontal position
pos = arrayfun(@(a) a.Position(1), ax);
[~, idx] = sort(pos);

% Leftmost axes
axLeft = ax(idx(1));

% Create a new figure
newFig = figure;

% Copy the axes into the new figure
axNew = copyobj(axLeft, newFig);

% Make it fill the new figure nicely
set(axNew,'Position',[0.13 0.11 0.775 0.815]);

% Optional: make it current for editing
axes(axNew)