%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nathan Wei
% Aerodynamics, Wind, And Renewable Energy (AWARE) Lab
% Mechanical Engineering and Applied Mechanics
% University of Pennsylvania
%
% Script for formatting figures for use in publications and presentations.
% Adapted from formatfigs.m from the Göttingen active-grid project (2017) 
% and an earlier script by Dr. Mike Krane (Penn State).
%
% Dependencies: export_fig 
%  (https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig)
% Created: 29.05.2018
% Updated: 14.11.2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ---------------------------- USER INPUTS ---------------------------- %%

pick_fig    = 1; % override fig_name / fig_dir and get figure manually if 1
save_fig    = 1; % export as .eps and save revised .fig if 1
fig_dir     = 'C:/Users/christopherbianco/Desktop/School_Code/LPT Data Analysis'; % where .fig lives

% These fields are only used if pick_fig = 0
fig_name    = 'all correlations paper.fig'; % figure to modify (.fig file)
save_name   = 'Fig9.eps'; % name to save .eps version under (include extension)

fsize_axesL = 18; % font size for axes labels [15] <-- example values given in brackets
fsize_axesT = 14; % font size for axes tick labels [12]
fsize_title = 14; % font size for plot title [15]
fsize_leg   = 14; % font size for legend [14]
linesize    = 2;  % line width [1]
markersize  = 8; % marker size [6]
resizeFig   = 1; % 1: resize figure according to figsize below; 0: don't resize
figsize     = [500, 500]; % dimensions for figure (pixels, W x H) [500, 400]
% 1000 width for full-width figures
% [330, 400] for 3 side-by-side plots (0.32\textwidth in LaTeX)
% [500, 400] for 2 side-by-side plots (0.48\textwidth in LaTeX)

renderer = 'painters'; %Renderer type

%% ---------------------------- ORGANIZATION --------------------------- %%

if pick_fig
    [fig_name, fig_dir] = uigetfile(fullfile(fig_dir, '*.fig'));
    save_name = [fig_name(1:end-4), '.eps'];
    if fig_name == 0
        return;
    end
end

%% ---------------------------- FORMATTING ----------------------------- %%

og_dir = pwd();

f1 = open(fullfile(fig_dir, fig_name)); % find and open figure
% Comment out below to forgo processing using %{

% Set parameters
set(f1, 'WindowStyle', 'normal');
set(f1, 'Renderer', renderer);
set(f1,'color','w'); % set background color to white

set(findall(f1,'-property','fontsize'),'fontsize',fsize_axesT);
axes = findall(gcf,'Type','Axes','-not','Type','AxesToolbar');
for aa = 1 : length(axes)
    axesLX = get(axes(aa), 'XLabel');
    set(axesLX, 'fontsize', fsize_axesL);
    axesLY = get(axes(aa), 'YLabel');
    set(axesLY, 'fontsize', fsize_axesL);
    axesLZ = get(axes(aa), 'ZLabel');
    set(axesLZ, 'fontsize', fsize_axesL);
    ttl = get(axes(aa), 'Title');
    set(ttl, 'fontsize', fsize_title, 'fontweight', 'bold', 'interpreter', 'latex');
    % delete(ttl); % uncomment if you don't want plot titles
end
cBar = findall(gcf,'Type','ColorBar');
if ~isempty(cBar)
    for ii = 1 : length(cBar)
        cBar(ii).Label.FontSize = fsize_axesL;
        % cBar(ii).Location = 'eastoutside'; % uncomment to change colorbar location
    end
end

objs = findall(f1,'-property','linewidth');
for ii = 1 : length(objs)
    if objs(ii).LineWidth < linesize
        set(objs(ii),'LineWidth',linesize); % avoid shrinking thicker lines
    end
end
% set(findall(f1,'-property','linewidth'),'linewidth',linesize);

set(findall(f1,'-property','MarkerSize'),'MarkerSize',markersize);
set(findall(f1,'-property','interpreter'),'interpreter','latex');
leg = findobj(f1, 'Type', 'Legend');
set(leg, 'fontsize', fsize_leg);
%set(leg, 'location', 'best'); % uncomment to change legend location

if resizeFig
    f1_size = [100, 100, figsize(1), figsize(2)]; % set new dimensions of figure
    set(f1,'Units','Pixels');
    set(f1,'OuterPosition',f1_size);
    set(f1,'PaperSize',f1_size(3:4));
end
grid on;

% Special handling instructions
q1 = findobj(f1, 'Type', 'quiver');
if ~isempty(q1)
    set(q1,'linewidth',0.5); % keep vectors in vector field small
end
%} % close comment for processing

if save_fig
    cd(fig_dir);
    saveas(gcf, fig_name); % save revised figure under same name
    export_fig(save_name, '-r300', '-eps', ['-', renderer]);
    export_fig(save_name, '-jpg', '-m3', ['-', renderer]);
    cd(og_dir);
    close(gcf);
end

%% Code to run export_fig over all MATLAB figures in a directory
%{
files = dir('*.fig');
for ii = 1 : length(files)
save_name = [files(ii).name(1:end-4), '.eps'];
f1 = open(fullfile(files(ii).name));
export_fig(save_name, '-r300', '-eps', ['-', renderer]);
close(f1);
end
%}