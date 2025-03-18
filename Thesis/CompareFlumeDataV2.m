%% Preamble
clear
close all

%% Set Working Directory
% Prompt the user to select a folder
folderPath = uigetdir;

% Check if the user selected a folder or canceled the selection
if folderPath ~= 0
    % Set the selected folder as the working directory
    cd(folderPath);
    fprintf('Working directory set to: %s\n', folderPath);
else
    fprintf('No folder selected. Working directory remains unchanged.\n');
end

%% Load and process data
prompt = 'How many experiments do you want to graph? ';
x = input(prompt);
All_CL = cell(1, x);
All_filtCL = cell(1, x);
All_times = cell(1, x);
All_sd = cell(1, x);
leg_vec = cell(1, x);

for i = 1:x
    % Ask if user wants to average datasets
    prompt = 'Do you want to average some datasets for this experiment? (y/n): ';
    avgChoice = input(prompt, 's');
    
    % Initialize variables
    dataset_CL = {};
    dataset_filtCL = {};
    t_common = [];

    if lower(avgChoice) == 'y'
        numDatasets = input('How many datasets to average? ');
    else
        numDatasets = 1;
    end
    
    for j = 1:numDatasets
        fprintf('Select file for experiment %d, dataset %d\n', i, j);
        [dataFile, dataPath] = uigetfile;
        
        % Load the selected data file
        allData = load(strcat(dataPath, dataFile));
        time = allData.out.time;
        Lift = allData.out.FT11_Lift_N;

        % Filter the lift data
        fc = 10; fs = 5000;
        [b, a] = butter(6, fc/(fs/2));
        filtLift = filtfilt(b, a, Lift);

        % Compute CL
        u = mean(abs(allData.out.Vectrino_x_ms));
        rho = 1000; L = 0.0751;
        CL = 2.*Lift/(rho*(u^2)*L);
        filtCL = filtfilt(b, a, CL); 

        % Convert time
        l = 0.1;
        t_conv = time .* u ./ l;

        % Ensure column vectors (important for correct averaging)
        t_conv = t_conv(:);
        CL = CL(:);
        filtCL = filtCL(:);

        % Define common time grid on first dataset
        if j == 1 || isempty(t_common)
            t_common = linspace(min(t_conv), max(t_conv), length(t_conv)).';
        end

        % Initialize resampled data
        CL_resampled = [];
        filtCL_resampled = [];

        % Resample only if enough data points exist
        if numel(t_conv) > 1 && numel(t_common) > 1
            CL_resampled = interp1(t_conv, CL, t_common, 'linear', 'extrap');
            filtCL_resampled = interp1(t_conv, filtCL, t_common, 'linear', 'extrap');
        else
            CL_resampled = CL;
            filtCL_resampled = filtCL;
        end

        % Store data in cell arrays
        dataset_CL{j} = CL_resampled;
        dataset_filtCL{j} = filtCL_resampled;
    end
    
    % Compute averages if needed
    if lower(avgChoice) == 'y'
        % Convert cell arrays to matrices where each column is a dataset
        valid_idx = ~cellfun(@isempty, dataset_CL);
        dataset_CL_matrix = cell2mat(dataset_CL(valid_idx));

        valid_idx = ~cellfun(@isempty, dataset_filtCL);
        dataset_filtCL_matrix = cell2mat(dataset_filtCL(valid_idx));

        % Ensure averaging across datasets (columns)
        avgCL = mean(dataset_CL_matrix, 2);  % Use mean along dim=2
        avg_filtCL = mean(dataset_filtCL_matrix, 2);
        std_filtCL = std(dataset_filtCL_matrix, 0, 2); % Standard deviation along dim=2

        All_CL{i} = avgCL;
        All_filtCL{i} = avg_filtCL;
        All_times{i} = t_common;
        All_sd{i} = std_filtCL;
    else
        All_CL{i} = CL;
        All_filtCL{i} = filtCL;
        All_times{i} = t_conv;
        All_sd{i} = zeros(size(filtCL));
    end

    % Legend entry
    prompt = 'Legend entry? ';
    leg_vec{i} = input(prompt, 's');
end


%% Sync up the times
for j = 1:x
    All_times{j} = All_times{j} - min(All_times{j});
end

% Align all datasets to the first one
ref_time = All_times{1};
ref_CL = All_filtCL{1};

for j = 2:x
    t = All_times{j};
    C = All_filtCL{j};
    
    % Resample the reference and current dataset to common time points
    t_common = linspace(min([ref_time(:); t(:)]), max([ref_time(:); t(:)]), max(length(ref_time), length(t)));
    ref_CL_resampled = interp1(ref_time, ref_CL, t_common);
    C_resampled = interp1(t, C, t_common, 'linear', 'extrap');
    
    % Compute cross-correlation
    [xc, lags] = xcorr(ref_CL_resampled, C_resampled);
    
    % Find the lag with maximum correlation
    [~, idx] = max(xc);
    time_lag = lags(idx) * mean(diff(t_common));
    
    % Shift the current dataset by the computed time lag
    All_times{j} = t + time_lag;
end

%% Plotting
close all

color_vec = ['b', 'r', 'y'];

figure(1);
hold on;
grid on;
xlabel('$t_c = \frac{tU}{L}$','Interpreter','Latex','FontSize',12);
ylabel('$C_L$','Interpreter','latex','FontSize',12);
xlim([10,40])

plot_handles = gobjects(1, x); % Store line handles for legend

% Plot each dataset and its corresponding legend
for j = 1:x
    t = All_times{j}(:).';  % Convert to row vector
    C = All_filtCL{j}(:).'; % Convert to row vector
    sd = All_sd{j}(:).';    % Convert to row vector
    
    C = C - mean(C(1:5000));
    % Plot main line and store handle for legend
    plot_handles(j) = plot(t, C,'Color', color_vec(j) , 'LineWidth', 2);

    % Shaded region (fill) for standard deviation, but exclude from legend
    fill([t, fliplr(t)], [C + sd, fliplr(C - sd)], color_vec(j), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Create legend using only line handles
legend(plot_handles, leg_vec(~cellfun('isempty', leg_vec)));

%Make a title
ifTitle = input('Do you want a title? (y/n): ', 's');
titlename = ''; % Initialize in case user says no
if lower(ifTitle) == 'y'
    titlename = input('Enter title: ', 's');
end
title(titlename,'Interpreter','latex','FontSize',12);

hold off;

saveGraph = input('Do you want to save the graph? (y/n): ', 's');

if lower(saveGraph) == 'y'
    filename = input('Enter filename (without extension): ', 's');
    saveas(gcf, [filename, '.png']); % Save the figure as a PNG file
    fprintf('Figure saved as %s.png\n', filename);
end
