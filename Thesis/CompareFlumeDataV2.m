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
% prompt = 'How many experiments do you want to graph? ';
% x = input(prompt);


x = 3;
numDatasets = 5;
avgChoice = 'y';

All_CL = cell(1, x);
All_filtCL = cell(1, x);
All_times = cell(1, x);
All_sd = cell(1, x);
leg_vec = cell(1, x);

leg_vec{1} = "High Flex";
leg_vec{2} = "Low Flex";
leg_vec{3} = "Rigid";

for i = 1:x
    % Ask if user wants to average datasets
    % prompt = 'Do you want to average some datasets for this experiment? (y/n): ';
    % avgChoice = input(prompt, 's');
    
    % Initialize variables
    dataset_CL = {};
    dataset_filtCL = {};
    t_common = [];
    
    %Asks number of datasets to average if the user wants to average
    % if lower(avgChoice) == 'y'
    %     numDatasets = input('How many datasets to average? ');
    % else
    %     numDatasets = 1;
    % end
    
    %Starts averaging datasets
    for j = 1:numDatasets

        %If it is the first dataset, asks the user to select the file.
        %Otherwise, increments the iteration by one and picks that as the
        %next file
        if j == 1
        fprintf('Select file for experiment %d, dataset %d\n', i, j);
        [dataFile, dataPath] = uigetfile;

        % Split the filename using underscores
        fileParts = split(dataFile, '_');

        % Extract the part after the 5th underscore and increment it
        baseName = strjoin(fileParts(1:5), '_'); % Keep the first 5 parts as base
        suffixName = strjoin(fileParts(7:end), '_'); % Keep the rest after incrementing

        % Extract and increment the number after the 5th underscore
        numPart = str2double(fileParts{6}); % Convert to number

        else
        newNumPart = numPart + (j-1);
        newFileName = sprintf('%s_%d_%s', baseName, newNumPart, suffixName);
        %Assign the new file name
        dataFile = newFileName;
        end
        
        % Load the selected data file
        allData = load(strcat(dataPath, dataFile));
        time = allData.out.time;
        Lift = allData.out.FT11_Lift_N;

        % Filter the lift data
        fc = 10; fs = 5000;
        [b, a] = butter(6, fc/(fs/2));
        filtLift = filtfilt(b, a, Lift);

        % Compute CL
        u=0.35;
        qinf = 0.5*1000*0.35*0.35*0.1*0.3;
        CL = Lift/qinf;
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
    % prompt = 'Legend entry? ';
    % leg_vec{i} = input(prompt, 's');
end


% %% Sync up the times
% for j = 1:x
%     All_times{j} = All_times{j} - min(All_times{j});
% end
% 
% % Align all datasets to the first one
% ref_time = All_times{1};
% ref_CL = All_filtCL{1};
% 
% for j = 2:x
%     t = All_times{j};
%     C = All_filtCL{j};
% 
%     % Resample the reference and current dataset to common time points
%     t_common = linspace(min([ref_time(:); t(:)]), max([ref_time(:); t(:)]), max(length(ref_time), length(t)));
%     ref_CL_resampled = interp1(ref_time, ref_CL, t_common);
%     C_resampled = interp1(t, C, t_common, 'linear', 'extrap');
% 
%     % Compute cross-correlation
%     [xc, lags] = xcorr(ref_CL_resampled, C_resampled);
% 
%     % Find the lag with maximum correlation
%     [~, idx] = max(xc);
%     time_lag = lags(idx) * mean(diff(t_common));
% 
%     % Shift the current dataset by the computed time lag
%     All_times{j} = t + time_lag;
% end

%% Getting maxima
posMax = zeros(1,x);
negMax = zeros(1,x);
posMax_er = zeros(1,x);
negMax_er = zeros(1,x);

for j = 1:x
    posMax(j) = max(All_filtCL{j});
    negMax(j) = min(All_filtCL{j});

    % Find index of posMax and negMax in the data
    idx_pos = find(All_filtCL{j} == posMax(j), 1, 'first');
    idx_neg = find(All_filtCL{j} == negMax(j), 1, 'first');

    % Get corresponding standard deviations
    posMax_er(j) = All_sd{j}(idx_pos)/sqrt(5);
    negMax_er(j) = All_sd{j}(idx_neg)/sqrt(5);


end
% disp(posMax)
% disp(negMax)
disp(['Positive error is '])
disp(posMax_er)
disp(['Negative error is '])
disp(negMax_er)
%% Plotting
close all

color_vec = ["#0072BD", "#D95319", "#EDB120"];  % String array

%Initialize the figure
figure(1);
% Get the current position of the figure
pos = get(gcf, 'Position');  
% Scale the width and height by a factor (e.g., 1.5 times larger)
scaleFactor = 1.5;  % Adjust this factor as needed
newWidth = pos(3) * scaleFactor;
newHeight = pos(4) * scaleFactor;
% Keep the figure centered by adjusting the position
newX = pos(1) - (newWidth - pos(3)) / 2;
newY = pos(2) - (newHeight - pos(4)) / 2;
% Set the new position
set(gcf, 'Position', [newX, newY, newWidth, newHeight]);

hold on;
grid on;
xlabel('$t_c = \frac{tU}{L}$','Interpreter','Latex','FontSize',18);
ylabel('$C_L$','Interpreter','latex','FontSize',18);
xlim([10,40])

plot_handles = gobjects(1, x); % Store line handles for legend

% Plot each dataset and its corresponding legend
for j = 1:x
    t = All_times{j}(:).';  % Convert to row vector
    C = All_filtCL{j}(:).'; % Convert to row vector
    sd = All_sd{j}(:).';    % Convert to row vector
    se = sd./sqrt(5); %Standard error for five trials 
     
    C = C - mean(C(1:5000));

    C = -C;
    % Plot main line and store handle for legend
    plot_handles(j) = plot(t, C,'Color', color_vec(j) , 'LineWidth', 2);

     % Convert color to char for fill()
    fill_color = hex2rgb(color_vec(j));
    
    % Shaded region (fill) for standard deviation
    fill([t, fliplr(t)], [C + se, fliplr(C - se)], fill_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Create legend using only line handles
legend(plot_handles, leg_vec(~cellfun('isempty', leg_vec)));

%Make a title
ifTitle = input('Do you want a title? (y/n): ', 's');
titlename = ''; % Initialize in case user says no
if lower(ifTitle) == 'y'
    titlename = input('Enter title: ', 's');
end
title(titlename,'Interpreter','latex','FontSize',18);

hold off;

saveGraph = input('Do you want to save the graph? (y/n): ', 's');

if lower(saveGraph) == 'y'
    filename = input('Enter filename (without extension): ', 's');
    saveas(gcf, [filename, '.png']); % Save the figure as a PNG file
    saveas(gcf, [filename, '.fig']); % Save the figure as a FIG file
    fprintf('Figure saved as %s.png\n', filename);
end
