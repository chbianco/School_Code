%% PREAMBLE
close all;
clear variables;
clc;

set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultFigureColor', 'white');

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

% Initialize the figure
figure(1)
xlabel('Angle of Attack');
ylabel('$C_L$')
legend('show')

figure(2)
xlabel('Angle of Attack');
ylabel('$C_D$')
legend('show')

figure(3)
xlabel('Angle of Attack');
ylabel('$C_M$')
legend('show')

for i = 1:x
    
    % Prompts user to select file
    fprintf('Select file for experiment %d', i);
    [dataFile, dataPath] = uigetfile;

    %Prompts user for a legend name
    prompt = 'Legend entry? ';
    leg = input(prompt, 's');

    % Loads data file
    allData = load(strcat(dataPath, dataFile));
    aoa = allData.aoa_vec;

    figure(1)
    hold on
    errorbar(aoa,allData.CL_Mean,allData.CL_Std, 'DisplayName', leg);

    figure(2)
    hold on
    errorbar(aoa,allData.CD_Mean,allData.CD_Std, 'DisplayName', leg);

    figure(3)
    hold on
    errorbar(aoa,allData.CM_Mean,allData.CM_Std,  'DisplayName', leg);

end

figure(1)
hold off
figure(2)
hold off
figure(3)
hold off
