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

qinf = 0.5*1000*0.35*0.35*0.1*0.3; % dynamic pressure
drag_at_zero = 0.1286/qinf; % save drag at zero AoA


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

    CD = allData.CD_Mean;
    
    drag_bias =  CD(1);


    figure(1)
    hold on
    errorbar(aoa,allData.CL_Mean,allData.CL_Std./sqrt(14), 'DisplayName', leg);

    figure(2)
    hold on
    errorbar(aoa,CD - drag_bias + drag_at_zero ,allData.CD_Std./sqrt(14), 'DisplayName', leg);

    figure(3)
    hold on
    errorbar(aoa,allData.CM_Mean,allData.CM_Std./sqrt(14),  'DisplayName', leg);

end

figure(1)
hold off
figure(2)
hold off
figure(3)
hold off
