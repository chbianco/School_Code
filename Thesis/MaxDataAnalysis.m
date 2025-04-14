%% PREAMBLE
close all;
clear variables;
clc;

% Set default formatting properties for figures using LaTeX formatting
set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultFigureColor', 'white');

%% Initializing and Loading Data
leg_vec = ["High Flex", "Low Flex", "Rigid"];
color_vec = ["#0072BD", "#D95319", "#EDB120"];  % String array


%Column 1 is scenario, 2 is encounter location, 3 is circulation, 4 ,is
%CL1, 5 is CL2, 6 is r
rigid = readmatrix('C:\Users\Christopher\Desktop\Good Thesis Data\MaxData\RigidMaxData.csv');
highFlex = readmatrix('C:\Users\Christopher\Desktop\Good Thesis Data\MaxData\HighFlexMaxData.csv');
lowFlex = readmatrix('C:\Users\Christopher\Desktop\Good Thesis Data\MaxData\LowFlexMaxData.csv');

%% Plotting
figure(1)
ar1 = subplot(1, 2, 1);
hold(ar1, 'on');
xlabel(ar1, 'Encounter Location $y^* = y/L$');
ylabel(ar1, '$C_{L, max}$');
title(ar1, 'Maximum Lift');
%xlim(ar1, [10 40]);
ylim(ar1, [-0.8 1.1]);
grid(ar1, 'on');
scatter(highFlex(:,2), highFlex(:,5), 50, 'filled', ...
    'MarkerFaceColor', color_vec(1), ...
    'MarkerEdgeColor', color_vec(1), ...
    'DisplayName', leg_vec(1));
scatter(lowFlex(:,2), lowFlex(:,5), 50, 'filled', ...
    'MarkerFaceColor', color_vec(2), ...
    'MarkerEdgeColor', color_vec(2), ...
    'DisplayName', leg_vec(2));
scatter(rigid(:,2), rigid(:,5), 50, 'filled', ...
    'MarkerFaceColor', color_vec(3), ...
    'MarkerEdgeColor', color_vec(3), ...
    'DisplayName', leg_vec(3));
legend('show')
hold(ar1, 'off')

ar2 = subplot(1, 2, 2);
hold(ar2, 'on');
xlabel(ar2, 'Encounter Location $y^* = y/L$');
ylabel(ar2, '$r^*_{max}$');
title(ar2, 'Maximum Laser Displacement');
%xlim(ar2, [10 40]);
%ylim(ar2, [-0.4 0.4]);
grid(ar2, 'on');
scatter(highFlex(:,2), highFlex(:,6), 50, 'filled', ...
    'MarkerFaceColor', color_vec(1), ...
    'MarkerEdgeColor', color_vec(1), ...
    'DisplayName', leg_vec(1));
scatter(lowFlex(:,2), lowFlex(:,6), 50, 'filled', ...
    'MarkerFaceColor', color_vec(2), ...
    'MarkerEdgeColor', color_vec(2), ...
    'DisplayName', leg_vec(2));
scatter(rigid(:,2), rigid(:,6), 50, 'filled', ...
    'MarkerFaceColor', color_vec(3), ...
    'MarkerEdgeColor', color_vec(3), ...
    'DisplayName', leg_vec(3));
legend('show')
hold(ar2, 'off')

figure(2)

ag1 = subplot(1, 2, 1);
hold(ag1, 'on');
xlabel(ag1, 'Circulation $\Gamma^* = \overline{\Gamma}/cU$');
ylabel(ag1, '$C_{L, max}$');
title(ag1, 'Maximum Lift');
%xlim(ag1, [10 40]);
ylim(ag1, [-0.8 1.1]);
grid(ag1, 'on');
scatter(highFlex(:,3), highFlex(:,5), 50, 'filled', ...
    'MarkerFaceColor', color_vec(1), ...
    'MarkerEdgeColor', color_vec(1), ...
    'DisplayName', leg_vec(1));
scatter(lowFlex(:,3), lowFlex(:,5), 50, 'filled', ...
    'MarkerFaceColor', color_vec(2), ...
    'MarkerEdgeColor', color_vec(2), ...
    'DisplayName', leg_vec(2));
scatter(rigid(:,3), rigid(:,5), 50, 'filled', ...
    'MarkerFaceColor', color_vec(3), ...
    'MarkerEdgeColor', color_vec(3), ...
    'DisplayName', leg_vec(3));
legend('show')
hold(ag1, 'off')

ag2 = subplot(1, 2, 2);
hold(ag2, 'on');
xlabel(ag2, 'Circulation $\Gamma^* = \overline{\Gamma}/cU$');
ylabel(ag2, '$r^*_{max}$');
title(ag2, 'Maximum Laser Displacement');
%xlim(ag2, [10 40]);
%ylim(ag2, [-0.4 0.4]);
grid(ag2, 'on');
scatter(highFlex(:,3), highFlex(:,6), 50, 'filled', ...
    'MarkerFaceColor', color_vec(1), ...
    'MarkerEdgeColor', color_vec(1), ...
    'DisplayName', leg_vec(1));
scatter(lowFlex(:,3), lowFlex(:,6), 50, 'filled', ...
    'MarkerFaceColor', color_vec(2), ...
    'MarkerEdgeColor', color_vec(2), ...
    'DisplayName', leg_vec(2));
scatter(rigid(:,3), rigid(:,6), 50, 'filled', ...
    'MarkerFaceColor', color_vec(3), ...
    'MarkerEdgeColor', color_vec(3), ...
    'DisplayName', leg_vec(3));
legend('show')
hold(ag2, 'off')
