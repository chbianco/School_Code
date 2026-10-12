%% Preamble
close all;
clear variables;
clc;

set(groot, 'defaultTextInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter', 'Latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextFontSize', 20);
set(groot, 'defaultAxesFontSize', 24);
set(groot, 'defaultLineLineWidth', 2);

%% Load Data
data = readmatrix('Turbulence_HW3_data.txt');

%% Zero out uw and vw

data(:, 13) = zeros(1, 541);
data(:, 14) = zeros(1, 541);

%% Deviatoric stress
