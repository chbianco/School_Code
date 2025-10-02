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

%% Loading Data 
%REMEMBER TO CHANGE THIS IF ON WINDOWS cause apples hates me 
fileDir = '/Users/christopherbianco/Desktop/School_Code/Wind Physics/HW1'; %Mac

%fileDir = 'C:\Users\Christopher\Desktop\School_Code\Wind Physics\HW1'; %Windows

data = load(fullfile(fileDir, '08_28_2019_22_00_00_000.mat'));