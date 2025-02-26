% Script to generate brain plots for illustrating the concept of the paper.
% Random values were chosen for illustration purposes of the concept behind
% the actual findings.

% Code for Simple-Brain-Plot can be obtained from
% https://github.com/dutchconnectomelab/Simple-Brain-Plot/tree/main. 
% Change sw_dir accordingly. 



clear; clc;
sw_dir = '/code/Simple-Brain-Plot';  % adapt path to Simple-Brain-Plot
addpath(sw_dir);

wd = pwd();
outdir = fullfile(wd, 'outputs', 'concept_figure');
mkdir(outdir);

%% define brain plot variables
% Create a gray color map
gray_map = [1, 1, 1; 0, 0, 0];   % Start with white and end with black
numColors = 34;
x = linspace(1, size(gray_map, 1), numColors);
cm = interp1(1:size(gray_map, 1), gray_map, x, 'linear');
colormap(cm);
set(gcf, 'Color', 'w');
saveas(gcf, fullfile(outdir, 'colorbar_grayscale.svg'));


% also create a red to green colorbar for the environment (taken from Simple-Brain-Plot git)
cm_env = [0.8431    0.1882    0.1529; ...
        0.9569    0.4275    0.2627; ...
        0.9922    0.6824    0.3804; ...
        0.9961    0.8784    0.5451; ...
        1.0000    1.0000    0.7490; ...
        0.8510    0.9373    0.5451; ...
        0.6510    0.8510    0.4157; ...
        0.4000    0.7412    0.3882; ...
        0.1020    0.5961    0.3137];
%cm_env_smooth = interp1(1:size(cm_env, 1), cm_env, linspace(1, size(cm_env, 1), 500), 'linear');
cm_env = interp1(cm_env, 1:0.01:size(cm_env,1));
colormap(cm_env);
c = colorbar;
c.Ticks = [];
c.Box = 'off';
set(gcf, 'Color', 'w');
saveas(gcf, fullfile(outdir, 'colorbar_env.svg'));

% load aparc regions
load(fullfile(sw_dir,'examples', 'regionDescriptions.mat'));

%% Group difference
% age birth
values = zeros(size(regionDescriptions.aparc_aseg));
values(55) = 0.1; % inferiorparietal
values(59) = 0.05; % lateralorbitofrontal
values(62) = 0.15; % middletemporal
values(74) = 0.15; % rostral middle frontal
values(75) = 0.3; % superior frontal
values(77) = 0.1; % superior temporal
values(78) = 0.05; % supramarginal
values(80) = 0.05; % temporal pole
values(82) = 0.3; % insula

outname = fullfile(outdir, 'group_diff_age_30w.svg');
plotBrain(regionDescriptions.aparc_aseg, values, cm, 'atlas', 'aparc', ...
    'limits', [0,1],'savePath', outname);

%--------------------------------------------
% age 7
values = zeros(size(regionDescriptions.aparc_aseg));
values(55) = 0.25; % inferiorparietal
values(59) = 0.15; % lateralorbitofrontal
values(62) = 0.25; % middletemporal
values(74) = 0.1; % rostral middle frontal
values(75) = 0.25; % superior frontal
values(77) = 0.1; % superior temporal
values(78) = 0.15; % supramarginal
values(80) = 0.1; % temporal pole
values(82) = 0.25; % insula

outname = fullfile(outdir, 'group_diff_age_7y.svg');
plotBrain(regionDescriptions.aparc_aseg, values, cm, 'atlas', 'aparc', ...
    'limits', [0,1],'savePath', outname);

%---
% age 25
values = zeros(size(regionDescriptions.aparc_aseg));
values(55) = 0.35; % inferiorparietal
values(59) = 0.2; % lateralorbitofrontal
values(62) = 0.5; % middletemporal
values(74) = 0.15; % rostral middle frontal
values(75) = 0.35; % superior frontal
values(77) = 0.1; % superior temporal
values(78) = 0.15; % supramarginal
values(80) = 0.2; % temporal pole
values(82) = 0.15; % insula

outname = fullfile(outdir, 'group_diff_age_25y.svg');
plotBrain(regionDescriptions.aparc_aseg, values, cm, 'atlas', 'aparc', ...
    'limits', [0,1], 'savePath', outname);



%% Individual plots IBAPs
% Subject 1 (red SES), age birth
values = zeros(size(regionDescriptions.aparc_aseg));
values(62) = 0.6; % middletemporal
values(77) = 0.4; % superior temporal
values(69) = 0.5; % postcentral
values(75) = 0.5; % superior frontal
values(74) = 0.6; % rostral middle frontal
values(66) = 0.7; % pars orbitalis
values(58) = 0.3; % lateral occipital
values(80) = 0.3; % temporal pole

outname = fullfile(outdir, 'sub1_age_birth.svg');
plotBrain(regionDescriptions.aparc_aseg, values, cm, 'atlas', 'aparc', ...
    'limits', [0,1],'savePath', outname);

%---
% Subject 1, age 7y
values = zeros(size(regionDescriptions.aparc_aseg));
values(62) = 0.9; % middletemporal
values(77) = 0.6; % superior temporal
values(69) = 0.8; % postcentral
values(75) = 0.7; % superior frontal
values(74) = 0.6; % rostral middle frontal
values(66) = 0.8; % pars orbitalis
values(58) = 0.3; % lateral occipital
values(80) = 0.5; % temporal pole

outname = fullfile(outdir, 'sub1_age_7y.svg');
plotBrain(regionDescriptions.aparc_aseg, values, cm, 'atlas', 'aparc', ...
    'limits', [0,1],'savePath', outname);

%---
% Subject 1, age 25y
values = zeros(size(regionDescriptions.aparc_aseg));
values(62) = 1; % middletemporal
values(77) = 0.8; % superior temporal
values(69) = 1; % postcentral
values(75) = 0.85; % superior frontal
values(74) = 0.7; % rostral middle frontal
values(66) = 0.9; % pars orbitalis
values(58) = 0.5; % lateral occipital
values(80) = 0.8; % temporal pole

outname = fullfile(outdir, 'sub1_age_25y.svg');
plotBrain(regionDescriptions.aparc_aseg, values, cm, 'atlas', 'aparc', ...
    'limits', [0,1],'savePath', outname);


%% Subject 2
% Subject 2, birth
values = zeros(size(regionDescriptions.aparc_aseg));
values(78) = 0.7; % supramarginal
values(65) = 0.7; % pars opercularis
values(51) = 0.5; % caudal middle frontal
values(56) = 0.9; % inferior temporal
values(59) = 0.6; % lateral orbitofrontal


outname = fullfile(outdir, 'sub2_age_birth.svg');
plotBrain(regionDescriptions.aparc_aseg, values, cm, 'atlas', 'aparc', ...
    'limits', [0,1],'savePath', outname);

%---
% Subject 2, 7 y
values = zeros(size(regionDescriptions.aparc_aseg));
values(78) = 0.3; % supramarginal
values(65) = 0.75; % pars opercularis
values(51) = 0.5; % caudal middle frontal
values(56) = 0.65; % inferior temporal
values(59) = 0.7; % lateral orbitofrontal


outname = fullfile(outdir, 'sub2_age_7y.svg');
plotBrain(regionDescriptions.aparc_aseg, values, cm, 'atlas', 'aparc', ...
    'limits', [0,1],'savePath', outname);

%---
% Subject 2, 25 y
values = zeros(size(regionDescriptions.aparc_aseg));
values(78) = 0.4; % supramarginal
values(65) = 0.75; % pars opercularis
values(51) = 0.45; % caudal middle frontal
values(56) = 0.75; % inferior temporal
values(59) = 0.7; % lateral orbitofrontal


outname = fullfile(outdir, 'sub2_age_25y.svg');
plotBrain(regionDescriptions.aparc_aseg, values, cm, 'atlas', 'aparc', ...
    'limits', [0,1],'savePath', outname);


%% Subject 3
% Subject 3, birth
values = zeros(size(regionDescriptions.aparc_aseg));
values(69) = 0.7; % postcentral
values(77) = 0.3; % superior temporal
values(81) = 0.5; % transverse temporal
values(62) = 0.2; % middletemporal


outname = fullfile(outdir, 'sub3_age_birth.svg');
plotBrain(regionDescriptions.aparc_aseg, values, cm, 'atlas', 'aparc', ...
    'limits', [0,1],'savePath', outname);

%---
% Subject 3, 7y
values = zeros(size(regionDescriptions.aparc_aseg));
values(69) = 0.6; % postcentral
values(77) = 0.25; % superior temporal
values(81) = 0.7; % transverse temporal


outname = fullfile(outdir, 'sub3_age_7y.svg');
plotBrain(regionDescriptions.aparc_aseg, values, cm, 'atlas', 'aparc', ...
    'limits', [0,1],'savePath', outname);

%---
% Subject 3, 25y
values = zeros(size(regionDescriptions.aparc_aseg));
values(69) = 0.6; % postcentral
values(77) = 0.5; % superior temporal
values(81) = 0.8; % transverse temporal

outname = fullfile(outdir, 'sub3_age_25y.svg');
plotBrain(regionDescriptions.aparc_aseg, values, cm, 'atlas', 'aparc', ...
    'limits', [0,1],'savePath', outname);


%----------
close all;