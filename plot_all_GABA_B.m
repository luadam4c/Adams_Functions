%% plot_all_GABA_B.m
% Last updated 2018-12-18

% Get all subdirectories in present working directory
[~, directoryPaths] = all_subdirs;

% Count the number of subdirectories to go through
nDirs = numel(directoryPaths);

% Loop through all subdirectories
for iPath = 1:nDirs
    % Get this directory
    directoryThis = directoryPaths{iPath};

    % Print message
    fprintf('Plotting all figures for .abf files in %s\n', directoryThis);

    % Plot all figures and time
    tic;
    plot_all_abfs('Directory', directoryThis, 'ExpMode', 'patch');
    toc;
    
    % Close all figures
    close all
end

%% Post-process special file #1
% Special file #1
specialPath1 = fullfile('2019_01_04/003/', ...
                '2019_01_04_003_GABAB-IPSCs_features_table_by_sweep.xlsx');

% Parameters
outFolder1 = fullfile('2018_01_04_003', 'GABAB-IPSCs');
varsToPlot = {'peakAmplitude', 'peakDelayMs'};
xLabel = 'sweepNumberExceptLast';

% Load the table
specialTable1 = readtable(specialPath1);

% Remove last entry
specialTable1 = specialTable1(1:end-1, :);

% Replot tuning curve
h = plot_table(specialTable1, 'VariableNames', varsToPlot, ...
                'OutFolder', outFolder1, 'XLabel', xLabel, ...
                'LineWidth', 1, 'MarkerEdgeColor', rgb('DarkOrchid'), ...
                'MarkerFaceColor', rgb('LightSkyBlue'));
