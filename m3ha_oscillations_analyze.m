%% Analyzes all GAT blocker oscillations data
%
% Requires:
%       cd/archive_dependent_scripts.m
%       cd/parse_all_multiunit.m
%       cd/plot_measures.m

% File History:
% 2019-10-17 Adapted from Glucose_analyze.m
% 2019-10-18 Changed minSpikeRateInBurstHz from 100 Hz to 200 Hz
% 2019-10-18 Changed minSpikeRateInBurstHz from 200 Hz to 100 Hz
% 2019-10-18 Changed maxInterBurstIntervalMs from 1500 ms to 1000 ms
% 2019-11-25 Added plotFigure1Population

%% Hard-coded parameters
figure01Dir = fullfile('/media', 'adamX', 'm3ha', 'manuscript', 'figures', 'Figure01');
parentDir = fullfile('/media', 'adamX', 'm3ha', 'oscillations');
% parentDir = fullfile('/media', 'shareX', 'Data_for_test_analysis', 'parse_multiunit_m3ha');
archiveDir = parentDir;
dirsToAnalyze = {'no711-final', 'snap5114-final', 'dual-final'};
% dirsToAnalyze = {'snap5114-final', 'dual-final'};
% dirsToAnalyze = {'no711-final'};
% dirsToAnalyze = {'dual-final'};
% dirsToAnalyze = {'no711-test', 'snap5114-test'};
% dirsToAnalyze = {'snap5114-test'};
specificSlicesToAnalyze = {};

plotFigure1Individual = false;
plotFigure1Population = false; % true;

parseIndividualFlag = true;
saveMatFlag = false; % true;
plotRawFlag = false; % true;
plotSpikeDetectionFlag = false; % true;
plotRasterFlag = false; % true;
plotSpikeDensityFlag = true;
plotSpikeHistogramFlag = true;
plotAutoCorrFlag = false; % true;
plotMeasuresFlag = true;
plotContourFlag = false; % true;
plotCombinedFlag = true;

parsePopulationRestrictedFlag = true;
plotChevronFlag = true;
plotByFileFlag = true;
plotByPhaseFlag = true;
plotNormByFileFlag = true;
plotNormByPhaseFlag = true;
plotPopAverageFlag = true;
plotSmoothNormPopAvgFlag = true;
parsePopulationAllFlag = false; %true;
plotAllMeasurePlotsFlag = false; %true;

archiveScriptsFlag = true; 

% For compute_default_signal2noise.m
relSnrThres2Max = 0.1;

% For detect_spikes_multiunit.m
filtFreq = [100, 1000];
minDelayMs = 25;

% For compute_spike_density.m
binWidthMs = 10;                % use a bin width of 10 ms by default
resolutionMs = 5;

% For compute_spike_histogram.m
% minBurstLengthMs = 20;          % bursts must be at least 20 ms by default
minBurstLengthMs = 100;          % bursts must be at least 100 ms by default
maxFirstInterBurstIntervalMs = 2000;
maxInterBurstIntervalMs = 1000; % bursts are no more than 
                                %   1 second apart
% maxInterBurstIntervalMs = 1500; % bursts are no more than 
%                                %   1.5 seconds apart
% minSpikeRateInBurstHz = 100;    % bursts must have a spike rate of 
%                                   at least 100 Hz by default
minSpikeRateInBurstHz = 30;    % bursts must have a spike rate of 
                                %   at least 30 Hz by default

% For compute_autocorrelogram.m
filterWidthMs = 100;
minRelProm = 0.02;

% For compute_phase_average.m & plot_measures.m
sweepsRelToPhase2 = -19:40;         % select between -20 & 40 min
% nSweepsLastOfPhase = 5;           % select from last 10 values of each phase
nSweepsLastOfPhase = 10;            % select from last 10 values of each phase
nSweepsToAverage = 5;               % select 5 values to average
% nSweepsToAverage = 10;            % select 10 values to average
selectionMethod = 'maxRange2Mean';  % average values within 40% of mean 
% selectionMethod = 'notNaN';       % average all values that are not NaNs
maxRange2Mean = 40;                 % range is not more than 40% of mean 
                                    %   by default
removeOutliersInPlot = false;

% For plot_measures.m
plotType = 'tuning';
sweepLengthSec = 60;
timeLabel = 'Time';
phaseLabel = 'Phase';
phaseStrings = {'Baseline', 'Wash-on', 'Wash-out'};
varsToPlot = {'oscIndex4'; 'oscPeriod2Ms'; ...
                    'oscDurationSec'; ...
                    'nSpikesTotal'; 'nSpikesInOsc'; ...
                    'nBurstsTotal'; 'nBurstsInOsc'; ...
                    'nSpikesPerBurst'; 'nSpikesPerBurstInOsc'};
varLabels = {'Oscillatory Index 4'; 'Oscillation Period 2 (ms)'; ...
                'Oscillation Duration (s)'; ...
                'Total Spike Count'; 'Number of Spikes in Oscillation'; ...
                'Total Number of Bursts'; 'Number of Bursts in Oscillation'; ...
                'Number of Spikes Per Burst'; ...
                'Number of Spikes Per Burst in Oscillation'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Chevron plots for Figure 01
if plotFigure1Population
    % Get all paths to Chevron tables
    [~, allSheetPaths] = all_files('Directory', figure01Dir, ...
                                'Suffix', 'chevron', 'Extension', 'csv');
                            
	% Read in all Chevron tables
    allChevronTables = cellfun(@(x) readtable(x, 'ReadRowNames', true), ...
                                allSheetPaths, 'UniformOutput', false);

    % Create figure names
    figPathBasesChevron = extract_fileparts(allSheetPaths, 'pathbase');
                            
	% Plot and save all Chevron tables
    cellfun(@(x, y) plot_and_save_chevron(x, y), ...
            allChevronTables, figPathBasesChevron);
end

% Run through all directories
for iDir = 1:numel(dirsToAnalyze)
    % Get the current directory to analyze
    dirThis = fullfile(parentDir, dirsToAnalyze{iDir});

    % Parse all slices in this directory
    if parseIndividualFlag
        parse_all_multiunit('Directory', dirThis, ...
                'SliceBases', specificSlicesToAnalyze, ...
                'SaveMatFlag', saveMatFlag, ...
                'PlotRawFlag', plotRawFlag, ...
                'PlotSpikeDetectionFlag', plotSpikeDetectionFlag, ...
                'PlotRasterFlag', plotRasterFlag, ...
                'PlotSpikeDensityFlag', plotSpikeDensityFlag, ...
                'PlotSpikeHistogramFlag', plotSpikeHistogramFlag, ...
                'PlotAutoCorrFlag', plotAutoCorrFlag, ...
                'PlotMeasuresFlag', plotMeasuresFlag, ...
                'PlotContourFlag', plotContourFlag, ...
                'PlotCombined', plotCombinedFlag, ...
                'RelSnrThres2Max', relSnrThres2Max, ...
                'FiltFreq', filtFreq, ...
                'MinDelayMs', minDelayMs, ...
                'BinWidthMs', binWidthMs, ...
                'ResolutionMs', resolutionMs, ...
                'MinBurstLengthMs', minBurstLengthMs, ...
                'MaxFirstInterBurstIntervalMs', maxFirstInterBurstIntervalMs, ...
                'MaxInterBurstIntervalMs', maxInterBurstIntervalMs, ...
                'MinSpikeRateInBurstHz', minSpikeRateInBurstHz, ...
                'FilterWidthMs', filterWidthMs, ...
                'MinRelProm', minRelProm, ...
                'NSweepsLastOfPhase', nSweepsLastOfPhase, ...
                'NSweepsToAverage', nSweepsToAverage, ...
                'MaxRange2Mean', maxRange2Mean);
    end
    
    if parsePopulationAllFlag
        % Plot measures for all phases
        plot_measures('Directory', dirThis, ...
                        'PlotAll', plotAllMeasurePlotsFlag, ...
                        'PlotChevronFlag', plotChevronFlag, ...
                        'PlotByFileFlag', plotByFileFlag, ...
                        'PlotByPhaseFlag', plotByPhaseFlag, ...
                        'PlotNormByFileFlag', plotNormByFileFlag, ...
                        'PlotNormByPhaseFlag', plotNormByPhaseFlag, ...
                        'PlotPopAverageFlag', plotPopAverageFlag, ...
                        'PlotSmoothNormPopAvgFlag', plotSmoothNormPopAvgFlag, ...
                        'RemoveOutliersInPlot', removeOutliersInPlot, ...
                        'NSweepsLastOfPhase', nSweepsLastOfPhase, ...
                        'NSweepsToAverage', nSweepsToAverage, ...
                        'SelectionMethod', 'notNaN', ...
                        'MaxRange2Mean', maxRange2Mean, ...
                        'PlotType', plotType, ...
                        'SweepLengthSec', sweepLengthSec, ...
                        'TimeLabel', timeLabel, ...
                        'PhaseLabel', phaseLabel, ...
                        'PhaseStrings', phaseStrings, ...
                        'VarsToPlot', varsToPlot, ...
                        'VarLabels', varLabels);
    end

    if parsePopulationRestrictedFlag
        % Plot measures for sweeps 1-lastSweepToMeasure only
        plot_measures('Directory', dirThis, ...
                        'SweepsRelToPhase2', sweepsRelToPhase2, ...
                        'PlotAll', plotAllMeasurePlotsFlag, ...
                        'PlotChevronFlag', plotChevronFlag, ...
                        'PlotByFileFlag', plotByFileFlag, ...
                        'PlotByPhaseFlag', plotByPhaseFlag, ...
                        'PlotNormByFileFlag', plotNormByFileFlag, ...
                        'PlotNormByPhaseFlag', plotNormByPhaseFlag, ...
                        'PlotPopAverageFlag', plotPopAverageFlag, ...
                        'PlotSmoothNormPopAvgFlag', plotSmoothNormPopAvgFlag, ...
                        'RemoveOutliersInPlot', removeOutliersInPlot, ...
                        'NSweepsLastOfPhase', nSweepsLastOfPhase, ...
                        'NSweepsToAverage', nSweepsToAverage, ...
                        'SelectionMethod', 'notNaN', ...
                        'MaxRange2Mean', maxRange2Mean, ...
                        'PlotType', plotType, ...
                        'SweepLengthSec', sweepLengthSec, ...
                        'TimeLabel', timeLabel, ...
                        'PhaseLabel', phaseLabel, ...
                        'PhaseStrings', phaseStrings, ...
                        'VarsToPlot', varsToPlot, ...
                        'VarLabels', varLabels);
    end
    
    close all
end

% Archive all scripts for this run
if archiveScriptsFlag
    archive_dependent_scripts(mfilename, 'OutFolder', archiveDir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_and_save_chevron(chevronTable, figPathBase)

% Extract figure base
figBase = extract_fileparts(figPathBase, 'dirbase');

% Extract drug name
drugName = extractBefore(figBase, '-');

% Make drug name all caps
drugNameAllCaps = upper(drugName);

% Create parameter tick labels
pTickLabels = {'Baseline'; drugNameAllCaps};

% Extract readout measure
measureName = extractAfter(extractBefore(figBase, '_chevron'), 'Phase2_');

% Decide on measure-dependent stuff
switch measureName
    case 'oscDurationSec'
        readoutLimits = [0, 25];
        readoutLabel = 'Oscillation Duration (s)';
    case 'oscPeriod2Ms'
        readoutLimits = [400, 1000];
        readoutLabel = 'Oscillation Period (ms)';
    otherwise
        error('measureName unrecognized!');
end

% Create figure
fig = set_figure_properties('AlwaysNew', true);

% Plot Chevron
plot_chevron(chevronTable, 'PlotMeanDifference', true, 'PlotErrorBars', true, ...
                'ColorMap', 'k', 'ReadoutLimits', readoutLimits, ...
                'PTickLabels', pTickLabels, ...
                'ReadoutLabel', readoutLabel, 'FigTitle', 'suppress', ...
                'LegendLocation', 'suppress');

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'inches', 'Width', 1.3, 'Height', 1, ...
                        'PlotMarkerSize', 2);
            
% Save figure
save_all_figtypes(fig, figPathBase, {'png', 'epsc'});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

parentDir = fullfile('/media', 'adamX', 'Glucose', 'oscillations', 'metformin');
lastSweepToMeasure = 45;        % select between sweeps 1:45

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%