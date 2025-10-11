%% Analyzes rat sniff whisk data 
% 
% This script must be run in a directory that contains the data directory
% specified below.
%
% Dataset: sniff_whisk_data
% Description from Jeff: 
%   - unilateral whisking and sniffing data, head restrained
%   - whisking data are from tracking a single vibrissa in full-field video
%   - breathing data are from a thermocouple
%   - some trials contain puffs of ammonia or control puffs of unscented air to the nose
%       time of application is indicated by the pulse in channel "piezo"
%   - check sign
%
% Requires:
%       cd/all_files.m
%       cd/array_fun.m
%       cd/check_dir.m
%       cd/compute_combined_trace.m
%       cd/count_vectors.m
%       cd/create_labels_from_numbers.m
%       cd/extract_fields.m
%       cd/extract_fileparts.m
%       cd/extract_subvectors.m
%       cd/find_matching_files.m
%       cd/parse_oscillation.m
%       cd/plot_traces.m
%       cd/plot_vertical_line.m
%       cd/plot_vertical_shade.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/virt_plot_amplitude_correlation.m
%       cd/virt_plot_log_decrement_jitter.m
%       cd/virt_plot_phase_response.m
%       cd/write_table.m
%       cd/virt_detect_whisk_analysis_windows.m
%
% Used by:

% File History:
% 2025-08-31 Created by Adam Lu
% 2025-09-05 Now detects sniff transitions
% 2025-09-05 Now detects and plots whisk peaks and valleys
% 2025-09-05 Now stores sniff and whisk fundamental frequencies in metadata
% 2025-09-06 Now finds and plots sniff start windows based on sniff/whisk criteria
% 2025-09-09 Added a third plot for individual sniff start windows
% 2025-09-09 Now calculates average whisk amplitude ratios
% 2025-09-10 Changed to compute whisk logarithmic decrements
% 2025-09-11 Now computes averages statistics for whisk logarithmic decrements
% 2025-09-11 Renamed analysis window as sniff start window
% 2025-09-12 Updated table outputs to be saved
% 2025-09-12 Fixed sniff vec filter cutoff at [1, 15] Hz, filter order changed to 3
% 2025-09-12 Fixed whisk vec filter cutoff at [3, 25] Hz, filter order changed to 3
% 2025-09-12 Added minPeakPromWhisk
% 2025-09-12 Added analysis of basal respiration cycles
% 2025-09-12 Changed promThresholdPerc from 10 to 5
% 2025-09-13 Added analysis of amplitude correlations
% 2025-09-13 Changed minPeakPromWhisk to 5 and Added maxWhiskDurationMs at 250 ms
%               and force analysis windows to have whisks without NaN amplitudes
% 2025-09-15 Attempted new first whisk in analysis window definition: first
%               prevalley (protraction) after 30 ms before the resp prevalley
% 2025-09-16 Changed definition of first whisk in analysis window back
% 2025-09-18 Now plots regression lines
% 2025-09-18 Now plots phase response detections and phase response curves
% 2025-09-19 Changed promThresholdPercResp from 5 to 15
% 2025-09-19 Changed phase detection to use whisk valleys (protraction) and breathOnsetLatencyMs
% 2025-10-06 Refactored to use virt_plot_* functions for aggregate plots
% 2025-10-06 Refactored to use virt_detect_whisk_analysis_windows.m by Gemini
% 2025-10-09 Modified by Gemini to update combine_tables function.
% 2025-10-09 Modified by Gemini to pass trialName to parse_whisk_vecs.

%% Hard-coded parameters
% Input Directory and file naming conventions
nameDataDir = 'data_sniff_whisk';       % Name of the directory containing the data files
extDataFile = 'mat';                    % Extension for the data files
prefixDataFile = 'swdata_';             % Prefix for the data files
suffixDataFile = '_angle';              % Suffix for the data files
suffixParamFile = '_parameters';        % Suffix for the parameter files

% Output Directory and file naming conventions
nameOutDir = 'output_sniff_whisk';      % Name of the output directory
fileNameMetaData = 'sniff_whisk_metadata.csv';               % File name of the metadata table file
fileNameSniffStartWinTable = 'sniff_start_window_table.csv'; % File name of the sniff start window table file
fileNameBasalRespCycleTable = 'basal_resp_cycle_table.csv';  % File name of the basal respiration cycle table file
fileNameAnalysisResults = 'sniff_whisk_analysis.mat';        % File name of the analysis results mat file
fileNameAlgorithmDiffs = 'sniff_whisk_algorithm_differences.txt'; % File name for algorithm differences log

% Detection parameters
ammoniaPuffString = 'ammpuff';
airPuffString = 'airpuff';
minPulseAmplitude = 2;                  % Minimum pulse amplitude for piezo trace

% Analysis parameters
%   Note: Keep this consistent with virt_moore.m
invertThermSignal = true;
invertWhiskSignal = true;
amplitudeDefinition = 'peak-to-avgvalley';
whiskDirForPhase = 'protraction';   % whisk direction for phase calculations
%whiskDirForPhase = 'retraction';   % whisk direction for phase calculations
fundFreqRange = [0.5, 20];          % range of possible fundamental frequencies to be detected
fCutoffResp = [1, 15];              % bandpass filter cutoff for resp trace (Moore et al 2013 used [1, 15] Hz)
fCutoffWhisk = [3, 25];             % bandpass filter cutoff for whisk trace (Moore et al 2013 used [3, 25] Hz)
fCutoffRelToFund = [];              % don't use this
filterOrderResp = 3;                % Butterworth filter order for resp trace (Moore et al 2013 used 3)
filterOrderWhisk = 3;               % Butterworth filter order for whisk trace (Moore et al 2013 used 3)
promThresholdPercResp = 15;         % Percentage of amplitude range for minimum peak prominence for resp peaks
promThresholdPercWhisk = 5;         % Percentage of amplitude range for minimum peak prominence for whisk peaks
minPeakPromWhisk = 5;               % Minimum whisk angle change (degrees) to detect as a peak (Moore et al 2013 used 5 degrees)
maxWhiskDurationMs = 250;           % Maximum whisk inter-valley interval (ms) (Moore et al 2013 used whisk duration < 250 ms)
minPeakDistanceMsResp = 30;         % Minimum peak distance (ms) for resp peaks
minPeakDistanceMsWhisk = 30;        % Minimum peak distance (ms) for whisk peaks
breathOnsetLatencyMs = 60; %30;     % Presumed latency (ms) for from PB neuron activation to breath onset (use for phase analysis, Moore et al 2013 used 30 ms)
sniffFreqThreshold = 4;             % Frequency threshold for sniffing in Hz
basalFreqThreshold = 3;             % Frequency threshold for basal respiration in Hz
nWhisksSniffStartToAnalyze = 5;     % Number of whisks at the start of a sniff period to be analyzed
minWhisksBasalRespToAnalyze = 3;    % Minimum number of whisks at the start of a basal respiration cycle to be analyzed
maxWhisksBasalRespToAnalyze = 7;    % Maximum number of whisks at the start of a basal respiration cycle to be analyzed
nCorrToAnalyze = 4;                 % Number of whisk amplitude correlations to analyze

% Hard-coded strings in file names to exclude from averaging
excludeStringsFromAverage = {'ammpuff', 'airpuff', 'baseline', 'eth'};

% Plotting parameters
%fileNumsToPlot = 10;                    % The file number(s) to plot (max 38)
%fileNumsToPlot = 38;                    % The file number(s) to plot (max 38)
fileNumsToPlot = 4;                    % The file number(s) to plot (max 38)
%fileNumsToPlot = 13;                    % The file number(s) to plot (max 38)
%fileNumsToPlot = (14:38)';               % The file number(s) to plot (max 38)
%fileNumsToPlot = (1:38)';               % The file number(s) to plot (max 38)
toSpeedUp = false;                      % Whether to use parpool and hide figures
%toSpeedUp = true;                      % Whether to use parpool and hide figures
whiskAngleLimits = [-75, 75];           % Whisk angle limits to be plotted
piezoLimits = [-1, 10];
colorWhisk = [0, 0, 1];                 % Color for whisk trace (Blue)
colorResp = [1, 0, 0];                  % Color for resp trace (Red)
colorStim = [0, 0, 0];                  % Color for stim trace (Black)
colorAmmoniaPuff = [0.6, 0.8, 0.2];     % Color for Ammonia Puff (Yellow Green)
colorAirPuff = 0.5 * [1, 1, 1];         % Color for Air Puff (Gray)
colorSniffStartWin = [0.5, 1.0, 0.8];   % Color for sniff start windows (Aquamarine)
colorBasalRespCycle = [1.0, 0.8, 0.5];  % Color for basal respiration cycles (Light Orange)
faceAlphaSniffStartWin = 0.8;           % Transparencies for sniff start windows
faceAlphaBasalRespCycle = 0.8;          % Transparencies for basal respiration cycles
markerRespPeaksValleys = 'o';           % Marker for resp peaks and valleys (circle)
colorRespPeaksValleys = [1, 0.7, 0];    % Color for resp peaks and valleys (Orange)
markerWhiskPeaksValleys = 'x';          % Marker for whisk peaks and valleys (cross)
colorWhiskPeaksValleys = [0, 0.8, 0.8]; % Color for whisk peaks and valleys (Cyan)
lineStyleWhiskAmplitudes = '-';         % Line Style for whisk peak amplitudes (solid line)
colorWhiskAmplitudes = [0, 0.8, 0];     % Color for whisk peak amplitudes (Green)
colorSniffStart = [1, 0.7, 0];          % Color for sniff start transition lines (Orange)
colorSniffEnd = [0.6, 0, 0.8];          % Color for sniff end transition lines (Dark Violet)
lineWidthForSample = 0.5;               % Line width for sample traces plots
lineWidthForAnalysis = 1;               % Line width for sniff start window plots
markerSizeForSample = 6;                % Marker size for sample traces plots
markerSizeForAnalysis = 12;             % Marker size for sniff start window plots
colorT0 = [0, 0, 0.5];                  % Color for T0 horizontal bar (Dark Blue)
colorT1 = [0, 1, 1];                    % Color for T1 horizontal bar (Cyan)
barHeight = 5;                          % Height of horizontal bars in degrees
barYLow = -50;                          % Y-position of the bars in degrees
barFaceAlpha = 0.8;                     % Transparency for the filled bars
jitterWidth = 0.5;
markerSizeJitter = 4;
colorScatter = 'b';
markerTypeScatter = '.';
markerSizeScatter = 4;
markerLineWidthScatter = 0.5;
colorBestFit = [];              % Use default: if significant 'r', otherwise 'k'
lineStyleBestFit = '--';
lineWidthBestFit = 1;
textLocBestFit = 'topleft';
colorThrOrig = [1, 1, 1] * 0.5; % Gray
lineStyleThrOrig = '--';
lineWidthThrOrig = 0.5;
textLocThrOrig = 'bottomright';
whiskAngleLabel = 'Whisk Angle (degrees)';
respToWhiskRangeRatio = 0.8;            % Sniff amplitude to whisk amplitude ratios for plotting
timeLabel = 'Time (s)';
whiskLabel = 'Whisking';
respLabel = 'Breathing';
legendLocation1 = 'suppress'; %'northeast';
legendLocation2 = 'suppress';
legendLocation3 = 'suppress';
legendLocation4 = 'suppress';
subplotOrder1 = 'twoCols';              % Plot subplots as two columns in all traces plot
centerPosition1 = [400, 200, 1000, 160];% Center position for center subplot in all traces plot
figTitlePrefix1 = 'Whisking (blue) and breathing (red) data';
figTitlePrefix2 = 'Stimulation (Puff) data';
figTitlePrefix3 = 'Sniff Start Windows';
figTitlePrefix4 = 'Basal Respiration Cycles';
figTitle5 = 'Sniff Start Whisk Logarithmic Decrements';
figTitle6 = 'Basal Respiration Whisk Logarithmic Decrements';
figTitle7 = 'Sniff Start Whisk Amplitude Correlations';
figTitle8 = 'Basal Respiration Whisk Amplitude Correlations';
figTitle9 = ['Whisk Phase Response Curve with breath onset latency of ', num2str(breathOnsetLatencyMs), 'ms'];
figPrefix1 = 'sniff_whisk_all_traces_';
figPrefix2 = 'sniff_whisk_all_stims_';
figPrefix3 = 'sniff_whisk_sniffstart_windows_';
figPrefix4 = 'sniff_whisk_basalresp_cycles_';
figName5 = 'sniffstart_whisk_log_decrements_jitter';
figName6 = 'basalresp_whisk_log_decrements_jitter';
figName7 = 'sniffstart_whisk_amplitudes_scatter';
figName8 = 'basalresp_whisk_amplitudes_scatter';
figName9 = 'phase_response_scatter';
figTypes = {'png'}; % {'eps', 'png'};

% Create a structure for aggregate plotting parameters
pPlot.jitterWidth = jitterWidth;
pPlot.markerSizeJitter = markerSizeJitter;
pPlot.colorScatter = colorScatter;
pPlot.markerTypeScatter = markerTypeScatter;
pPlot.markerSizeScatter = markerSizeScatter;
pPlot.markerLineWidthScatter = markerLineWidthScatter;
pPlot.colorBestFit = colorBestFit;
pPlot.lineStyleBestFit = lineStyleBestFit;
pPlot.lineWidthBestFit = lineWidthBestFit;
pPlot.textLocBestFit = textLocBestFit;
pPlot.colorThrOrig = colorThrOrig;
pPlot.lineStyleThrOrig = lineStyleThrOrig;
pPlot.lineWidthThrOrig = lineWidthThrOrig;
pPlot.textLocThrOrig = textLocThrOrig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Get current directory
pathParentDir = pwd;

% Check if the data directory exists
% Construct the full path to the data directory
pathDataDir = fullfile(pathParentDir, nameDataDir);

% Check if the data directory exists. If not, display an error and stop.
fprintf('Checking for data directory: %s\n', pathDataDir);
if ~exist(pathDataDir, 'dir')
    fprintf('Error: Data directory not found.\n');
    fprintf('Please ensure the directory "%s" exists in the current path.\n', nameDataDir);
    return; % Stop the script
end
fprintf('Data directory found.\n\n');

% Create an output directory if it does not exist
pathOutDir = fullfile(pathParentDir, nameOutDir);
check_dir(pathOutDir);

% Open the log file for writing algorithm differences
pathAlgorithmDiffs = fullfile(pathOutDir, fileNameAlgorithmDiffs);
fileIDDiffs = fopen(pathAlgorithmDiffs, 'w');
if fileIDDiffs == -1
    error('Could not open log file for writing: %s', pathAlgorithmDiffs);
end

%% Find and Match Files
% Get a list of all data files in the data directory
fprintf('Searching for data and parameter files...\n');
[~, pathDataFiles] = ...
    all_files('Directory', pathDataDir, 'Extension', extDataFile, ...
                    'Prefix', prefixDataFile, 'Suffix', suffixDataFile, ...
                    'Sortby', 'date');

% Extract the file names
nameDataFiles = extract_fileparts(pathDataFiles, 'filebase');

% Extract the trial names by removing the prefix and suffix
trialNames = cellfun(@(x) erase(x, prefixDataFile), nameDataFiles, 'UniformOutput', false);
trialNames = cellfun(@(x) erase(x, suffixDataFile), trialNames, 'UniformOutput', false);

% Find the matching parameter files
[~, pathParamFiles] = ...
    find_matching_files(trialNames, 'Directory', pathDataDir, ...
                        'Extension', extDataFile, 'Suffix', suffixParamFile);
fprintf('\nFinished searching.\n\n');

% Count files
nFiles = numel(trialNames);

% Get all file numbers
fileNumbers = (1:nFiles)';

% Check if any matching pairs were found
if nFiles == 0
    fprintf('No matching data and parameter file pairs were found in the directory.\n');
    return
end

%% Load data
% Load data from each file: the data struct contains a single field
% 'sniffwhiskdata'
dataStructs = cellfun(@(x) load(x), pathDataFiles);
  
% Extract the sniffwhiskdata into a cell array
%   Note: Cannot concatenate into a structure array since not all
%   structures have therm and piezo fields
sniffWhiskDataCell = extract_fields(dataStructs, 'sniffwhiskdata');

% Extract the time vectors
tVecs = extract_fields(sniffWhiskDataCell, 't');

% Extract the whisk angle vectors from the camera detection
whiskVecs = extract_fields(sniffWhiskDataCell, 'svangle');

% Extract the sniff/breath vectors from the thermocouple channel
thermVecs = extract_fields(sniffWhiskDataCell, 'therm');

% Extract the external air pulses (puffs of ammonia or control puffs of unscented air to the nose)
pulseVecs = extract_fields(sniffWhiskDataCell, 'piezo');

% Count the number of vectors
nSweepsEachFile = count_vectors(tVecs);

% Check whether is ammonia puff or air puff trials
isAmmoniaPuff = contains(trialNames, ammoniaPuffString);

% Check whether is ammonia puff or air puff trials
isAirPuff = contains(trialNames, airPuffString);

%% Reformat data for ease of plotting
% Invert whisk angle vectors and thermocouple vectors
if invertThermSignal
    thermVecs = cellfun(@(x) -x, thermVecs, 'UniformOutput', false);
end
if invertWhiskSignal
    whiskVecs = cellfun(@(x) -x, whiskVecs, 'UniformOutput', false);
end

% Find the maximum and minimum whisk angles
maxWhiskAngle = apply_iteratively(@max, whiskVecs); % == 71.9
minWhiskAngle = apply_iteratively(@min, whiskVecs); % == -72.5

% Find the maximum and minimum thermocouple values for each file
maxThermValueEachFile = cellfun(@(x) apply_iteratively(@max, x), thermVecs);
minThermValueEachFile = cellfun(@(x) apply_iteratively(@min, x), thermVecs);

% Calculate plotting limits for thermocouple vectors for each file
meanThermLimitsEachFile = (maxThermValueEachFile + minThermValueEachFile) / 2;
rangeThermLimitsEachFile = (maxThermValueEachFile - minThermValueEachFile) / respToWhiskRangeRatio;
lowerThermLimitEachFile = meanThermLimitsEachFile - rangeThermLimitsEachFile / 2;

% Create sniff/breath vectors scaled to the order of whisk angle vectors
rangeWhiskAngleLimits = diff(whiskAngleLimits);
respVecs = cellfun(@(a, b, c) whiskAngleLimits(1) + (a - b) * ...
                    rangeWhiskAngleLimits / c, thermVecs, ...
                    num2cell(lowerThermLimitEachFile), num2cell(rangeThermLimitsEachFile), ...
                    'UniformOutput', false);

%% Detect the first pulse start and end times
pulseParams = cellfun(@(x, y) parse_pulse(x, 'TimeVecs', y, ...
        'MinPulseAmplitude', minPulseAmplitude), pulseVecs, tVecs, 'UniformOutput', false);
pulseStartTimes = cellfun(@(x) x.timeAfterStartMs, pulseParams, 'UniformOutput', false);
pulseEndTimes = cellfun(@(x) x.timeBeforeEndMs, pulseParams, 'UniformOutput', false);

%% Detect sniff transitions and basal respiration peaks
fprintf('Detecting sniff transitions and basal respiration peaks...\n');

% Use array_fun to apply the function to each file's data
[respPeakTablesAll, respValleyTablesAll, sniffStartTimesAll, sniffEndTimesAll, ...
    sniffFreqsFundamentalAll, basalRespPeakTablesAll] = ...
    array_fun(@(x, y, z) parse_resp_vecs(x, y, z, ...
                    amplitudeDefinition, fundFreqRange, fCutoffResp, fCutoffRelToFund, filterOrderResp, ...
                    promThresholdPercResp, minPeakDistanceMsResp, ...
                    sniffFreqThreshold, basalFreqThreshold), ...
                respVecs, tVecs, num2cell(nSweepsEachFile), 'UniformOutput', false);

fprintf('Finished detecting sniff transitions and basal respiration peaks.\n\n');

%% Detect whisk peaks and valleys and define analysis windows
fprintf('Detecting whisk peaks and valleys and defining analysis windows ...\n');
[whiskPeakTablesAll, whiskValleyTablesAll, whiskFreqsFundamentalAll, ...
    sniffStartWinTablesAll, basalRespCycleTablesAll] = ...
    array_fun(@(a, b, c, d, e, f, g, h) parse_whisk_vecs(a, b, c, d, e, f, g, h, ...
                    amplitudeDefinition, whiskDirForPhase, fundFreqRange, fCutoffWhisk, fCutoffRelToFund, filterOrderWhisk, ...
                    promThresholdPercWhisk, minPeakDistanceMsWhisk, minPeakPromWhisk, maxWhiskDurationMs, ...
                    breathOnsetLatencyMs, nWhisksSniffStartToAnalyze, minWhisksBasalRespToAnalyze, ...
                    sniffFreqThreshold, basalFreqThreshold, maxWhisksBasalRespToAnalyze, fileIDDiffs), ...
                num2cell(fileNumbers), trialNames, whiskVecs, tVecs, num2cell(nSweepsEachFile), ...
                sniffStartTimesAll, sniffEndTimesAll, basalRespPeakTablesAll, ...
                'UniformOutput', false, 'UseParpool', false);
fprintf('Finished detecting whisk peaks and valleys and defining analysis windows.\n\n');

%% Augment Sniff Start Windows with Sniff Data
fprintf('Augmenting sniff start windows with sniff data...\n');
sniffStartWinTablesAll = ...
    array_fun(@(a, b, c) augment_sniffstart_windows(a, b, c), ...
                sniffStartWinTablesAll, respPeakTablesAll, respValleyTablesAll, ...
                'UniformOutput', false);
fprintf('Finished augmenting sniff start windows.\n\n');

%% Plot all data
[handles1, handles2, handles3, handles4] = ...
    array_fun(@(a) plot_one_file (a, trialNames, tVecs, whiskVecs, respVecs, ...
        pulseVecs, pulseStartTimes, pulseEndTimes, ...
        respPeakTablesAll, respValleyTablesAll, sniffStartTimesAll, sniffEndTimesAll, ...
        whiskPeakTablesAll, whiskValleyTablesAll, sniffStartWinTablesAll, basalRespCycleTablesAll, ...
        isAmmoniaPuff, isAirPuff, whiskLabel, respLabel, ...
        colorWhisk, colorResp, colorStim, colorAmmoniaPuff, colorAirPuff, ...
        colorSniffStartWin, colorBasalRespCycle, faceAlphaSniffStartWin, faceAlphaBasalRespCycle, ...
        markerRespPeaksValleys, colorRespPeaksValleys, ...
        markerWhiskPeaksValleys, colorWhiskPeaksValleys, ...
        lineStyleWhiskAmplitudes, colorWhiskAmplitudes, ...
        colorSniffStart, colorSniffEnd, ...
        lineWidthForSample, lineWidthForAnalysis, ...
        markerSizeForSample, markerSizeForAnalysis, ...
        colorT0, colorT1, barHeight, barYLow, barFaceAlpha, ...
        pathOutDir, timeLabel, whiskAngleLimits, whiskAngleLabel, piezoLimits, ...
        legendLocation1, legendLocation2, legendLocation3, legendLocation4, ...
        subplotOrder1, centerPosition1, ...
        figTitlePrefix1, figTitlePrefix2, figTitlePrefix3, figTitlePrefix4, ...
        figPrefix1, figPrefix2, figPrefix3, figPrefix4, figTypes, toSpeedUp), ...
        fileNumsToPlot, 'UseParpool', toSpeedUp);

%% Calculate whisk logarithmic decrement statistics per file
fprintf('Calculating statistics for whisk logarithmic decrements for each file ...\n');

% Calculate the number of analysis windows per file
nSniffStartWindowsPerFile = cellfun(@height, sniffStartWinTablesAll);
nBasalRespCyclesPerFile = cellfun(@height, basalRespCycleTablesAll);

% Calculate the average whisk logarithmic decrements per file
[meanSniffWhiskLogDecrementsPerFile, stderrSniffWhiskLogDecrementsPerFile, ...
    lower95SniffWhiskLogDecrementsPerFile, upper95SniffWhiskLogDecrementsPerFile] = ...
    array_fun(@(x) compute_stats_for_cellnumeric(x.whiskLogDecrements), ...
            sniffStartWinTablesAll, 'UniformOutput', false);
[meanBasalWhiskLogDecrementsPerFile, stderrBasalWhiskLogDecrementsPerFile, ...
    lower95BasalWhiskLogDecrementsPerFile, upper95BasalWhiskLogDecrementsPerFile] = ...
    array_fun(@(x) compute_stats_for_cellnumeric(x.whiskLogDecrements), ...
            basalRespCycleTablesAll, 'UniformOutput', false);

fprintf('Finished calculating statistics for each file.\n\n');

%% Analyze and plot aggregate data
% Identify files to use for averaging
[isUsedForAverage, fileNumsToAverage] = ...
    identify_files_for_averaging(trialNames, excludeStringsFromAverage);

% Combine sniff start windows from selected files
sniffStartWinTableToAverage = combine_tables(sniffStartWinTablesAll, fileNumsToAverage);

% Combine basal respiration cycles from selected files
basalRespCycleTableToAverage = combine_tables(basalRespCycleTablesAll, fileNumsToAverage);

% Average the whisk logarithmic decrements in sniff start windows and plot as a grouped jitter plot
fprintf('Plotting aggregated sniff-start whisk log decrements jitter plot...\n');
[sniffWhiskLogDecrementResults, handles5] = ...
    virt_plot_log_decrement_jitter(sniffStartWinTableToAverage, pPlot, ...
                                'GroupingColumn', 'fileNumber', ...
                                'OutDir', pathOutDir, 'FigTypes', figTypes, ...
                                'FigTitle', figTitle5, 'FigName', figName5);
fprintf('Finished plotting aggregated sniff-start whisk log decrements jitter plot with mean/CI overlay.\n\n');

% Average the whisk logarithmic decrements in basal respiration cycles and plot as a grouped jitter plot
fprintf('Plotting aggregated basal-cycle whisk log decrements jitter plot...\n');
[basalWhiskLogDecrementResults, handles6] = ...
    virt_plot_log_decrement_jitter(basalRespCycleTableToAverage, pPlot, ...
                                'GroupingColumn', 'fileNumber', ...
                                'MaxDecrements', maxWhisksBasalRespToAnalyze, ...
                                'OutDir', pathOutDir, 'FigTypes', figTypes, ...
                                'FigTitle', figTitle6, 'FigName', figName6);
fprintf('Finished plotting aggregated basal-cycle whisk log decrements jitter plot with mean/CI overlay.\n\n');

% Plot successive whisk amplitudes in sniff start windows against each
%   other and compute the correlation coefficients
[sniffWhiskAmpCorrelationResults, handles7] = ...
    virt_plot_amplitude_correlation(sniffStartWinTableToAverage, pPlot, ...
                                'GroupingColumn', 'fileNumber', ...
                                'NCorrelations', nCorrToAnalyze, ...
                                'OutDir', pathOutDir, 'FigTypes', figTypes, ...
                                'FigTitle', figTitle7, 'FigName', figName7);

% Plot successive whisk amplitudes in basal respiration cycles against each
%   other and compute the correlation coefficients
[basalWhiskAmpCorrelationResults, handles8] = ...
    virt_plot_amplitude_correlation(basalRespCycleTableToAverage, pPlot, ...
                                'GroupingColumn', 'fileNumber', ...
                                'NCorrelations', nCorrToAnalyze, ...
                                'OutDir', pathOutDir, 'FigTypes', figTypes, ...
                                'FigTitle', figTitle8, 'FigName', figName8);

% Plot the phase response curve from aggregate data
fprintf('Plotting aggregated phase response curve...\n');
[phaseResponseResults, handles9] = ...
    virt_plot_phase_response(basalRespCycleTableToAverage, pPlot, ...
                            'WhiskDir', whiskDirForPhase, ...
                            'GroupingColumn', 'fileNumber', ...
                            'OutDir', pathOutDir, 'FigTypes', figTypes, ...
                            'FigTitle', figTitle9, 'FigName', figName9);
fprintf('Finished plotting aggregated phase response curve.\n\n');

%% Put all metadata together

% Load parameters from each file
paramStructs = cellfun(@load, pathParamFiles);
paramTable = struct2table(paramStructs);

% Create a table for all meta data
metaDataFiles = table(trialNames, pathDataFiles, pathParamFiles, ...
                 'VariableNames', {'TrialName', 'DataFilePath', 'ParameterFilePath'});
metaDataTable = horzcat(metaDataFiles, paramTable);

% Add detection results
metaDataTable.isUsedForAverage = isUsedForAverage;
metaDataTable.nSweeps = nSweepsEachFile;
metaDataTable.isAmmoniaPuff = isAmmoniaPuff;
metaDataTable.isAirPuff = isAirPuff;
metaDataTable.sniffFreqFundamental = sniffFreqsFundamentalAll;
metaDataTable.whiskFreqFundamental = whiskFreqsFundamentalAll;
metaDataTable.sniffStartTimes = sniffStartTimesAll;
metaDataTable.sniffEndTimes = sniffEndTimesAll;
metaDataTable.nSniffStartWindows = nSniffStartWindowsPerFile;
metaDataTable.meanSniffWhiskLogDecrements = meanSniffWhiskLogDecrementsPerFile;
metaDataTable.stderrSniffWhiskLogDecrements = stderrSniffWhiskLogDecrementsPerFile;
metaDataTable.lower95SniffWhiskLogDecrements = lower95SniffWhiskLogDecrementsPerFile;
metaDataTable.upper95SniffWhiskLogDecrements = upper95SniffWhiskLogDecrementsPerFile;
metaDataTable.nBasalRespCycles = nBasalRespCyclesPerFile;
metaDataTable.meanBasalWhiskLogDecrements = meanBasalWhiskLogDecrementsPerFile;
metaDataTable.stderrBasalWhiskLogDecrements = stderrBasalWhiskLogDecrementsPerFile;
metaDataTable.lower95BasalWhiskLogDecrements = lower95BasalWhiskLogDecrementsPerFile;
metaDataTable.upper95BasalWhiskLogDecrements = upper95BasalWhiskLogDecrementsPerFile;

%% Save results
% Save metadata table
pathMetaData = fullfile(pathOutDir, fileNameMetaData);
metaDataToPrint = write_table(metaDataTable, pathMetaData);
fprintf('Metadata table saved to %s!\n', pathMetaData);

% Save sniff start window info
pathSniffStartWinTable = fullfile(pathOutDir, fileNameSniffStartWinTable);
sniffStartWinTableToPrint = write_table(sniffStartWinTableToAverage, pathSniffStartWinTable);
fprintf('Sniff Start Window table saved to %s!\n', pathSniffStartWinTable);

% Save basal respiration cycle info
pathBasalRespCycleTable = fullfile(pathOutDir, fileNameBasalRespCycleTable);
basalRespCycleTableToPrint = write_table(basalRespCycleTableToAverage, pathBasalRespCycleTable);
fprintf('Basal Respiration Cycle table saved to %s!\n', pathBasalRespCycleTable);

% Save analysis results in a mat file
pathMatFile = fullfile(pathOutDir, fileNameAnalysisResults);
save(pathMatFile, 'metaDataTable', 'sniffStartWinTableToAverage', ...
     'basalRespCycleTableToAverage', 'sniffWhiskLogDecrementResults', ...
     'basalWhiskLogDecrementResults', 'respPeakTablesAll', ...
     'respValleyTablesAll', 'whiskPeakTablesAll', 'whiskValleyTablesAll');
fprintf('Analysis results saved to %s!\n', pathMatFile);

% Close the log file
fclose(fileIDDiffs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [handles1, handles2, handles3, handles4] = ...
    plot_one_file (sampleFileNum, trialNames, tVecs, whiskVecs, respVecs, ...
        pulseVecs, pulseStartTimes, pulseEndTimes, ...
        respPeakTablesAll, respValleyTablesAll, sniffStartTimesAll, sniffEndTimesAll, ...
        whiskPeakTablesAll, whiskValleyTablesAll, sniffStartWinTablesAll, basalRespCycleTablesAll, ...
        isAmmoniaPuff, isAirPuff, whiskLabel, respLabel, ...
        colorWhisk, colorResp, colorStim, colorAmmoniaPuff, colorAirPuff, ...
        colorSniffStartWin, colorBasalRespCycle, faceAlphaSniffStartWin, faceAlphaBasalRespCycle, ...
        markerRespPeaksValleys, colorRespPeaksValleys, ...
        markerWhiskPeaksValleys, colorWhiskPeaksValleys, ...
        lineStyleWhiskAmplitudes, colorWhiskAmplitudes, ...
        colorSniffStart, colorSniffEnd, ...
        lineWidthForSample, lineWidthForAnalysis, ...
        markerSizeForSample, markerSizeForAnalysis, ...
        colorT0, colorT1, barHeight, barYLow, barFaceAlpha, ...
        pathOutDir, timeLabel, whiskAngleLimits, whiskAngleLabel, piezoLimits, ...
        legendLocation1, legendLocation2, legendLocation3, legendLocation4, ...
        subplotOrder1, centerPosition1, ...
        figTitlePrefix1, figTitlePrefix2, figTitlePrefix3, figTitlePrefix4, ...
        figPrefix1, figPrefix2, figPrefix3, figPrefix4, figTypes, toSpeedUp)

% Extract data to plot
trialName = trialNames{sampleFileNum};
tVecsThis = tVecs{sampleFileNum};
whiskVecsThis = whiskVecs{sampleFileNum};
respVecsThis = respVecs{sampleFileNum};
pulseVecsThis = pulseVecs{sampleFileNum};
pulseStartTimesThis = pulseStartTimes{sampleFileNum};
pulseEndTimesThis = pulseEndTimes{sampleFileNum};
isAmmoniaPuffThis = isAmmoniaPuff(sampleFileNum);
isAirPuffThis = isAirPuff(sampleFileNum);
respPeakTablesThis = respPeakTablesAll{sampleFileNum};
respValleyTablesThis = respValleyTablesAll{sampleFileNum};
sniffStartTimesThis = sniffStartTimesAll{sampleFileNum};
sniffEndTimesThis = sniffEndTimesAll{sampleFileNum};
whiskPeakTablesThis = whiskPeakTablesAll{sampleFileNum};
whiskValleyTablesThis = whiskValleyTablesAll{sampleFileNum};
sniffStartWinTableForFile = sniffStartWinTablesAll{sampleFileNum};
basalRespCycleTableForFile = basalRespCycleTablesAll{sampleFileNum};

% Count the number of sweeps to plot
nSweeps = size(tVecsThis, 2);
recordNumbers = (1:nSweeps)';

% Create labels for legend
whiskLabels = create_labels_from_numbers(recordNumbers, 'Prefix', whiskLabel);
sniffLabels = create_labels_from_numbers(recordNumbers, 'Prefix', respLabel);

% Decide on puff window colors and titles
if isAmmoniaPuffThis
    figTitlePrefix1 = replace(figTitlePrefix1, 'data', 'data with Ammonia Puffs (Yellow Green)');
    figTitlePrefix2 = replace(figTitlePrefix2, 'Puff', 'Ammonia Puff');
    puffWindowColor = colorAmmoniaPuff;
elseif isAirPuffThis
    figTitlePrefix1 = replace(figTitlePrefix1, 'data', 'data with Control Air Puffs (gray)');
    figTitlePrefix2 = replace(figTitlePrefix2, 'Puff', 'Control Air Puff');
    puffWindowColor = colorAirPuff; 
else
    puffWindowColor = colorAirPuff; 
end

% Create figure titles
figTitle1 = [figTitlePrefix1, ' for File # ', num2str(sampleFileNum)];
figTitle2 = [figTitlePrefix2, ' for File # ', num2str(sampleFileNum)];
figTitle3 = [figTitlePrefix3, ' for File # ', num2str(sampleFileNum)];
figTitle4 = [figTitlePrefix4, ' for File # ', num2str(sampleFileNum)];

% Create paths for saving
figName1 = [figPrefix1, 'File', num2str(sampleFileNum), '_', trialName];
figName2 = [figPrefix2, 'File', num2str(sampleFileNum), '_', trialName];
figName3 = [figPrefix3, 'File', num2str(sampleFileNum), '_', trialName];
figName4 = [figPrefix4, 'File', num2str(sampleFileNum), '_', trialName];
figPath1 = fullfile(pathOutDir, figName1);
figPath2 = fullfile(pathOutDir, figName2);
figPath3 = fullfile(pathOutDir, figName3);
figPath4 = fullfile(pathOutDir, figName4);

% Create all traces plot
fprintf('Plotting all traces for file %d ...\n', sampleFileNum);
if toSpeedUp
    set_figure_properties('AlwaysNew', true, 'Visible', 'off');
else
    set_figure_properties('FigNumber', 1, 'AlwaysNew', false, 'ClearFigure', true);
end
if ~isempty(respVecsThis)
    handles1 = plot_traces(tVecsThis, whiskVecsThis, 'ColorMap', colorWhisk, ...
                'DataToCompare', respVecsThis, 'ColorMapToCompare', colorResp, ...
                'XBoundaries', [pulseStartTimesThis, pulseEndTimesThis], ...
                'XBoundaryColor', puffWindowColor, ...
                'TraceLabels', whiskLabels, 'TraceLabelsToCompare', sniffLabels, ...
                'PlotMode', 'parallel', 'FigTitle', figTitle1, ...
                'TightInset', true, 'FigExpansion', 'auto', ...
                'SubPlotOrder', subplotOrder1, 'CenterPosition', centerPosition1, ...
                'XLabel', timeLabel, 'LegendLocation', legendLocation1, ...
                'YLimits', whiskAngleLimits, 'YLabel', whiskAngleLabel, ...
                'LineWidth', lineWidthForSample);

    % Overlay detected peaks, valleys, and transitions
    % Get all subplots
    allSubPlots = handles1.subPlots;

    % Plot annotations for each sweep
    for iSwp = 1:nSweeps
        % Hold on to current subplot
        subplot(allSubPlots(iSwp));
        hold on;

        % Get the sniff start windows for this specific sweep
        if ~isempty(sniffStartWinTableForFile)
            sniffStartWinTableForSweep = ...
                sniffStartWinTableForFile(sniffStartWinTableForFile.sweepNumber == iSwp, :);
        else
            sniffStartWinTableForSweep = [];
        end

        % Plot sniff start windows if they exist
        if ~isempty(sniffStartWinTableForSweep)
            % Get all window boundaries for this sweep
            winBoundaries = [sniffStartWinTableForSweep.sniffStartWinStartTime, ...
                             sniffStartWinTableForSweep.sniffStartWinEndTime];

            % Plot all shades for this sweep
            plot_vertical_shade(winBoundaries', 'Color', colorSniffStartWin, 'FaceAlpha', faceAlphaSniffStartWin);
        end

        % Get the basal respiration cycles for this specific sweep
        if ~isempty(basalRespCycleTableForFile)
            basalRespCycleTableForSweep = ...
                basalRespCycleTableForFile(basalRespCycleTableForFile.sweepNumber == iSwp, :);
        else
            basalRespCycleTableForSweep = [];
        end

        % Plot basal respiration cycles if they exist
        if ~isempty(basalRespCycleTableForSweep)
            % Get all window boundaries for this sweep
            winBoundaries = [basalRespCycleTableForSweep.basalRespCycleStartTime, ...
                             basalRespCycleTableForSweep.basalRespCycleEndTime];

            % Plot all shades for this sweep
            plot_vertical_shade(winBoundaries', 'Color', colorBasalRespCycle, 'FaceAlpha', faceAlphaBasalRespCycle);
        end

        % Extract the peak and valley tables for this sweep
        respPeakTable = respPeakTablesThis{iSwp};
        respValleyTable = respValleyTablesThis{iSwp};
        whiskPeakTable = whiskPeakTablesThis{iSwp};
        whiskValleyTable = whiskValleyTablesThis{iSwp};
        
        % Plot peaks if they exist
        if ~isempty(respPeakTable) && ismember('peakIndex', respPeakTable.Properties.VariableNames)
            plot(respPeakTable.peakTime, respPeakTable.peakValue, ...
                markerRespPeaksValleys, 'Color', colorRespPeaksValleys, ...
                'MarkerSize', markerSizeForSample, 'LineWidth', lineWidthForSample);
        end
        if ~isempty(whiskPeakTable) && ismember('peakIndex', whiskPeakTable.Properties.VariableNames)
            plot(whiskPeakTable.peakTime, whiskPeakTable.peakValue, ...
                markerWhiskPeaksValleys, 'Color', colorWhiskPeaksValleys, ...
                'MarkerSize', markerSizeForSample, 'LineWidth', lineWidthForSample);
        end
     
        % Plot valleys if they exist
        if ~isempty(respValleyTable) && ismember('valleyIndex', respValleyTable.Properties.VariableNames)
            plot(respValleyTable.valleyTime, respValleyTable.valleyValue, ...
                markerRespPeaksValleys, 'Color', colorRespPeaksValleys, ...
                'MarkerSize', markerSizeForSample, 'LineWidth', lineWidthForSample);
        end
        if ~isempty(whiskValleyTable) && ismember('valleyIndex', whiskValleyTable.Properties.VariableNames)
            plot(whiskValleyTable.valleyTime, whiskValleyTable.valleyValue, ...
                markerWhiskPeaksValleys, 'Color', colorWhiskPeaksValleys, ...
                'MarkerSize', markerSizeForSample, 'LineWidth', lineWidthForSample);
        end
        
        % Get the transition times for this sweep
        sniffStartTimesThisSweep = sniffStartTimesThis{iSwp};
        sniffEndTimesThisSweep = sniffEndTimesThis{iSwp};

        % Plot basal-to-sniff transitions
        if ~isempty(sniffStartTimesThisSweep)
            plot_vertical_line(sniffStartTimesThisSweep, 'Color', colorSniffStart, ...
                'LineStyle', '--', 'LineWidth', lineWidthForSample);
        end

        % Plot sniff-to-basal transitions
        if ~isempty(sniffEndTimesThisSweep)
            plot_vertical_line(sniffEndTimesThisSweep, 'Color', colorSniffEnd, ...
                'LineStyle', '--', 'LineWidth', lineWidthForSample);
        end
    end

    % Save figure
    save_all_figtypes(handles1.fig, figPath1, figTypes);
else
    handles1 = plot_traces(tVecsThis, whiskVecsThis, 'ColorMap', colorWhisk, ...
                'TraceLabels', whiskLabels, ...
                'XBoundaries', [pulseStartTimesThis, pulseEndTimesThis], ...
                'XBoundaryColor', puffWindowColor, ...
                'PlotMode', 'parallel', 'FigTitle', figTitle1, ...
                'XLabel', timeLabel, 'LegendLocation', legendLocation1, ...
                'YLimits', whiskAngleLimits, 'YLabel', whiskAngleLabel, ...
                'FigName', figPath1, 'FigTypes', figTypes, ...
                'LineWidth', lineWidthForSample);
end

% Plot stim traces with detection
if ~isempty(pulseVecsThis)
    fprintf('Plotting stim traces with detection for file %d ...\n', sampleFileNum);
    if toSpeedUp
        set_figure_properties('AlwaysNew', true, 'Visible', 'off');
    else
        set_figure_properties('FigNumber', 2, 'AlwaysNew', false, 'ClearFigure', true);
    end
    handles2 = plot_traces(tVecsThis, pulseVecsThis, 'ColorMap', colorStim, ...
                'XBoundaries', [pulseStartTimesThis, pulseEndTimesThis], ...
                'XBoundaryColor', puffWindowColor, ...
                'PlotMode', 'parallel', 'FigTitle', figTitle2, ...
                'XLabel', timeLabel, 'LegendLocation', legendLocation2, ...
                'YLimits', piezoLimits, ...
                'FigName', figPath2, 'FigTypes', figTypes);
else
    handles2.fig = gobjects;
    handles2.subPlots = gobjects;
    handles2.plotsData = gobjects;
    handles2.plotsDataToCompare = gobjects;
end

% Plot each sniff start window in a separate subplot
if ~isempty(sniffStartWinTableForFile)
    fprintf('Plotting sniff start windows for file %d ...\n', sampleFileNum);

    % Get the number of sniff start windows
    nWindows = height(sniffStartWinTableForFile);

    % Extract variables from the table
    sweepNums = sniffStartWinTableForFile.sweepNumber;
    winStarts = sniffStartWinTableForFile.sniffStartWinStartTime;
    winEnds = sniffStartWinTableForFile.sniffStartWinEndTime;
    sniffStarts = sniffStartWinTableForFile.sniffStartTime;
    respPeakTimesForWin = sniffStartWinTableForFile.respPeakTimes;
    respPeakValuesForWin = sniffStartWinTableForFile.respPeakValues;
    respValleyTimesForWin = sniffStartWinTableForFile.respValleyTimes;
    respValleyValuesForWin = sniffStartWinTableForFile.respValleyValues;
    whiskPeakTimesForWin = sniffStartWinTableForFile.whiskPeakTimes;
    whiskPeakValuesForWin = sniffStartWinTableForFile.whiskPeakValues;
    whiskPeakAmplitudesForWin = sniffStartWinTableForFile.whiskPeakAmplitudes;
    whiskValleyTimesForWin = sniffStartWinTableForFile.whiskValleyTimes;
    whiskValleyValuesForWin = sniffStartWinTableForFile.whiskValleyValues;
    whiskLogDecrementsForWin = sniffStartWinTableForFile.whiskLogDecrements;
    
    % Compute window durations
    windowDurations = winEnds - winStarts;

    % Compute maximum window duration
    maxWindowDuration = max(windowDurations);    

    % Compute window ends to plot
    winEndsToPlot = winStarts + maxWindowDuration;

    % Loop through each sniff start window to extract data segments
    tVecsForWin = cell(nWindows, 1);
    whiskVecsForWin = cell(nWindows, 1);
    respVecsForWin = cell(nWindows, 1);
    for iWin = 1:nWindows
        % Get info for the current window
        swpNum = sweepNums(iWin);
        winStart = winStarts(iWin);
        winEnd = winEnds(iWin);
        
        % Get the full traces for the corresponding sweep
        tFull = tVecsThis(:, swpNum);
        whiskFull = whiskVecsThis(:, swpNum);
        sniffFull = respVecsThis(:, swpNum);
        
        % Find indices within the current sniff start window
        indicesInWin = find(tFull >= winStart & tFull <= winEnd);
        
        % Use extract_subvectors to get the data segments
        tVecsForWin{iWin} = extract_subvectors(tFull, 'Indices', indicesInWin);
        whiskVecsForWin{iWin} = extract_subvectors(whiskFull, 'Indices', indicesInWin);
        respVecsForWin{iWin} = extract_subvectors(sniffFull, 'Indices', indicesInWin);
    end

    % Set figure properties for the new plot
    if toSpeedUp
        set_figure_properties('AlwaysNew', true, 'Visible', 'off');
    else
        set_figure_properties('FigNumber', 3, 'AlwaysNew', false, 'ClearFigure', true);
    end
    
    % Plot all extracted windows in parallel subplots
    handles3 = plot_traces(tVecsForWin, whiskVecsForWin, 'ColorMap', colorWhisk, ...
        'DataToCompare', respVecsForWin, 'ColorMapToCompare', colorResp, ...
        'PlotMode', 'parallel', 'FigTitle', figTitle3, ...
        'XLabel', timeLabel, 'LegendLocation', legendLocation3, ...
        'YLimits', whiskAngleLimits, 'YLabel', whiskAngleLabel, ...
        'LineWidth', lineWidthForAnalysis);

    % Overlay detections on each sniff start window subplot
    for iWin = 1:nWindows
        % Get info for the current window
        winStart = winStarts(iWin);
        winEndToPlot = winEndsToPlot(iWin);
        sniffStart = sniffStarts(iWin);
        respPeakTimes = respPeakTimesForWin{iWin};
        respPeakValues = respPeakValuesForWin{iWin};
        respValleyTimes = respValleyTimesForWin{iWin};
        respValleyValues = respValleyValuesForWin{iWin};
        whiskPeakTimes = whiskPeakTimesForWin{iWin};
        whiskPeakValues = whiskPeakValuesForWin{iWin};
        whiskPeakAmplitudes = whiskPeakAmplitudesForWin{iWin};
        whiskValleyTimes = whiskValleyTimesForWin{iWin};
        whiskValleyValues = whiskValleyValuesForWin{iWin};
        whiskLogDecrements = whiskLogDecrementsForWin{iWin};

        % Select the correct subplot
        subplot(handles3.subPlots(iWin));
        hold on;

        % Restrict x-axis limits
        xlim([winStart, winEndToPlot]);

        % Plot sniff start times
        if ~isnan(sniffStart)
            plot_vertical_line(sniffStart, 'Color', colorSniffStart, ...
                'LineWidth', lineWidthForAnalysis, 'LineStyle', '--');
        end
        
        % Plot resp peaks used for analysis
        if ~isempty(respPeakTimes)
            plot(respPeakTimes, respPeakValues, ...
                markerRespPeaksValleys, 'Color', colorRespPeaksValleys, ...
                'LineWidth', lineWidthForAnalysis, 'MarkerSize', markerSizeForAnalysis);
        end

        % Plot resp valleys used for analysis
        if ~isempty(respValleyTimes)
            plot(respValleyTimes, respValleyValues, ...
                markerRespPeaksValleys, 'Color', colorRespPeaksValleys, ...
                'LineWidth', lineWidthForAnalysis, 'MarkerSize', markerSizeForAnalysis);
        end
        
        % Plot whisk peaks used for analysis
        if ~isempty(whiskPeakTimes)
            plot(whiskPeakTimes, whiskPeakValues, ...
                markerWhiskPeaksValleys, 'Color', colorWhiskPeaksValleys, ...
                'LineWidth', lineWidthForAnalysis, 'MarkerSize', markerSizeForAnalysis);

            % Prepare coordinate matrices for vectorized plotting
            % Each column represents a line: [x_start; x_end], [y_start; y_end]
            xCoords = [whiskPeakTimes'; whiskPeakTimes'];
            yCoords = [whiskPeakValues'; whiskPeakValues' - whiskPeakAmplitudes'];
            
            % Plot whisk amplitudes
            plot(xCoords, yCoords, lineStyleWhiskAmplitudes, 'Color', colorWhiskAmplitudes);
        end

        % Plot whisk valleys used for analysis
        if ~isempty(whiskValleyTimes)
            plot(whiskValleyTimes, whiskValleyValues, ...
                markerWhiskPeaksValleys, 'Color', colorWhiskPeaksValleys, ...
                'LineWidth', lineWidthForAnalysis, 'MarkerSize', markerSizeForAnalysis);
        end

        % Display the successive logarithmic decrements (deltas)
        if ~isempty(whiskLogDecrements)
            axLimits = axis;
            xPos = axLimits(1) + 0.05 * (axLimits(2) - axLimits(1));
            yPos = axLimits(4) - 0.1 * (axLimits(4) - axLimits(3));
            text(xPos, yPos, ['Deltas:', sprintf(' %.2f', whiskLogDecrements)]);
        end
    end
  
    % Save the figure
    save_all_figtypes(handles3.fig, figPath3, figTypes);
else
    % If there are no sniff start windows, create empty handles
    handles3.fig = gobjects;
    handles3.subPlots = gobjects;
    handles3.plotsData = gobjects;
    handles3.plotsDataToCompare = gobjects;
end

% Plot each basal respiration cycle in a separate subplot
if ~isempty(basalRespCycleTableForFile)
    fprintf('Plotting basal respiration cycles for file %d ...\n', sampleFileNum);

    % Get the number of cycles
    nCycles = height(basalRespCycleTableForFile);
    
    % Extract variables from the table
    sweepNums = basalRespCycleTableForFile.sweepNumber;
    cycleStarts = basalRespCycleTableForFile.basalRespCycleStartTime;
    cycleEnds = basalRespCycleTableForFile.basalRespCycleEndTime;
    respPeakTimesForCycle = basalRespCycleTableForFile.basalRespPeakTime;
    respPeakValuesForCycle = basalRespCycleTableForFile.basalRespPeakValue;
    respPreValleyTimesForCycle = basalRespCycleTableForFile.basalRespPreValleyTime;
    respPostValleyTimesForCycle = basalRespCycleTableForFile.basalRespPostValleyTime;
    whiskPeakTimesForCycle = basalRespCycleTableForFile.whiskPeakTimes;
    whiskPeakValuesForCycle = basalRespCycleTableForFile.whiskPeakValues;
    whiskPeakAmplitudesForCycle = basalRespCycleTableForFile.whiskPeakAmplitudes;
    whiskValleyTimesForCycle = basalRespCycleTableForFile.whiskValleyTimes;
    whiskValleyValuesForCycle = basalRespCycleTableForFile.whiskValleyValues;
    whiskLogDecrementsForCycle = basalRespCycleTableForFile.whiskLogDecrements;
    breathOnsetTimesForCycle = basalRespCycleTableForFile.breathOnsetTime;
    eventTimeWhiskBeforeForCycle = basalRespCycleTableForFile.eventTimeWhiskBefore;
    eventTimeTwoWhisksBeforeForCycle = basalRespCycleTableForFile.eventTimeTwoWhisksBefore;
    preIEIWhiskAfterForCycle = basalRespCycleTableForFile.preIEIWhiskAfter;
    
    % Compute window starts to plot
    %   This is the earlier of two whisk peaks before respiration peak,
    %       the pre-valley of the inspiratory whisk, 
    %       or the respiration onset
    winStartsToPlot = min([eventTimeTwoWhisksBeforeForCycle, cycleStarts, respPreValleyTimesForCycle], [], 2);

    % Compute required window durations
    windowDurations = cycleEnds - winStartsToPlot;

    % Compute maximum window duration required
    maxWinDuration = max(windowDurations);    

    % Compute window ends to plot
    winEndsToPlot = winStartsToPlot + maxWinDuration;

    % Loop through each cycle to extract data segments
    tVecsForCycle = cell(nCycles, 1);
    whiskVecsForCycle = cell(nCycles, 1);
    respVecsForCycle = cell(nCycles, 1);
    for iCycle = 1:nCycles
        % Get info for the current cycle
        swpNum = sweepNums(iCycle);
        winStart = winStartsToPlot(iCycle);
        winEnd = winEndsToPlot(iCycle);
        
        % Get the full traces for the corresponding sweep
        tFull = tVecsThis(:, swpNum);
        whiskFull = whiskVecsThis(:, swpNum);
        sniffFull = respVecsThis(:, swpNum);
        
        % Find indices within the current cycle
        indicesInWin = find(tFull >= winStart & tFull <= winEnd);
        
        % Use extract_subvectors to get the data segments
        tVecsForCycle{iCycle} = extract_subvectors(tFull, 'Indices', indicesInWin);
        whiskVecsForCycle{iCycle} = extract_subvectors(whiskFull, 'Indices', indicesInWin);
        respVecsForCycle{iCycle} = extract_subvectors(sniffFull, 'Indices', indicesInWin);
    end

    % Set figure properties for the new plot
    if toSpeedUp
        set_figure_properties('AlwaysNew', true, 'Visible', 'off');
    else
        set_figure_properties('FigNumber', 4, 'AlwaysNew', false, 'ClearFigure', true);
    end
    
    % Plot all extracted cycles in parallel subplots
    handles4 = plot_traces(tVecsForCycle, whiskVecsForCycle, 'ColorMap', colorWhisk, ...
        'DataToCompare', respVecsForCycle, 'ColorMapToCompare', colorResp, ...
        'PlotMode', 'parallel', 'FigTitle', figTitle4, ...
        'XLabel', timeLabel, 'LegendLocation', legendLocation4, ...
        'YLimits', whiskAngleLimits, 'YLabel', whiskAngleLabel, ...
        'LineWidth', lineWidthForAnalysis);

    % Overlay detections on each basal respiration cycle subplot
    for iCycle = 1:nCycles
        % Get info for the current cycle
        winStartToPlot = winStartsToPlot(iCycle);
        winEndToPlot = winEndsToPlot(iCycle);
        respPeakTime = respPeakTimesForCycle(iCycle);
        respPeakValue = respPeakValuesForCycle(iCycle);
        respPreValleyTime = respPreValleyTimesForCycle(iCycle);
        respPostValleyTime = respPostValleyTimesForCycle(iCycle);
        whiskPeakTimes = whiskPeakTimesForCycle{iCycle};
        whiskPeakValues = whiskPeakValuesForCycle{iCycle};
        whiskPeakAmplitudes = whiskPeakAmplitudesForCycle{iCycle};
        whiskValleyTimes = whiskValleyTimesForCycle{iCycle};
        whiskValleyValues = whiskValleyValuesForCycle{iCycle};
        whiskLogDecrements = whiskLogDecrementsForCycle{iCycle};
        breathOnsetTime = breathOnsetTimesForCycle(iCycle);
        eventTimeWhiskBefore = eventTimeWhiskBeforeForCycle(iCycle);
        eventTimeTwoWhisksBefore = eventTimeTwoWhisksBeforeForCycle(iCycle);
        preIEIWhiskAfter = preIEIWhiskAfterForCycle(iCycle);

        % Select the correct subplot
        ax = handles4.subPlots(iCycle);
        subplot(ax);
        hold on;

        % Restrict x-axis limits
        xlim([winStartToPlot, winEndToPlot]);
        
        % Plot phase response detection with horizontal bars
        % Check if the necessary data for phase response is available (not NaN)
        if ~isnan(eventTimeWhiskBefore) && ~isnan(eventTimeTwoWhisksBefore) && ...
           ~isnan(preIEIWhiskAfter) && ~isnan(breathOnsetTime)
            % Set the y coordinates for the horizontal bars
            yT0T1Rect = barYLow + barHeight * [-1, 1, 1, -1] / 2;
            yResetRect = barYLow + barHeight + barHeight * [-1, 1, 1, -1] / 2;

            % 1. Plot preIEIWhiskBefore with colorT0
            % This is the unperturbed inter-whisk-interval (T0)
            xT0Rect = [eventTimeTwoWhisksBefore, eventTimeTwoWhisksBefore, eventTimeWhiskBefore, eventTimeWhiskBefore];
            fill(xT0Rect, yT0T1Rect, colorT0, 'FaceAlpha', barFaceAlpha);
            
            % 2. Plot preIEIWhiskAfter with colorT1
            % This is the perturbed inter-whisk-interval (T1)
            peakTimeAfter = eventTimeWhiskBefore + preIEIWhiskAfter;
            xT1Rect = [eventTimeWhiskBefore, eventTimeWhiskBefore, peakTimeAfter, peakTimeAfter];
            fill(xT1Rect, yT0T1Rect, colorT1, 'FaceAlpha', barFaceAlpha);

            % 3. Plot relativeResetTime with colorResp
            % This represents the timing of the breath onset within the unperturbed cycle
            xResetRect = [eventTimeWhiskBefore, eventTimeWhiskBefore, breathOnsetTime, breathOnsetTime];
            fill(xResetRect, yResetRect, colorResp, 'FaceAlpha', barFaceAlpha);
        end

        % Plot resp peaks and valley times used for analysis
        if ~isempty(respPeakTime)
            plot(respPeakTime, respPeakValue, ...
                markerRespPeaksValleys, 'Color', colorRespPeaksValleys, ...
                'LineWidth', lineWidthForAnalysis, 'MarkerSize', markerSizeForAnalysis);
            plot_vertical_line(respPreValleyTime, 'Color', colorRespPeaksValleys, ...
                'LineStyle', '--', ...
                'LineWidth', lineWidthForAnalysis, 'MarkerSize', markerSizeForAnalysis);
            plot_vertical_line(respPostValleyTime, 'Color', colorRespPeaksValleys, ...
                'LineStyle', '--', ...
                'LineWidth', lineWidthForAnalysis, 'MarkerSize', markerSizeForAnalysis);
        end
       
        % Plot whisk peaks used for analysis
        if ~isempty(whiskPeakTimes)
            plot(whiskPeakTimes, whiskPeakValues, ...
                markerWhiskPeaksValleys, 'Color', colorWhiskPeaksValleys, ...
                'LineWidth', lineWidthForAnalysis, 'MarkerSize', markerSizeForAnalysis);

            % Prepare coordinate matrices for vectorized plotting
            xCoords = [whiskPeakTimes'; whiskPeakTimes'];
            yCoords = [whiskPeakValues'; whiskPeakValues' - whiskPeakAmplitudes'];
            
            % Plot whisk amplitudes
            plot(xCoords, yCoords, lineStyleWhiskAmplitudes, 'Color', colorWhiskAmplitudes);
        end

        % Plot whisk valleys used for analysis
        if ~isempty(whiskValleyTimes)
            plot(whiskValleyTimes, whiskValleyValues, ...
                markerWhiskPeaksValleys, 'Color', colorWhiskPeaksValleys, ...
                'LineWidth', lineWidthForAnalysis, 'MarkerSize', markerSizeForAnalysis);
        end

        % Display the successive logarithmic decrements (deltas)
        if ~isempty(whiskLogDecrements)
            axLimits = axis;
            xPos = axLimits(1) + 0.05 * (axLimits(2) - axLimits(1));
            yPos = axLimits(4) - 0.1 * (axLimits(4) - axLimits(3));
            text(xPos, yPos, ['Deltas:', sprintf(' %.2f', whiskLogDecrements)]);
        end
    end

    % Save the figure
    save_all_figtypes(handles4.fig, figPath4, figTypes);
else
    % If there are no basal respiration cycles, create empty handles
    handles4.fig = gobjects;
    handles4.subPlots = gobjects;
    handles4.plotsData = gobjects;
    handles4.plotsDataToCompare = gobjects;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [respPeakTables, valleyTables, sniffStartTimes, sniffEndTimes, sniffFreqFundamental, basalRespPeakTable] = ...
    parse_resp_vecs(respVecsThisFile, tVecsThisFile, nSweeps, ...
                    amplitudeDefinition, fundFreqRange, fCutoff, fCutoffRelToFund, filterOrder, ...
                    promThresholdPerc, minPeakDistanceMs, ...
                    sniffFreqThreshold, basalFreqThreshold)
%% Parses Sniff vectors for each file

% Hard-coded parameters for this function's logic
minRespPeaksForBasal = 2;      % Minimum number of peaks to detect basal respiration
minRespPeaksForTransition = 3; % Minimum number of peaks to detect transitions

% Compute inter peak interval thresholds
sniffIpiThresholdSec = 1 / sniffFreqThreshold;
basalIpiThresholdSec = 1 / basalFreqThreshold;

% Initialize cell arrays to store results for each sweep in this file.
respPeakTables = cell(nSweeps, 1);
valleyTables = cell(nSweeps, 1);
sniffStartTimes = cell(nSweeps, 1);
sniffEndTimes = cell(nSweeps, 1);
sniffFreqFundamental = cell(nSweeps, 1);
basalRespPeakTablesPerSweep = cell(nSweeps, 1);

% If sniff vector is empty, return
if isempty(respVecsThisFile)
    basalRespPeakTable = table();
    return
end

% Loop through each sweep within the current file.
for iSwp = 1:nSweeps
    % Extract the sniff vector for the current sweep.
    respVec = respVecsThisFile(:, iSwp);
    % Extract the time vector for the current sweep.
    timeVec = tVecsThisFile(:, iSwp);

    % Find all peak times and preceding valley times for the current sniff vector.
    [respPeakTable, respValleyTable, otherResults] = ...
        parse_oscillation(respVec, ...
                        'TimeVec', timeVec, 'TimeUnits', 's', ...
                        'AmpMode', amplitudeDefinition, ...
                        'FundFreqRange', fundFreqRange, ...
                        'FilterCutoffs', fCutoff, ...
                        'FilterCutoffsRelToFund', fCutoffRelToFund, ...
                        'FilterOrder', filterOrder, ...
                        'PromThresholdPerc', promThresholdPerc, ...
                        'MinPeakDistanceMs', minPeakDistanceMs);

    % Store the tables
    respPeakTables{iSwp} = respPeakTable;
    valleyTables{iSwp} = respValleyTable;

    % Store the fundamental frequency if it was computed
    if isfield(otherResults, 'freqFundamental')
        sniffFreqFundamental{iSwp} = otherResults.freqFundamental;
    else
        sniffFreqFundamental{iSwp} = NaN;
    end

    % If fewer than the minimum required peaks are found, we can't compute IPIs, so we skip.
    if height(respPeakTable) < minRespPeaksForBasal
        sniffStartTimes{iSwp} = []; % Store an empty result for this sweep.
        sniffEndTimes{iSwp} = []; % Store an empty result for this sweep.
        basalRespPeakTablesPerSweep{iSwp} = table(); % Store an empty table for this sweep.
        continue; % Move to the next sweep.
    end

    % Extract the peak times from the resulting table.
    peakTimes = respPeakTable.peakTime;

    % Compute the inter-peak intervals (IPIs) between all resp peaks.
    interPeakIntervals = diff(peakTimes);

    % Detect basal respiration peaks. A basal respiration peak is a resp trace peak
    % with a succeeding inter-peak interval >= basalIpiThresholdSec.
    succeedingIPIs = interPeakIntervals; % for peaks 1 to end-1
    isBasalRespPeak = succeedingIPIs >= basalIpiThresholdSec;
    respPeakNumber = find(isBasalRespPeak);
    basalRespSucceedingIPIs = succeedingIPIs(isBasalRespPeak);

    if ~isempty(respPeakNumber)
        % Get the corresponding peak and valley info from respPeakTable
        basalRespPeaksData = respPeakTable(respPeakNumber, :);

        % Count the number of basal respiration peaks
        nBasalRespPeaks = height(basalRespPeaksData);

        % Create the basalRespPeakTable for this sweep
        sweepNumber = repmat(iSwp, nBasalRespPeaks, 1);
        basalRespPeakNumber = (1:nBasalRespPeaks)';
        basalRespPeakTimes = basalRespPeaksData.peakTime;
        basalRespPeakValues = basalRespPeaksData.peakValue;
        basalRespPreValleyTimes = basalRespPeaksData.preValleyTime;
        basalRespPostValleyTimes = basalRespPeaksData.postValleyTime;

        % Determine whether each basal respiration peak is 
        %   preceded by a basal respiration peak
        isWithinBasal = false(nBasalRespPeaks, 1);
        if nBasalRespPeaks > 1
            % The current peak number must be exact 1 more than the previous peak number
            isWithinBasal(2:end) = (diff(respPeakNumber) == 1);
        end

        basalRespPeakTablesPerSweep{iSwp} = ...
            table(sweepNumber, respPeakNumber, basalRespPeakNumber, ...
                  basalRespPeakTimes, basalRespPeakValues, basalRespPreValleyTimes, ...
                  basalRespPostValleyTimes, basalRespSucceedingIPIs, ...
                  isWithinBasal);
    else
        basalRespPeakTablesPerSweep{iSwp} = table(); % No basal peaks found
    end

    % If fewer than the minimum required peaks are found, we can't compute IPIs, so we skip.
    if height(respPeakTable) < minRespPeaksForTransition
        sniffStartTimes{iSwp} = []; % Store an empty result for this sweep.
        sniffEndTimes{iSwp} = []; % Store an empty result for this sweep.
        continue; % Move to the next sweep.
    end

    % Extract the times of the valleys that precede each peak.
    preValleyTimes = respPeakTable.preValleyTime;

    % Remove first and last peak from peaks to test
    preValleyTimesToTest = preValleyTimes(2:end-1);

    % Find the preceding inter-peak intervals for each peak 2:end
    preIPIs = interPeakIntervals(1:end-1);

    % Find the succeeding inter-peak intervals for each peak 2:end
    postIPIs = interPeakIntervals(2:end);

    % Detect start of sniffing transition times.
    % This is defined as a long IPI (basal) followed by a short IPI (sniffing).
    isSniffStart = preIPIs > sniffIpiThresholdSec & postIPIs <= sniffIpiThresholdSec;
    sniffStartTimes{iSwp} = preValleyTimesToTest(isSniffStart);

    % Detect end of sniffing transition times.
    % This is defined as a short IPI (sniffing) followed by a long IPI (basal).
    isSniffEnd = preIPIs <= sniffIpiThresholdSec & postIPIs > sniffIpiThresholdSec;
    sniffEndTimes{iSwp} = preValleyTimesToTest(isSniffEnd);
end

% Vertically concatenate all sweep tables for this file into a single table
basalRespPeakTable = vertcat(basalRespPeakTablesPerSweep{:});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [whiskPeakTables, whiskValleyTables, whiskFreqFundamental, sniffStartWinTable, basalRespCycleTable] = ...
    parse_whisk_vecs(fileNumber, trialName, whiskVecsThisFile, tVecsThisFile, nSweeps, ...
                    sniffStartTimesThisFile, sniffEndTimesThisFile, basalRespPeakTableThisFile, ...
                    amplitudeDefinition, whiskDirForPhase, fundFreqRange, fCutoff, fCutoffRelToFund, filterOrder, ...
                    promThresholdPerc, minPeakDistanceMs, minPeakProm, maxWhiskDurationMs, ...
                    breathOnsetLatencyMs, nWhisksSniffStartToAnalyze, minWhisksBasalRespToAnalyze, ...
                    sniffFreqThreshold, basalFreqThreshold, maxWhisksBasalRespToAnalyze, fileIDDiffs)
%% Parses whisk vectors for a single file and generates sniff start window and basal resp cycle tables

% Initialize cell arrays to store results for each sweep in this file.
whiskPeakTables = cell(nSweeps, 1);
whiskValleyTables = cell(nSweeps, 1);
whiskFreqFundamental = cell(nSweeps, 1);

% Detect whisk peaks in parallel
parfor iSwp = 1:nSweeps
    % Find all peak times and preceding valley times for the current whisk vector.
    [whiskPeakTable, whiskValleyTable, otherResults] = ...
        parse_oscillation(whiskVecsThisFile(:, iSwp), ...
                        'TimeVec', tVecsThisFile(:, iSwp), 'TimeUnits', 's', ...
                        'AmpMode', amplitudeDefinition, ...
                        'FundFreqRange', fundFreqRange, ...
                        'FilterCutoffs', fCutoff, ...
                        'FilterCutoffsRelToFund', fCutoffRelToFund, ...
                        'FilterOrder', filterOrder, ...
                        'PromThresholdPerc', promThresholdPerc, ...
                        'MinPeakProminence', minPeakProm, ...
                        'MinValleyProminence', minPeakProm, ...
                        'MinPeakDistanceMs', minPeakDistanceMs, ...
                        'MaxDurationForAmplitudeMs', maxWhiskDurationMs);
    
    % Store the tables
    whiskPeakTables{iSwp} = whiskPeakTable;
    whiskValleyTables{iSwp} = whiskValleyTable;
    
    % Store the fundamental frequency if it was computed
    if isfield(otherResults, 'freqFundamental')
        whiskFreqFundamental{iSwp} = otherResults.freqFundamental;
    else
        whiskFreqFundamental{iSwp} = NaN;
    end
end

% Package analysis parameters into a structure
analysisParams.nWhisksSniffStartToAnalyze = nWhisksSniffStartToAnalyze;
analysisParams.sniffFreqThreshold = sniffFreqThreshold;
analysisParams.basalFreqThreshold = basalFreqThreshold;
analysisParams.minWhisksBasalRespToAnalyze = minWhisksBasalRespToAnalyze;
analysisParams.maxWhisksBasalRespToAnalyze = maxWhisksBasalRespToAnalyze;
analysisParams.whiskDirForPhase = whiskDirForPhase;
analysisParams.breathOnsetLatencyMs = breathOnsetLatencyMs;

% Detect Sniff and Basal Analysis Windows using the modular function
% Note: Pass empty tables for peakTable and valleyTable as they are not used
%       in 'transitionTimes' or 'respPeak' modes.
[sniffStartWinTable, basalRespCycleTable] = ...
    virt_detect_whisk_analysis_windows(table(), table(), analysisParams, ...
        'SniffDefinitionMode', 'transitionTimes', ...
        'BasalDefinitionMode', 'respPeak', ...
        'TVecs', tVecsThisFile, ...
        'SniffStartTimes', sniffStartTimesThisFile, ...
        'SniffEndTimes', sniffEndTimesThisFile, ...
        'NSweeps', nSweeps, ...
        'WhiskPeakTables', whiskPeakTables, ...
        'WhiskValleyTables', whiskValleyTables, ...
        'BasalRespPeakTable', basalRespPeakTableThisFile, ...
        'FileNumber', fileNumber, ...
        'TrialName', trialName, ...
        'FileIDDiffs', fileIDDiffs, ...
        'ToPostMessage', false); % Suppress messages as this is in a loop

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sniffStartWinTable = augment_sniffstart_windows(sniffStartWinTable, ...
                                        respPeakTablesThisFile, ...
                                        respValleyTablesThisFile)
%% Augments a single sniff start window table with sniff data for one file

% If no sniff start windows were found for this file, return the empty table
if isempty(sniffStartWinTable)
    % Define empty columns for the case where the input table is empty
    newCols = {'respPeakTimes', 'respPeakValues', 'respPeakAmplitudes', ...
                'respPreValleyTimes', 'respPostValleyTimes', ...
                'respValleyTimes', 'respValleyValues'};
    for i = 1:numel(newCols)
        sniffStartWinTable.(newCols{i}) = cell(0, 1);
    end

    return;
end

% Get the number of sniff start windows
nWindows = height(sniffStartWinTable);

% Initialize new columns as cell arrays
respPeakTimesInWin = cell(nWindows, 1);
respPeakValuesInWin = cell(nWindows, 1);
respPeakAmpsInWin = cell(nWindows, 1);
respPreValleyTimesInWin = cell(nWindows, 1);
respPostValleyTimesInWin = cell(nWindows, 1);
respValleyTimesInWin = cell(nWindows, 1);
respValleyValuesInWin = cell(nWindows, 1);

% Loop through each sniff start window in the table
for iWin = 1:nWindows
    % Get window start/end times and sweep number
    winStart = sniffStartWinTable.sniffStartWinStartTime(iWin);
    winEnd = sniffStartWinTable.sniffStartWinEndTime(iWin);
    swpNum = sniffStartWinTable.sweepNumber(iWin);

    % Get the sniff data for the corresponding sweep
    respPeakTable = respPeakTablesThisFile{swpNum};
    respValleyTable = respValleyTablesThisFile{swpNum};

    % Find resp peaks within the window
    if ~isempty(respPeakTable)
        isPeakInWin = respPeakTable.peakTime >= winStart & ...
                      respPeakTable.peakTime <= winEnd;
        respPeakTimesInWin{iWin} = respPeakTable.peakTime(isPeakInWin);
        respPeakValuesInWin{iWin} = respPeakTable.peakValue(isPeakInWin);
        respPeakAmpsInWin{iWin} = respPeakTable.amplitude(isPeakInWin);
        respPreValleyTimesInWin{iWin} = respPeakTable.preValleyTime(isPeakInWin);
        respPostValleyTimesInWin{iWin} = respPeakTable.postValleyTime(isPeakInWin);
    end

    % Find resp valleys within the window
    if ~isempty(respValleyTable)
        isValleyInWin = respValleyTable.valleyTime >= winStart & ...
                        respValleyTable.valleyTime <= winEnd;
        respValleyTimesInWin{iWin} = respValleyTable.valleyTime(isValleyInWin);
        respValleyValuesInWin{iWin} = respValleyTable.valleyValue(isValleyInWin);
    end
end

% Add the new data as columns to the table
sniffStartWinTable.respPeakTimes = respPeakTimesInWin;
sniffStartWinTable.respPeakValues = respPeakValuesInWin;
sniffStartWinTable.respPeakAmplitudes = respPeakAmpsInWin;
sniffStartWinTable.respPreValleyTimes = respPreValleyTimesInWin;
sniffStartWinTable.respPostValleyTimes = respPostValleyTimesInWin;
sniffStartWinTable.respValleyTimes = respValleyTimesInWin;
sniffStartWinTable.respValleyValues = respValleyValuesInWin;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [means, stderrs, lower95s, upper95s] = compute_stats_for_cellnumeric (vecs)
%% Computes the statistics for a cell array of numeric vectors

means = compute_combined_trace(vecs, 'mean');
stderrs = compute_combined_trace(vecs, 'stderr');
lower95s = compute_combined_trace(vecs, 'lower95');
upper95s = compute_combined_trace(vecs, 'upper95');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [isUsedForAverage, fileNumsToAverage] = ...
                identify_files_for_averaging (trialNames, excludeStringsFromAverage)
%% Identifies files for averaging by excluding specific trial names.

% Initialize as all true
isUsedForAverage = true(size(trialNames));

% Iteratively apply exclusion criteria
for i = 1:numel(excludeStringsFromAverage)
    isUsedForAverage = isUsedForAverage & ~contains(trialNames, excludeStringsFromAverage{i});
end

% Find the file numbers (indices) to be used for averaging
fileNumsToAverage = find(isUsedForAverage);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combinedTable = combine_tables (tablesAll, tableNumsToCombine)
%% Combines selected tables into one large table.

% Select the tables to be combined
tablesToCombine = tablesAll(tableNumsToCombine);

% Vertically concatenate all tables in the list into a single table
combinedTable = vertcat(tablesToCombine{:});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Expand the cell arrays to a single string
metaDataToPrint = metaDataTable;
[metaDataToPrint.roiswitch, metaDataToPrint.records, ...
    metaDataToPrint.anesrecords, metaDataToPrint.sniffStartTimes, ...
    metaDataToPrint.sniffEndTimes, metaDataToPrint.sniffFreqFundamental, ...
    metaDataToPrint.whiskFreqFundamental, ...
    metaDataTable.meanWhiskLogDecrements, metaDataTable.stderrWhiskLogDecrements, ...
    metaDataTable.lower95WhiskLogDecrements, metaDataTable.upper95WhiskLogDecrements] = ...
    argfun(@(x) cellfun(@(a) print_cellstr(a, 'Delimiter', ' ', ...
                            'ToPrint', false, 'OmitQuotes', true, ...
                            'OmitBraces', true, 'OmitNewline', true), ...
                            x, 'UniformOutput', false), ...
        metaDataTable.roiswitch, metaDataTable.records, ...
        metaDataTable.anesrecords, metaDataTable.sniffStartTimes, ...
        metaDataTable.sniffEndTimes, metaDataTable.sniffFreqFundamental, ...
        metaDataTable.whiskFreqFundamental, ...
        metaDataTable.meanWhiskLogDecrements, metaDataTable.stderrWhiskLogDecrements, ...
        metaDataTable.lower95WhiskLogDecrements, metaDataTable.upper95WhiskLogDecrements);

% Display the resulting table in the command window
fprintf('Generated Metadata Table:\n');
disp(metaDataToPrint);

fCutoffRelToFund = [0.1, 10];  % ratio of Butterworth bandpass filter cutoff for whisk trace
                                % relative to the fundamental frequency

% Remove any empty tables from the cell array
tablesToCombine(cellfun('isempty', tablesToCombine)) = [];

% Vertically concatenate all tables in the list into a single table
if ~isempty(tablesToCombine)
    combinedTable = vertcat(tablesToCombine{:});
else
    % Return an empty table if no data was found
    combinedTable = table();
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

