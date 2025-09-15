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
%       cd/create_subplots.m
%       cd/create_labels_from_numbers.m
%       cd/extract_fields.m
%       cd/extract_fileparts.m
%       cd/extract_subvectors.m
%       cd/find_matching_files.m
%       cd/force_matrix.m
%       cd/parse_oscillation.m
%       cd/plot_correlation_coefficient.m
%       cd/plot_grouped_jitter.m
%       cd/plot_grouped_scatter.m
%       cd/plot_test_result.m
%       cd/plot_traces.m
%       cd/plot_vertical_line.m
%       cd/plot_vertical_shade.m
%       cd/print_cellstr.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/resize_subplots_for_labels.m
%       cd/test_difference.m
%       cd/vecfun.m
%       cd/write_table.m
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
amplitudeDefinition = 'peak-to-avgvalley';
fundFreqRange = [0.5, 20];  % range of possible fundamental frequencies to be detected
fCutoffResp = [1, 15];      % bandpass filter cutoff for resp trace (Moore et al 2013 used [1, 15] Hz)
fCutoffWhisk = [3, 25];     % bandpass filter cutoff for whisk trace (Moore et al 2013 used [3, 25] Hz)
fCutoffRelToFund = [];      % don't use this
filterOrderResp = 3;        % Butterworth filter order for resp trace (Moore et al 2013 used 3)
filterOrderWhisk = 3;       % Butterworth filter order for whisk trace (Moore et al 2013 used 3)
promThresholdPercResp = 5;  % Percentage of amplitude range for minimum peak prominence for resp peaks
promThresholdPercWhisk = 5; % Percentage of amplitude range for minimum peak prominence for whisk peaks
minPeakPromWhisk = 5;       % Minimum whisk angle change (degrees) to detect as a peak (Moore et al 2013 used 5 degrees)
maxWhiskDurationMs = 250;   % Maximum whisk inter-valley interval (ms) (Moore et al 2013 used whisk duration < 250 ms)
minPeakDistanceMsResp = 30;     % Minimum peak distance (ms) for resp peaks
minPeakDistanceMsWhisk = 30;    % Minimum peak distance (ms) for whisk peaks
breathOnsetLatencyMs = 30;    % Presumed latency (ms) for from PB neuron activation to breath onset 
sniffFreqThreshold = 4;     % Frequency threshold for sniffing in Hz
basalFreqThreshold = 3;     % Frequency threshold for basal respiration in Hz
nWhisksSniffStartToAnalyze = 5; % Number of whisks at the start of a sniff period to be analyzed
minWhisksBasalRespToAnalyze = 3;  % Minimum number of whisks at the start of a basal respiration cycle to be analyzed
maxWhisksBasalRespToAnalyze = 7;  % Maximum number of whisks at the start of a basal respiration cycle to be analyzed
nCorrToAnalyze = 4;         % Number of whisk amplitude correlations to analyze

% Hard-coded strings in file names to exclude from averaging
excludeStringsFromAverage = {'ammpuff', 'airpuff', 'baseline', 'eth'};

% Plotting parameters
%fileNumsToPlot = 10;                    % The file number(s) to plot (max 38)
%fileNumsToPlot = 38;                    % The file number(s) to plot (max 38)
%fileNumsToPlot = 4;                    % The file number(s) to plot (max 38)
fileNumsToPlot = (1:38)';               % The file number(s) to plot (max 38)
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
markerTypeScatter = '.';
markerSizeScatter = 4;
markerLineWidthScatter = 0.5;
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
figPrefix1 = 'sniff_whisk_all_traces_';
figPrefix2 = 'sniff_whisk_all_stims_';
figPrefix3 = 'sniff_whisk_sniffstart_windows_';
figPrefix4 = 'sniff_whisk_basalresp_cycles_';
figName5 = 'sniffstart_whisk_log_decrements_jitter';
figName6 = 'basalresp_whisk_log_decrements_jitter';
figName7 = 'sniffstart_whisk_amplitudes_scatter';
figName8 = 'basalresp_whisk_amplitudes_scatter';
figTypes = {'png'}; % {'eps', 'png'};

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
thermVecs = cellfun(@(x) -x, thermVecs, 'UniformOutput', false);
whiskVecs = cellfun(@(x) -x, whiskVecs, 'UniformOutput', false);

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
    array_fun(@(a, b, c, d, e, f, g) parse_whisk_vecs(a, b, c, d, e, f, g, ...
                    fileIDDiffs, amplitudeDefinition, fundFreqRange, fCutoffWhisk, fCutoffRelToFund, filterOrderWhisk, ...
                    promThresholdPercWhisk, minPeakDistanceMsWhisk, minPeakPromWhisk, maxWhiskDurationMs, ...
                    breathOnsetLatencyMs, nWhisksSniffStartToAnalyze, minWhisksBasalRespToAnalyze), ...
                num2cell(fileNumbers), whiskVecs, tVecs, num2cell(nSweepsEachFile), ...
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
            sniffStartWinTablesAll, 'UniformOutput', false);

fprintf('Finished calculating statistics for each file.\n\n');

%% Analyze and plot aggregate data
% Identify files to use for averaging
[isUsedForAverage, fileNumsToAverage] = ...
    identify_files_for_averaging(trialNames, excludeStringsFromAverage);

% Combine sniff start windows from selected files
sniffStartWinTableToAverage = ...
    combine_across_files(sniffStartWinTablesAll, fileNumsToAverage, trialNames);

% Combine basal respiration cycles from selected files
basalRespCycleTableToAverage = ...
    combine_across_files(basalRespCycleTablesAll, fileNumsToAverage, trialNames);

% Average the whisk logarithmic decrements in sniff start windows and plot as a grouped jitter plot
fprintf('Plotting aggregated sniff-start whisk log decrements jitter plot...\n');
[sniffWhiskLogDecrementResults, handles5] = ...
    average_log_decrement_jitter(sniffStartWinTableToAverage, [], ...
                                pathOutDir, figTypes, figTitle5, figName5);
fprintf('Finished plotting aggregated sniff-start whisk log decrements jitter plot with mean/CI overlay.\n\n');

% Average the whisk logarithmic decrements in basal respiration cycles and plot as a grouped jitter plot
fprintf('Plotting aggregated basal-cycle whisk log decrements jitter plot...\n');
[basalWhiskLogDecrementResults, handles6] = ...
    average_log_decrement_jitter(basalRespCycleTableToAverage, maxWhisksBasalRespToAnalyze, ...
                                pathOutDir, figTypes, figTitle6, figName6);
fprintf('Finished plotting aggregated basal-cycle whisk log decrements jitter plot with mean/CI overlay.\n\n');

% Plot successive whisk amplitudes in sniff start windows against each
%   other and compute the correlation coefficients
[sniffWhiskAmpCorrelationResults, handles7] = ...
    correlate_successive_whiskamp(sniffStartWinTableToAverage, nCorrToAnalyze, ...
                                pathOutDir, figTypes, figTitle7, figName7, ...
                                markerTypeScatter, markerSizeScatter, markerLineWidthScatter);

% Plot successive whisk amplitudes in basal respiration cycles against each
%   other and compute the correlation coefficients
[basalWhiskAmpCorrelationResults, handles8] = ...
    correlate_successive_whiskamp(basalRespCycleTableToAverage, nCorrToAnalyze, ...
                                pathOutDir, figTypes, figTitle8, figName8, ...
                                markerTypeScatter, markerSizeScatter, markerLineWidthScatter);

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
    sweepNums = basalRespCycleTableForFile.sweepNumber;
    cycleStarts = basalRespCycleTableForFile.basalRespCycleStartTime;
    cycleEnds = basalRespCycleTableForFile.basalRespCycleEndTime;
    respPeakTimesForCycle = basalRespCycleTableForFile.basalRespPeakTime;
    respPreValleyTimesForCycle = basalRespCycleTableForFile.basalRespPreValleyTime;
    respPostValleyTimesForCycle = basalRespCycleTableForFile.basalRespPostValleyTime;
    whiskPeakTimesForCycle = basalRespCycleTableForFile.whiskPeakTimes;
    whiskPeakValuesForCycle = basalRespCycleTableForFile.whiskPeakValues;
    whiskPeakAmplitudesForCycle = basalRespCycleTableForFile.whiskPeakAmplitudes;
    whiskValleyTimesForCycle = basalRespCycleTableForFile.whiskValleyTimes;
    whiskValleyValuesForCycle = basalRespCycleTableForFile.whiskValleyValues;
    whiskLogDecrementsForCycle = basalRespCycleTableForFile.whiskLogDecrements;
    
    % Compute cycle durations
    cycleDurations = cycleEnds - cycleStarts;

    % Compute maximum cycle duration
    maxCycleDuration = max(cycleDurations);    

    % Compute cycle ends to plot
    winEndsToPlot = cycleStarts + maxCycleDuration;

    % Loop through each cycle to extract data segments
    tVecsForCycle = cell(nCycles, 1);
    whiskVecsForCycle = cell(nCycles, 1);
    respVecsForCycle = cell(nCycles, 1);
    for iCycle = 1:nCycles
        % Get info for the current cycle
        swpNum = sweepNums(iCycle);
        winStart = cycleStarts(iCycle);
        winEnd = cycleEnds(iCycle);
        
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
        winStart = cycleStarts(iCycle);
        winEndToPlot = winEndsToPlot(iCycle);
        respPeakTime = respPeakTimesForCycle(iCycle);
        respPreValleyTime = respPreValleyTimesForCycle(iCycle);
        respPostValleyTime = respPostValleyTimesForCycle(iCycle);
        whiskPeakTimes = whiskPeakTimesForCycle{iCycle};
        whiskPeakValues = whiskPeakValuesForCycle{iCycle};
        whiskPeakAmplitudes = whiskPeakAmplitudesForCycle{iCycle};
        whiskValleyTimes = whiskValleyTimesForCycle{iCycle};
        whiskValleyValues = whiskValleyValuesForCycle{iCycle};
        whiskLogDecrements = whiskLogDecrementsForCycle{iCycle};

        % Select the correct subplot
        subplot(handles4.subPlots(iCycle));
        hold on;

        % Restrict x-axis limits
        xlim([winStart, winEndToPlot]);
        
        % Plot resp peaks and valley times used for analysis
        if ~isempty(respPeakTime)
            plot_vertical_line(respPeakTime, 'Color', colorRespPeaksValleys, ...
                'LineStyle', '--', ...
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
    basalRespPeakIndices = find(isBasalRespPeak);
    basalRespSucceedingIPIs = succeedingIPIs(isBasalRespPeak);

    if ~isempty(basalRespPeakIndices)
        % Get the corresponding peak and valley info from respPeakTable
        basalRespPeaksData = respPeakTable(basalRespPeakIndices, :);

        % Create the basalRespPeakTable for this sweep
        sweepNumber = repmat(iSwp, height(basalRespPeaksData), 1);
        basalRespPeakNumber = (1:height(basalRespPeaksData))';
        basalRespPeakTimes = basalRespPeaksData.peakTime;
        basalRespPreValleyTimes = basalRespPeaksData.preValleyTime;
        basalRespPostValleyTimes = basalRespPeaksData.postValleyTime;

        basalRespPeakTablesPerSweep{iSwp} = table(sweepNumber, basalRespPeakNumber, ...
                                    basalRespPeakTimes, basalRespPreValleyTimes, ...
                                    basalRespPostValleyTimes, basalRespSucceedingIPIs);
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
    parse_whisk_vecs(fileNumber, whiskVecsThisFile, tVecsThisFile, nSweeps, ...
                    sniffStartTimesThisFile, sniffEndTimesThisFile, basalRespPeakTableThisFile, ...
                    fileID, amplitudeDefinition, fundFreqRange, fCutoff, fCutoffRelToFund, filterOrder, ...
                    promThresholdPerc, minPeakDistanceMs, minPeakProm, maxWhiskDurationMs, ...
                    breathOnsetLatencyMs, nWhisksSniffStartToAnalyze, minWhisksBasalRespToAnalyze)
%% Parses whisk vectors for a single file and generates sniff start window and basal resp cycle tables

% Convert to seconds
breathOnsetLatency = breathOnsetLatencyMs / 1000;

% Initialize cell arrays to store results for each sweep in this file.
whiskPeakTables = cell(nSweeps, 1);
whiskValleyTables = cell(nSweeps, 1);
whiskFreqFundamental = cell(nSweeps, 1);

% Define empty table structure for sniff start windows
emptySniffTable = table('Size', [0, 14], 'VariableTypes', ...
    {'double', 'double', 'double', 'double', 'double', 'double', ...
     'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'}, ...
    'VariableNames', {'sweepNumber', 'windowNumber', 'sniffStartWinStartTime', ...
                      'sniffStartWinEndTime', 'sniffStartTime', 'sniffEndTime', ...
                      'whiskPeakTimes', 'whiskPeakValues', 'whiskPeakAmplitudes', ...
                      'whiskPreValleyTimes', 'whiskPostValleyTimes', ...
                      'whiskValleyTimes', 'whiskValleyValues', 'whiskLogDecrements'});

% Define empty table structure for basal respiration cycles
emptyBasalTable = table('Size', [0, 17], 'VariableTypes', ...
    {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', ...
     'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'}, ...
    'VariableNames', {'sweepNumber', 'basalRespCycleNumber', 'basalRespPeakNumber', ...
                      'basalRespCycleStartTime', 'basalRespCycleEndTime', ...
                      'basalRespPeakTime', 'basalRespPreValleyTime', ...
                      'basalRespPostValleyTime', 'basalRespSucceedingIPI', ...
                      'whiskPeakTimes', 'whiskPeakValues', 'whiskPeakAmplitudes', ...
                      'whiskPreValleyTimes', 'whiskPostValleyTimes', ...
                      'whiskValleyTimes', 'whiskValleyValues', 'whiskLogDecrements'});

% If whisk vector is empty, return empty results
if isempty(whiskVecsThisFile)
    sniffStartWinTable = emptySniffTable;
    basalRespCycleTable = emptyBasalTable;
    return
end

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

% Now, process sweeps serially to build the sniff start window table
allWindowsCell = {};
for iSwp = 1:nSweeps
    % Obtain whisk peak and valley tables for this sweep
    whiskPeakTable = whiskPeakTables{iSwp};
    whiskValleyTable = whiskValleyTables{iSwp};

    % Obtain whisk peak times and sniff period start and end times
    peakTimes = whiskPeakTable.peakTime;
    preValleyTimes = whiskPeakTable.preValleyTime;
    sniffStartTimes = sniffStartTimesThisFile{iSwp};
    sniffEndTimes = sniffEndTimesThisFile{iSwp};
    valleyTimes = whiskValleyTable.valleyTime;

    % If there are no sniff periods or no whisks, skip this sweep
    if isempty(sniffStartTimes) || isempty(whiskPeakTable)
        continue;
    end
    
    % If the first sniff period end time is before
    %   the first sniff period start time, remove that first end time
    if ~isempty(sniffEndTimes) && sniffEndTimes(1) < sniffStartTimes (1)
        sniffEndTimes(1) = [];
    end

    % If there are more sniff period start times then end times, 
    %   add the end of time vector as an sniff end time
    if numel(sniffStartTimes) > numel(sniffEndTimes)
        sniffEndTimes(end+1) = tVecsThisFile(end, iSwp);
    end

    % Loop through all sniff periods for this sweep
    iWindow = 0;            % Counter for sniff start windows for this sweep
    for iSniff = 1:numel(sniffStartTimes)
        % Get current sniff period start and end times
        currentSniffStart = sniffStartTimes(iSniff);
        currentSniffEnd = sniffEndTimes(iSniff);

        % Find all peaks within the sniff period
        isPeakInSniffPeriod = peakTimes >= currentSniffStart & ...
                            peakTimes <= currentSniffEnd;

        % Moore et al 2013 definition: Find all peaks with prevalley times 
        %   after breathOnsetLatency before sniff start time
        %   but peaks before sniff end time
        isPreValleyInSniffPeriod = ...
            preValleyTimes >= currentSniffStart - breathOnsetLatency & ...
            peakTimes <= currentSniffEnd;

        % If these definitions do not match, print warning
        peakNumInSniffPeriod = find(isPeakInSniffPeriod);
        preValleyNumInSniffPeriod = find(isPreValleyInSniffPeriod);
        preValleyTooEarly = setdiff(peakNumInSniffPeriod, preValleyNumInSniffPeriod);
        peakToLate = setdiff(preValleyNumInSniffPeriod, peakNumInSniffPeriod);
        if ~isempty(preValleyTooEarly)
            fprintf(fileID, ['Whisk protraction for first whisk peak after sniff start occurred ', ...
                     'too early for sniffing period #%d in sweep %d of file %d ', ...
                     ' for peak numbers (for the sweep): %s!!\n\n'], ...
                     iSniff, iSwp, fileNumber, num2str(preValleyTooEarly));
        end

        % Define wether a whisk is in a sniff period
        isInSniffPeriod = isPeakInSniffPeriod;
%        isInSniffPeriod = isPreValleyInSniffPeriod;

        % Restrict whisk peak table to those in the sniff period
        whiskPeaksInSniffPeriod = whiskPeakTable(isInSniffPeriod, :);

        % If there are at least nWhisksSniffStartToAnalyze peaks
        %   within this sniff period and the first of those have all valid amplitudes, 
        %   add to sniff start window with appropriate boundaries
        if height(whiskPeaksInSniffPeriod) >= nWhisksSniffStartToAnalyze
            % Extract whisk peaks to analyze
            whiskPeaksToAnalyze = whiskPeaksInSniffPeriod(1:nWhisksSniffStartToAnalyze, :);

            % Skip this window if some amplitudes not valid
            if any(isnan(whiskPeaksToAnalyze.amplitude))
                continue;
            end

            % Increment window count
            iWindow = iWindow + 1;

            % Start of sniff start window is the pre-valley time of the 
            %   1st peak to analyze, or if doesn't exist, the sniff start time 
            firstPeakToAnalyze = whiskPeaksToAnalyze(1, :);
            firstPreValleyTime = firstPeakToAnalyze.preValleyTime;
            if ~isnan(firstPreValleyTime)
                sniffStartWinStartTime = firstPreValleyTime;
            else
                sniffStartWinStartTime = currentSniffStart;
            end

            % End of sniff start window is the post-valley time of the 
            %   last peak to analyze, or if doesn't exist, the sniff end time
            lastPeakToAnalyze = whiskPeaksToAnalyze(end, :);
            lastPostValleyTime = lastPeakToAnalyze.postValleyTime;
            if ~isnan(lastPostValleyTime)
                sniffStartWinEndTime = lastPostValleyTime;
            else
                sniffStartWinEndTime = currentSniffEnd;
            end

            % Find all valleys within the sniff start window
            isValleyInWin = valleyTimes >= sniffStartWinStartTime & ...
                          valleyTimes <= sniffStartWinEndTime;
            whiskValleysToAnalyze = whiskValleyTable(isValleyInWin, :);

            % Compute the logarithmic decrements of successive whisk peak amplitudes
            peakAmplitudes = whiskPeaksToAnalyze.amplitude;
            if numel(peakAmplitudes) > 1
                logDecrements = log(peakAmplitudes(2:end) ./ peakAmplitudes(1:end-1));
            else
                logDecrements = nan(size(nWhisksSniffStartToAnalyze - 1, 1));
            end
                        
            % Create a one-row table for this window
            newWindow = table(iSwp, iWindow, sniffStartWinStartTime, sniffStartWinEndTime, ...
                currentSniffStart, currentSniffEnd, ...
                {whiskPeaksToAnalyze.peakTime}, {whiskPeaksToAnalyze.peakValue}, ...
                {peakAmplitudes}, {whiskPeaksToAnalyze.preValleyTime}, {whiskPeaksToAnalyze.postValleyTime}, ...
                {whiskValleysToAnalyze.valleyTime}, {whiskValleysToAnalyze.valleyValue}, {logDecrements}, ...
                'VariableNames', {'sweepNumber', 'windowNumber', 'sniffStartWinStartTime', 'sniffStartWinEndTime', ...
                                  'sniffStartTime', 'sniffEndTime', ...
                                  'whiskPeakTimes', 'whiskPeakValues', 'whiskPeakAmplitudes', ...
                                  'whiskPreValleyTimes', 'whiskPostValleyTimes', ...
                                  'whiskValleyTimes', 'whiskValleyValues', 'whiskLogDecrements'});

            allWindowsCell{end+1} = newWindow;
        end
    end
end

% Vertically concatenate all found windows into a single table
if ~isempty(allWindowsCell)
    sniffStartWinTable = vertcat(allWindowsCell{:});
else
    sniffStartWinTable = emptySniffTable;
end

% Find basal respiration cycles
allCyclesCell = {};
for iSwp = 1:nSweeps
    % Filter the file-level basal peak table for the current sweep
    if ~isempty(basalRespPeakTableThisFile)
        basalRespPeaks = basalRespPeakTableThisFile(basalRespPeakTableThisFile.sweepNumber == iSwp, :);
    else
        basalRespPeaks = table();
    end
    
    % Extract whisk data for this sweep
    whiskPeaks = whiskPeakTables{iSwp};
    whiskValleys = whiskValleyTables{iSwp};
    whiskPeakTimesThisSwp = whiskPeaks.peakTime;
    whiskPreValleyTimesThisSwp = whiskPeaks.preValleyTime;
    whiskValleyTimesThisSwp = whiskValleys.valleyTime;

    % Skip this sweep if no whisks or basal respiration peaks
    if isempty(basalRespPeaks) || isempty(whiskPeaks)
        continue;
    end

    % Count the number of basal respiration peaks for this sweep
    nBasalRespThisSwp = height(basalRespPeaks);

    % 1. Find inspiratory whisks within each basal respiration peak
    whiskPeakNumber = nan(nBasalRespThisSwp, 1);
    for iResp = 1:nBasalRespThisSwp
        respPeakTime = basalRespPeaks.basalRespPeakTimes(iResp);
        respPreValley = basalRespPeaks.basalRespPreValleyTimes(iResp);
        respPostValley = basalRespPeaks.basalRespPostValleyTimes(iResp);

        % Handle NaN in postValleyTime
        if isnan(respPostValley)
            respPostValley = tVecsThisFile(end, iSwp);
        end

        % Moore et al 2013 definition: Find the first whisk peak with pre-valley 
        %   after breathOnsetLatency before respiration pre-valley and
        %   before respiration peak
        firstPreValleyNumber = ...
            find(whiskPreValleyTimesThisSwp >= respPreValley - breathOnsetLatency & ...
                    whiskPreValleyTimesThisSwp < respPeakTime, 1 , 'first');
        
        % Find whisk peaks within the respiration valley-to-valley window
        possibleWhiskPeakNumbers = ...
            find(whiskPeakTimesThisSwp >= respPreValley & ...
                    whiskPeakTimesThisSwp < respPostValley);
    
        % Find the whisk peak closest to the basal respiration peak
        [~, closestIdxInPossible] = min(abs(whiskPeakTimesThisSwp(possibleWhiskPeakNumbers) - respPeakTime));
        closestWhiskPeakNumber = possibleWhiskPeakNumbers(closestIdxInPossible);

        % If the closest whisk peak is different from the one with 
        %   the first prevalley, print message
        if firstPreValleyNumber ~= closestWhiskPeakNumber
            fprintf(fileID, ['First whisk protraction (peak number %d) within basal respiration is different ', ...
                     'from closest peak (peak number %d) for basal respiration cycle #%d in sweep %d of file %d!!\n\n'], ...
                     firstPreValleyNumber, closestWhiskPeakNumber, iResp, iSwp, fileNumber);
        end

        % Skip this respiration if no whisks within the basal respiration peak
        if isempty(closestWhiskPeakNumber)
%        if isempty(firstPreValleyNumber)
            continue;
        end

        % Set first whisk peak as 'inspiratory whisk'
        whiskPeakNumber(iResp) = closestWhiskPeakNumber;
%        whiskPeakNumber(iResp) = firstPreValleyNumber;
    end

    % Add inspiratory whisk peak number to basal respiratory peak table for this sweep
    basalRespPeaks.whiskPeakNumber = whiskPeakNumber;

    % Check if inspiratory whisk found for each row
    hasInspiratoryWhisk = ~isnan(whiskPeakNumber);

    % Skip this sweep if no inspiratory whisk found
    if ~any(hasInspiratoryWhisk)
        continue;
    end

    % Get the valid whisk peak numbers
    whiskPeakNumberFound = whiskPeakNumber(hasInspiratoryWhisk);

    % Remove rows with no inspiratory whisk found)
    inspWhiskTableTemp = basalRespPeaks(hasInspiratoryWhisk, :);

    % Add associated whisk peak information
    inspWhiskTable = horzcat(inspWhiskTableTemp, whiskPeaks(whiskPeakNumberFound, :));

    % Count the number of inspiratory whisks
    nInspWhisk = height(inspWhiskTable);

    % Get all inspiratory whisk times this sweep
    inspWhiskTimesThisSwp = inspWhiskTable.peakTime;

    % 2. Define basal respiration cycles
    iCycle = 0;             % Counter for basal respiration cycles for this sweep
    for iInsp = 1:nInspWhisk       
        % Get the current inspiratory whisk info
        currentInspWhisk = inspWhiskTable(iInsp, :);
        inspWhiskPeakNum = currentInspWhisk.whiskPeakNumber;
        inspWhiskPeakTime = currentInspWhisk.peakTime;
        inspWhiskPreValleyTime = currentInspWhisk.preValleyTime;
        basalRespPeakNumber = currentInspWhisk.basalRespPeakNumber;
        basalRespPeakTime = currentInspWhisk.basalRespPeakTimes;
        basalRespPreValleyTime = currentInspWhisk.basalRespPreValleyTimes;
        basalRespPostValleyTime = currentInspWhisk.basalRespPostValleyTimes;
        basalRespSucceedingIPI = currentInspWhisk.basalRespSucceedingIPIs;

        % Find the next inspiratory whisk time
        if iInsp < nInspWhisk
            nextInspWhiskTime = inspWhiskTimesThisSwp(iInsp + 1);
        else
            % No next inspiratory whisk, set at end of sweep
            nextInspWhiskTime = tVecsThisFile(end, iSwp);
        end
              
        % The search end time is the earliest of the next insp whisk or 
        %   the end of the current respiration (the succeeding valley)
        searchEndTime = min(nextInspWhiskTime, basalRespPostValleyTime);
        
        % Find intervening whisks after the current inspiratory whisk
        %   and before the search end time
        interWhiskPeakNumbers = ...
            find(whiskPeakTimesThisSwp > inspWhiskPeakTime & ...
                whiskPeakTimesThisSwp < searchEndTime);
        
        % Define a basal respiration cycle to be analyzed as 
        %   one in which there are enough intervening whisk for a cycle
        %   and that all whisk peaks have valid amplitudes
        if numel(interWhiskPeakNumbers) >= minWhisksBasalRespToAnalyze - 1
            % Collect all whisks for this cycle (inspiratory + all intervening)
            whiskIndicesForCycle = [inspWhiskPeakNum; interWhiskPeakNumbers];
            whiskPeaksToAnalyze = whiskPeaks(whiskIndicesForCycle, :);

            % Skip this inspiratory whisk if some whisk peak
            %   does not have a valid amplitude
            if any(isnan(whiskPeaksToAnalyze.amplitude))
                continue;
            end

            % Increment cycle number
            iCycle = iCycle + 1;

            % Cycle start time is the preceding valley of the inspiratory whisk
            %   or if not present, the preceding valley of the basal respiration
            if ~isnan(inspWhiskPreValleyTime)
                basalRespCycleStartTime = inspWhiskPreValleyTime;
            else
                basalRespCycleStartTime = basalRespPreValleyTime;
            end

            % Cycle end time is the succeeding valley of the last intervening whisk
            %   or if not present, the search end time
            lastWhiskInCycle = whiskPeaksToAnalyze(end, :);
            lastWhiskPostValleyTime = lastWhiskInCycle.postValleyTime;
            if ~isnan(lastWhiskPostValleyTime)
                basalRespCycleEndTime = lastWhiskPostValleyTime;
            else
                basalRespCycleEndTime = searchEndTime;
            end
            
            % Find valleys within the cycle
            isValleyInCycle = whiskValleyTimesThisSwp >= basalRespCycleStartTime & ...
                              whiskValleyTimesThisSwp <= basalRespCycleEndTime;
            whiskValleysToAnalyze = whiskValleys(isValleyInCycle, :);

            % Compute log decrements
            peakAmplitudes = whiskPeaksToAnalyze.amplitude;
            if numel(peakAmplitudes) > 1
                logDecrements = log(peakAmplitudes(2:end) ./ peakAmplitudes(1:end-1));
            else
                logDecrements = [];
            end
                       
            % Create a one-row table for this cycle
            newCycle = table(iSwp, iCycle, basalRespPeakNumber, ...
                basalRespCycleStartTime, basalRespCycleEndTime, ...
                basalRespPeakTime, basalRespPreValleyTime, basalRespPostValleyTime, basalRespSucceedingIPI, ...
                {whiskPeaksToAnalyze.peakTime}, {whiskPeaksToAnalyze.peakValue}, ...
                {peakAmplitudes}, {whiskPeaksToAnalyze.preValleyTime}, {whiskPeaksToAnalyze.postValleyTime}, ...
                {whiskValleysToAnalyze.valleyTime}, {whiskValleysToAnalyze.valleyValue}, {logDecrements}, ...
                'VariableNames', {'sweepNumber', 'basalRespCycleNumber', 'basalRespPeakNumber', ...
                      'basalRespCycleStartTime', 'basalRespCycleEndTime', ...
                      'basalRespPeakTime', 'basalRespPreValleyTime', 'basalRespPostValleyTime', 'basalRespSucceedingIPI', ...
                      'whiskPeakTimes', 'whiskPeakValues', 'whiskPeakAmplitudes', ...
                      'whiskPreValleyTimes', 'whiskPostValleyTimes', ...
                      'whiskValleyTimes', 'whiskValleyValues', 'whiskLogDecrements'});

            allCyclesCell{end+1} = newCycle;
        end
    end
end

% Vertically concatenate all found cycles into a single table
if ~isempty(allCyclesCell)
    basalRespCycleTable = vertcat(allCyclesCell{:});
else
    basalRespCycleTable = emptyBasalTable;
end

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
                identify_files_for_averaging(trialNames, excludeStringsFromAverage)
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

function combinedTable = combine_across_files(tablesAll, ...
                                        fileNumsToAverage, trialNames)
%% Combines sniff start window tables from selected files into one large table.

% Initialize a cell array to hold tables that will be concatenated
tablesToCombine = {};

% Loop through each file number selected for averaging
for i = 1:numel(fileNumsToAverage)
    fileNum = fileNumsToAverage(i);
    
    % Get the sniff start window table for the current file
    currentTable = tablesAll{fileNum};
    
    % Proceed only if the table is not empty
    if ~isempty(currentTable)
        % Add new columns for the file number and trial name
        fileNumColumn = repelem(fileNum, height(currentTable), 1);
        trialNameColumn = repelem(trialNames(fileNum), height(currentTable), 1);
        
        currentTable.fileNumber = fileNumColumn;
        currentTable.trialName = trialNameColumn;
        
        % Add the modified table to our list
        tablesToCombine{end+1} = currentTable;
    end
end

% Vertically concatenate all tables in the list into a single table
if ~isempty(tablesToCombine)
    combinedTable = vertcat(tablesToCombine{:});
else
    % Return an empty table if no data was found
    combinedTable = table();
end

% Reorder columns
if ~isempty(combinedTable)
    combinedTable = movevars(combinedTable, 'fileNumber', 'Before', 1);
    combinedTable = movevars(combinedTable, 'trialName', 'Before', 1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results, handles] = average_log_decrement_jitter(T, maxDecrementsToAnalyze, pathOutDir, figTypes, figTitle, figName)
%% Averages and plots whisk logarithmic decrements as a grouped jitter plot,
%  with mean, 95% CI, and statistical test results overlaid.

% Initialize output
results = struct;
handles.fig = gobjects;

% If there's no data to plot, exit early
if isempty(T) || ~ismember('whiskLogDecrements', T.Properties.VariableNames)
    fprintf('No whisk log decrements to average. Skipping jitter plot!\n');
    return;
end

% Extract the necessary columns from the combined table
whiskLogDecrementsCell = T.whiskLogDecrements;
fileNumbers = T.fileNumber;

% Check if there is anything to average after unpacking
if isempty(whiskLogDecrementsCell) || all(cellfun(@isempty, whiskLogDecrementsCell))
    fprintf('No log decrement data found to plot.\n');
    return;
end

% Convert cell array to a matrix, with each column being a decrement order
%   and each row being an analysis window
allLogDecrementsMatrix = transpose(force_matrix(whiskLogDecrementsCell, 'CombineMethod', 'leftAdjustPad'));

% Count the number of analysis windows
nAnalysisWindows = size(allLogDecrementsMatrix, 1);

% Count the number of log decrements
nWhiskLogDecrements = size(allLogDecrementsMatrix, 2);

% Restrict to maxDecrementsToAnalyze
if nWhiskLogDecrements > maxDecrementsToAnalyze; 
    nDecrementsToAnalyze = maxDecrementsToAnalyze;
    allLogDecrementsMatrix = allLogDecrementsMatrix(:, 1:nDecrementsToAnalyze);
else
    nDecrementsToAnalyze = nWhiskLogDecrements;
end

% Create matching decrement orders and file numbers
allDecrementOrdersMatrix = repmat((1:nDecrementsToAnalyze), nAnalysisWindows, 1);
allFileNumbersMatrix = repmat(fileNumbers, 1, nDecrementsToAnalyze);

%% Compute statistics
% Calculate mean and 95% confidence intervals
[meanWhiskLogDecrements, stderrWhiskLogDecrements, ...
    lower95WhiskLogDecrements, upper95WhiskLogDecrements] = ...
    compute_stats_for_cellnumeric(allLogDecrementsMatrix');

% Perform significance tests for each decrement
stats = vecfun(@test_difference, allLogDecrementsMatrix);

% Extract results
pValues = extract_fields(stats, 'pValue');
testFunctions = extract_fields(stats, 'testFunction');
symbols = extract_fields(stats, 'symbol');

% Get the geometric mean of amplitude ratios
avgWhiskAmpRatios = exp(meanWhiskLogDecrements);

%% Save results
results.allLogDecrements = allLogDecrementsMatrix;
results.allDecrementOrders = allDecrementOrdersMatrix;
results.allFileNumbers = allFileNumbersMatrix;
results.nAnalysisWindows = nAnalysisWindows;
results.nWhiskLogDecrements = nWhiskLogDecrements;
results.nDecrementsToAnalyze = nDecrementsToAnalyze;
results.meanWhiskLogDecrements = meanWhiskLogDecrements;
results.stderrWhiskLogDecrements = stderrWhiskLogDecrements;
results.lower95WhiskLogDecrements = lower95WhiskLogDecrements;
results.upper95WhiskLogDecrements = upper95WhiskLogDecrements;
results.pValues = pValues;
results.testFunctions = testFunctions;
results.symbols = symbols;
results.avgWhiskAmpRatios = avgWhiskAmpRatios;

%% Prepare data for the jitter plot
allLogDecrementsVec = allLogDecrementsMatrix(:);
allDecrementOrdersVec = allDecrementOrdersMatrix(:);
allFileNumbersVec = allFileNumbersMatrix(:);
xTickLabels = arrayfun(@(x) sprintf('ln(A%d/A%d)', x+1, x), 1:nDecrementsToAnalyze, 'UniformOutput', false);

%% Plot the jitter plot and overlay stats
% Set up figure properties
fig = set_figure_properties('AlwaysNew', true, 'ClearFigure', true);
figPath = fullfile(pathOutDir, figName);

% 1. Generate the base jitter plot, ensuring it doesn't plot its own stats
handles = plot_grouped_jitter(allLogDecrementsVec, allFileNumbersVec, allDecrementOrdersVec, ...
    'XTickLabels', xTickLabels, ...
    'YLabel', 'Log Decrement (ln(A_{n+1}/A_{n}))', ...
    'LegendLocation', 'suppress', ...
    'PlotMeanValues', false, 'PlotErrorBars', false, ...
    'RunTTest', false, 'RunRankTest', false);

% 2. Hold the plot to overlay new elements
hold on;

% 3. Plot the mean and 95% CI error bars with a distinct style
xValues = (1:nDecrementsToAnalyze)';
errLower = meanWhiskLogDecrements - lower95WhiskLogDecrements;
errUpper = upper95WhiskLogDecrements - meanWhiskLogDecrements;

handles.errorbar = ...
    errorbar(xValues, meanWhiskLogDecrements, errLower, errUpper, '_', ...
         'Color', 'k', 'LineWidth', 2.5, 'CapSize', 20, 'Marker', 'none');
handles.mean = ...
    plot(xValues, meanWhiskLogDecrements, 'o', 'MarkerEdgeColor', 'k', ...
     'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 1.5);

% 4. Add a horizontal line at y=0 for the null hypothesis
handles.null = yline(0, '--k', 'LineWidth', 1);

% 5. Annotate plot with symbol and p-values near the top of the plot area
handles.testResult = ...
    plot_test_result(pValues, 'TestFunction', testFunctions, 'Symbol', symbols, ...        
                 'XLocText', xValues, 'YLocTextRel', 0.90, 'YLocStarRel', 0.95);

% 6. Add geometric mean of amplitude ratios below the p-values
% Get y-axis limits to position the text
yLimits = ylim;
% Position text at 85% of the y-axis height (below p-values at 90%)
yPosText = yLimits(1) + 0.85 * diff(yLimits); 
% Create the text labels for each point
ratioLabels = arrayfun(@(x) sprintf('Ratio: %.2f', x), avgWhiskAmpRatios, 'UniformOutput', false);
% Add the text to the plot
handles.ratioText = text(xValues, repmat(yPosText, size(xValues)), ratioLabels, ...
                          'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% 7. Finalize plot
title(figTitle);
grid on;
hold off;

% 8. Save the figure
save_all_figtypes(fig, figPath, figTypes);

% 9. Return in handles
handles.fig = fig;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results, handles] = ...
    correlate_successive_whiskamp(T, nCorrToAnalyze, pathOutDir, ...
                            figTypes, figTitle, figName, ...
                            markerTypeScatter, markerSizeScatter, markerLineWidthScatter)
%% This function correlates the amplitude of successive whisks 
% (e.g., A1 vs A2, A2 vs A3, etc.) and plots these correlations in 
% separate subplots, with points colored by fileNumber.

%% Preparation
% Initialize output
results = struct;
handles.fig = gobjects;

% If there's no data to plot, exit early
if isempty(T) || ~ismember('whiskPeakAmplitudes', T.Properties.VariableNames)
    fprintf('No whisk peak amplitudes to correlate. Skipping scatter plot!\n');
    return;
end

% Extract the necessary columns from the combined table
whiskPeakAmplitudesCell = T.whiskPeakAmplitudes;
fileNumbers = T.fileNumber;

% Check if there is anything to average after unpacking
if isempty(whiskPeakAmplitudesCell) || all(cellfun(@isempty, whiskPeakAmplitudesCell))
    fprintf('No whisk peak amplitudes data found to plot.\n');
    return;
end

% Convert cell array to a matrix, with each column being a peak order
%   and each row being an analysis window
whiskAmplitudesMatrix = transpose(force_matrix(whiskPeakAmplitudesCell, 'CombineMethod', 'leftAdjustPad'));

% Count the number of analysis windows
nAnalysisWindows = size(whiskAmplitudesMatrix, 1);

% Count the number of whisk peak orders per window
nWhiskPeakOrders = size(whiskAmplitudesMatrix, 2);

% Determine the number of correlations to compute
if nCorrToAnalyze > nWhiskPeakOrders - 1
    disp('Not enough whisk data to correlate successive amplitudes.');
    return;
end

%% Plot scatter plots and compute correlations
% Create a figure with as many subplots as correlations
[fig, ax] = create_subplots(nCorrToAnalyze, 'AlwaysNew', true, ...
                        'ClearFigure', true, 'FigExpansion', [1, 1]);
figPath = fullfile(pathOutDir, figName);

% Plot and compute each correlation and plot for each successive pair
isSignificant = nan(nCorrToAnalyze, 1);
corrCoeffs = nan(nCorrToAnalyze, 1);
pValues = nan(nCorrToAnalyze, 1);
for iCorr = 1:nCorrToAnalyze
    % Select the current subplot axes
    axes(ax(iCorr));

    % Extract amplitude data for whisk iCorr and whisk iCorr+1
    ampCurrent = whiskAmplitudesMatrix(:, iCorr);
    ampNext = whiskAmplitudesMatrix(:, iCorr + 1);

    % Set labels
    xLabel = ['Amplitude of Whisk Peak #', num2str(iCorr), ' (deg)'];
    yLabel = ['Amplitude of Whisk Peak #', num2str(iCorr + 1), ' (deg)'];

    % Plot the amplitudes against each other, colored by fileNumber
    handles = plot_grouped_scatter(ampCurrent, ampNext, fileNumbers, ...
                        'PlotEllipse', false, 'MarkerType', markerTypeScatter, ...
                        'MarkerSize', markerSizeScatter, ...
                        'MarkerLineWidth', markerLineWidthScatter, ...
                        'XLabel', xLabel, 'YLabel', yLabel, ...
                        'FigTitle', 'suppress', 'LegendLocation', 'suppress');

    % Compute and plot the correlation coefficient, ignoring NaN values
    [textObjects, isSignificant(iCorr), corrCoeffs(iCorr), pValues(iCorr)] = ...
        plot_correlation_coefficient('XData', ampCurrent, 'YData', ampNext);

    % Improve aesthetics
    grid on;
    hold off;
end

% Finalize plot
[supAx, ax] = resize_subplots_for_labels('FigTitle', figTitle);

% Save the figure
save_all_figtypes(fig, figPath, figTypes);

% 9. Return in handles
handles.fig = fig;
handles.ax = ax;
handles.supAx = supAx;
handles.textObjects = textObjects;

%% Save results
results.whiskAmplitudesMatrix = whiskAmplitudesMatrix;
results.fileNumbers = fileNumbers;
results.nAnalysisWindows = nAnalysisWindows;
results.nWhiskPeakOrders = nWhiskPeakOrders;
results.nCorrToAnalyze = nCorrToAnalyze;
results.corrCoeffs = corrCoeffs;
results.pValues = pValues;
results.isSignificant = isSignificant;

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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%