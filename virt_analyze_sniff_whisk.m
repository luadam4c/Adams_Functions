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
%       cd/argfun.m
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
%       cd/plot_test_result.m
%       cd/plot_traces.m
%       cd/plot_vertical_line.m
%       cd/plot_vertical_shade.m
%       cd/print_cellstr.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/test_difference.m
%       cd/vecfun.m
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
% TODO: Analyze basal respiration windows
% TODO: Fix sniff vec filter cutoff to [1, 15] Hz, filter order to 3
% TODO: Fix whisk vec filter cutoff to [3, 25] Hz, filter order to 3

%% Hard-coded parameters
% Directory and file naming conventions
nameDataDir = 'data_sniff_whisk';       % Name of the directory containing the data files
nameOutDir = 'output_sniff_whisk';      % Name of the output directory
fileNameMetaData = 'sniff_whisk_metadata.csv';      % File name of the meta data file
extDataFile = 'mat';                    % Extension for the data files
prefixDataFile = 'swdata_';             % Prefix for the data files
suffixDataFile = '_angle';              % Suffix for the data files
suffixParamFile = '_parameters';        % Suffix for the parameter files

% Detection parameters
ammoniaPuffString = 'ammpuff';
airPuffString = 'airpuff';
minPulseAmplitude = 2;                  % Minimum pulse amplitude for piezo trace

% Analysis parameters
%   Note: Keep this consistent with virt_moore.m
amplitudeDefinition = 'peak-to-avgvalley';
fundFreqRange = [0.5, 20];      % range of possible fundamental frequencies to be detected
fCutoffRelToFund = [0.1, 10];   % ratio of Butterworth bandpass filter cutoff for effector trace
                                % relative to the fundamental frequency
filterOrder = 2;            % filter order for effector trace
promThresholdPerc = 10;     % Percentage of amplitude range for minimum peak prominence
minPeakDistanceMs = 30;     % Minimum peak distance in ms
sniffIpiThresholdMs = 250;  % Inter-peak interval threshold for sniffing in ms
nWhisksToAnalyze = 5;       % Number of whisks at the start of a sniff period to be analyzed

% Hard-coded strings in file names to exclude from averaging
excludeStringsFromAverage = {'ammpuff', 'airpuff', 'baseline', 'eth'};

% Plotting parameters
%fileNumsToPlot = 38;                    % The file number(s) to plot (max 38)
%fileNumsToPlot = 4;                    % The file number(s) to plot (max 38)
fileNumsToPlot = (1:38)';              % The file number(s) to plot (max 38)
toSpeedUp = false;                      % Whether to use parpool and hide figures
%toSpeedUp = true;                      % Whether to use parpool and hide figures
whiskAngleLimits = [-75, 75];           % Whisk angle limits to be plotted
piezoLimits = [-1, 10];
colorWhisk = [0, 0, 1];                 % Color for whisk trace (Blue)
colorSniff = [1, 0, 0];                 % Color for sniff trace (Red)
colorStim = [0, 0, 0];                  % Color for stim trace (Black)
colorAmmoniaPuff = [0.6, 0.8, 0.2];     % Color for Ammonia Puff (Yellow Green)
colorAirPuff = 0.5 * [1, 1, 1];         % Color for Air Puff (Gray)
colorAnalysisWin = [0.5, 1.0, 0.8];     % Color for sniff start windows (Aquamarine)
faceAlphaAnalysisWin = 0.8;             % Transparencies for sniff start windows
markerSniffPeaksValleys = 'o';          % Marker for sniff peaks and valleys (circle)
colorSniffPeaksValleys = [1, 0.7, 0];   % Color for sniff peaks and valleys (Orange)
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
whiskAngleLabel = 'Whisk Angle (degrees)';
sniffToWhiskRangeRatio = 0.8;           % Sniff amplitude to whisk amplitude ratios for plotting
timeLabel = 'Time (s)';
whiskLabel = 'Whisking';
sniffLabel = 'Breathing';
legendLocation1 = 'northeast';
legendLocation2 = 'suppress';
legendLocation3 = 'suppress';
figTitlePrefix1 = 'Whisking (blue) and breathing (red) data';
figTitlePrefix2 = 'Stimulation (Puff) data';
figTitlePrefix3 = 'Sniff Start Windows';
figTitle4 = 'Whisk Logarithmic Decrements';
figPrefix1 = 'sniff_whisk_all_traces_';
figPrefix2 = 'sniff_whisk_all_stims_';
figPrefix3 = 'sniff_whisk_sniffstart_windows_';
figName4 = 'whisk_log_decrements_jitter';
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
rangeThermLimitsEachFile = (maxThermValueEachFile - minThermValueEachFile) / sniffToWhiskRangeRatio;
lowerThermLimitEachFile = meanThermLimitsEachFile - rangeThermLimitsEachFile / 2;

% Create sniff/breath vectors scaled to the order of whisk angle vectors
rangeWhiskAngleLimits = diff(whiskAngleLimits);
sniffVecs = cellfun(@(a, b, c) whiskAngleLimits(1) + (a - b) * ...
                    rangeWhiskAngleLimits / c, thermVecs, ...
                    num2cell(lowerThermLimitEachFile), num2cell(rangeThermLimitsEachFile), ...
                    'UniformOutput', false);

%% Detect the first pulse start and end times
pulseParams = cellfun(@(x, y) parse_pulse(x, 'TimeVecs', y, ...
        'MinPulseAmplitude', minPulseAmplitude), pulseVecs, tVecs, 'UniformOutput', false);
pulseStartTimes = cellfun(@(x) x.timeAfterStartMs, pulseParams, 'UniformOutput', false);
pulseEndTimes = cellfun(@(x) x.timeBeforeEndMs, pulseParams, 'UniformOutput', false);

%% Detect sniff transitions
fprintf('Detecting sniff transitions...\n');

% Use array_fun to apply the function to each file's data
[sniffPeakTablesAll, sniffValleyTablesAll, sniffStartTimesAll, sniffEndTimesAll, sniffFreqsFundamentalAll] = ...
    array_fun(@(x, y, z) parse_sniff_vecs(x, y, z, ...
                    amplitudeDefinition, fundFreqRange, fCutoffRelToFund, filterOrder, ...
                    promThresholdPerc, minPeakDistanceMs, ...
                    sniffIpiThresholdMs), ...
                sniffVecs, tVecs, num2cell(nSweepsEachFile), 'UniformOutput', false);


fprintf('Finished detecting sniff transitions.\n\n');

%% Detect whisk peaks and valleys
fprintf('Detecting whisk peaks and valleys...\n');
[whiskPeakTablesAll, whiskValleyTablesAll, whiskFreqsFundamentalAll, sniffStartWinTablesAll] = ...
    array_fun(@(a, b, c, d, e) parse_whisk_vecs(a, b, c, d, e, ...
                    amplitudeDefinition, fundFreqRange, fCutoffRelToFund, filterOrder, ...
                    promThresholdPerc, minPeakDistanceMs, nWhisksToAnalyze), ...
                whiskVecs, tVecs, num2cell(nSweepsEachFile), ...
                sniffStartTimesAll, sniffEndTimesAll, 'UniformOutput', false, ...
                'UseParpool', false);
fprintf('Finished detecting whisk peaks and valleys.\n\n');

%% Augment Sniff Start Windows with Sniff Data
fprintf('Augmenting sniff start windows with sniff data...\n');
sniffStartWinTablesAll = ...
    array_fun(@(a, b, c) augment_sniffstart_windows(a, b, c), ...
                sniffStartWinTablesAll, sniffPeakTablesAll, sniffValleyTablesAll, ...
                'UniformOutput', false);
fprintf('Finished augmenting sniff start windows.\n\n');

%% Calculate whisk logarithmic decrement statistics per file
fprintf('Calculating statistics for whisk logarithmic decrements for each file ...\n');

% Calculate the number of sniff start windows per file
nAnalysisWindowsPerFile = cellfun(@height, sniffStartWinTablesAll);

% Calculate the average whisk logarithmic decrements per file
meanWhiskLogDecrementsPerFile = ...
    cellfun(@(x) compute_combined_trace(x.whiskLogDecrements, 'mean'), ...
            sniffStartWinTablesAll, 'UniformOutput', false);

% Calculate the standard error of whisk logarithmic decrements per file
stderrWhiskLogDecrementsPerFile = ...
    cellfun(@(x) compute_combined_trace(x.whiskLogDecrements, 'stderr'), ...
            sniffStartWinTablesAll, 'UniformOutput', false);

% Calculate the lower bound of 95% confidence interval of whisk logarithmic decrements per file
lower95WhiskLogDecrementsPerFile = ...
    cellfun(@(x) compute_combined_trace(x.whiskLogDecrements, 'lower95'), ...
            sniffStartWinTablesAll, 'UniformOutput', false);

% Calculate the upper bound of 95% confidence interval of whisk logarithmic decrements per file
upper95WhiskLogDecrementsPerFile = ...
    cellfun(@(x) compute_combined_trace(x.whiskLogDecrements, 'upper95'), ...
            sniffStartWinTablesAll, 'UniformOutput', false);

fprintf('Finished calculating ratios.\n\n');

%% Plot aggregate data
% Identify files to use for averaging
[~, fileNumsToAverage] = identify_files_for_averaging(trialNames, excludeStringsFromAverage);

% Combine sniff start windows from selected files
sniffStartWinTableToAverage = combine_sniffstart_windows(sniffStartWinTablesAll, ...
                                    fileNumsToAverage, trialNames);

% Plot the whisk logarithmic decrements as a grouped jitter plot
average_log_decrement_jitter(sniffStartWinTableToAverage, pathOutDir, figTypes, figTitle4, figName4);

%% Plot all data
[handles1, handles2, handles3] = ...
    array_fun(@(a) plot_one_file (a, trialNames, tVecs, whiskVecs, sniffVecs, ...
        pulseVecs, pulseStartTimes, pulseEndTimes, ...
        sniffPeakTablesAll, sniffValleyTablesAll, sniffStartTimesAll, sniffEndTimesAll, ...
        whiskPeakTablesAll, whiskValleyTablesAll, sniffStartWinTablesAll, ...
        isAmmoniaPuff, isAirPuff, whiskLabel, sniffLabel, ...
        colorWhisk, colorSniff, colorStim, colorAmmoniaPuff, colorAirPuff, ...
        colorAnalysisWin, faceAlphaAnalysisWin, ...
        markerSniffPeaksValleys, colorSniffPeaksValleys, ...
        markerWhiskPeaksValleys, colorWhiskPeaksValleys, ...
        lineStyleWhiskAmplitudes, colorWhiskAmplitudes, ...
        colorSniffStart, colorSniffEnd, ...
        lineWidthForSample, lineWidthForAnalysis, ...
        markerSizeForSample, markerSizeForAnalysis, ...
        pathOutDir, timeLabel, whiskAngleLimits, whiskAngleLabel, piezoLimits, ...
        legendLocation1, legendLocation2, legendLocation3, ...
        figTitlePrefix1, figTitlePrefix2, figTitlePrefix3, ...
        figPrefix1, figPrefix2, figPrefix3, figTypes, toSpeedUp), ...
        fileNumsToPlot, 'UseParpool', toSpeedUp);

%% Put all metadata together and save

% Load parameters from each file
paramStructs = cellfun(@load, pathParamFiles);
paramTable = struct2table(paramStructs);

% Check if any matching pairs were found
if isempty(trialNames)
    fprintf('No matching data and parameter file pairs were found in the directory.\n');
else
    % Create a table for all meta data
    metaDataFiles = table(trialNames, pathDataFiles, pathParamFiles, ...
                     'VariableNames', {'TrialName', 'DataFilePath', 'ParameterFilePath'});
    metaData = horzcat(metaDataFiles, paramTable);

    % Add detection results
    metaData.nSweeps = nSweepsEachFile;
    metaData.isAmmoniaPuff = isAmmoniaPuff;
    metaData.isAirPuff = isAirPuff;
    metaData.sniffStartTimes = sniffStartTimesAll;
    metaData.sniffEndTimes = sniffEndTimesAll;
    metaData.sniffFreqFundamental = sniffFreqsFundamentalAll;
    metaData.whiskFreqFundamental = whiskFreqsFundamentalAll;
    metaData.nAnalysisWindows = nAnalysisWindowsPerFile;
    metaData.meanWhiskLogDecrements = meanWhiskLogDecrementsPerFile;
    metaData.stderrWhiskLogDecrements = stderrWhiskLogDecrementsPerFile;
    metaData.lower95WhiskLogDecrements = lower95WhiskLogDecrementsPerFile;
    metaData.upper95WhiskLogDecrements = upper95WhiskLogDecrementsPerFile;

    % Expand the cell arrays to a single string
    metaDataToPrint = metaData;
    [metaDataToPrint.roiswitch, metaDataToPrint.records, ...
        metaDataToPrint.anesrecords, metaDataToPrint.sniffStartTimes, ...
        metaDataToPrint.sniffEndTimes, metaDataToPrint.sniffFreqFundamental, ...
        metaDataToPrint.whiskFreqFundamental, ...
        metaData.meanWhiskLogDecrements, metaData.stderrWhiskLogDecrements, ...
        metaData.lower95WhiskLogDecrements, metaData.upper95WhiskLogDecrements] = ...
        argfun(@(x) cellfun(@(a) print_cellstr(a, 'Delimiter', ' ', ...
                                'ToPrint', false, 'OmitQuotes', true, ...
                                'OmitBraces', true, 'OmitNewline', true), ...
                                x, 'UniformOutput', false), ...
            metaData.roiswitch, metaData.records, ...
            metaData.anesrecords, metaData.sniffStartTimes, ...
            metaData.sniffEndTimes, metaData.sniffFreqFundamental, ...
            metaData.whiskFreqFundamental, ...
            metaData.meanWhiskLogDecrements, metaData.stderrWhiskLogDecrements, ...
            metaData.lower95WhiskLogDecrements, metaData.upper95WhiskLogDecrements);

    % Display the resulting table in the command window
    fprintf('Generated Metadata Table:\n');
    disp(metaDataToPrint);

    % Save metadata table
    pathMetaData = fullfile(pathOutDir, fileNameMetaData);
    writetable(metaDataToPrint, pathMetaData);
    fprintf('Metadata table saved to %s\n', pathMetaData);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [handles1, handles2, handles3] = ...
    plot_one_file (sampleFileNum, trialNames, tVecs, whiskVecs, sniffVecs, ...
        pulseVecs, pulseStartTimes, pulseEndTimes, ...
        sniffPeakTablesAll, sniffValleyTablesAll, sniffStartTimesAll, sniffEndTimesAll, ...
        whiskPeakTablesAll, whiskValleyTablesAll, sniffStartWinTablesAll, ...
        isAmmoniaPuff, isAirPuff, whiskLabel, sniffLabel, ...
        colorWhisk, colorSniff, colorStim, colorAmmoniaPuff, colorAirPuff, ...
        colorAnalysisWin, faceAlphaAnalysisWin, ...
        markerSniffPeaksValleys, colorSniffPeaksValleys, ...
        markerWhiskPeaksValleys, colorWhiskPeaksValleys, ...
        lineStyleWhiskAmplitudes, colorWhiskAmplitudes, ...
        colorSniffStart, colorSniffEnd, ...
        lineWidthForSample, lineWidthForAnalysis, ...
        markerSizeForSample, markerSizeForAnalysis, ...
        pathOutDir, timeLabel, whiskAngleLimits, whiskAngleLabel, piezoLimits, ...
        legendLocation1, legendLocation2, legendLocation3, ...
        figTitlePrefix1, figTitlePrefix2, figTitlePrefix3, ...
        figPrefix1, figPrefix2, figPrefix3, figTypes, toSpeedUp)

% Extract data to plot
trialName = trialNames{sampleFileNum};
tVecsThis = tVecs{sampleFileNum};
whiskVecsThis = whiskVecs{sampleFileNum};
sniffVecsThis = sniffVecs{sampleFileNum};
pulseVecsThis = pulseVecs{sampleFileNum};
pulseStartTimesThis = pulseStartTimes{sampleFileNum};
pulseEndTimesThis = pulseEndTimes{sampleFileNum};
isAmmoniaPuffThis = isAmmoniaPuff(sampleFileNum);
isAirPuffThis = isAirPuff(sampleFileNum);
sniffPeakTablesThis = sniffPeakTablesAll{sampleFileNum};
sniffValleyTablesThis = sniffValleyTablesAll{sampleFileNum};
sniffStartTimesThis = sniffStartTimesAll{sampleFileNum};
sniffEndTimesThis = sniffEndTimesAll{sampleFileNum};
whiskPeakTablesThis = whiskPeakTablesAll{sampleFileNum};
whiskValleyTablesThis = whiskValleyTablesAll{sampleFileNum};
sniffStartWinTableForFile = sniffStartWinTablesAll{sampleFileNum};

% Count the number of records to plot
nSweeps = size(tVecsThis, 2);
recordNumbers = (1:nSweeps)';

% Create labels for legend
whiskLabels = create_labels_from_numbers(recordNumbers, 'Prefix', whiskLabel);
sniffLabels = create_labels_from_numbers(recordNumbers, 'Prefix', sniffLabel);

% Decide on puff window colors and titles
if isAmmoniaPuffThis
    figTitlePrefix1 = replace(figTitlePrefix1, 'data', 'data with Ammonia Puffs (green)');
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

% Create paths for saving
figName1 = [figPrefix1, 'File', num2str(sampleFileNum), '_', trialName];
figName2 = [figPrefix2, 'File', num2str(sampleFileNum), '_', trialName];
figName3 = [figPrefix3, 'File', num2str(sampleFileNum), '_', trialName];
figPath1 = fullfile(pathOutDir, figName1);
figPath2 = fullfile(pathOutDir, figName2);
figPath3 = fullfile(pathOutDir, figName3);

% Create sample traces plot
fprintf('Plotting sample traces for file %d ...\n', sampleFileNum);
if toSpeedUp
    set_figure_properties('AlwaysNew', true, 'Visible', 'off');
else
    set_figure_properties('FigNumber', 1, 'AlwaysNew', false, 'ClearFigure', true);
end
if ~isempty(sniffVecsThis)
    handles1 = plot_traces(tVecsThis, whiskVecsThis, 'ColorMap', colorWhisk, ...
                'DataToCompare', sniffVecsThis, 'ColorMapToCompare', colorSniff, ...
                'XBoundaries', [pulseStartTimesThis, pulseEndTimesThis], ...
                'XBoundaryColor', puffWindowColor, ...
                'TraceLabels', whiskLabels, 'TraceLabelsToCompare', sniffLabels, ...
                'PlotMode', 'parallel', 'FigTitle', figTitle1, ...
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
            plot_vertical_shade(winBoundaries', 'Color', colorAnalysisWin, 'FaceAlpha', faceAlphaAnalysisWin);
        end

        % Extract the peak and valley tables for this sweep
        sniffPeakTable = sniffPeakTablesThis{iSwp};
        sniffValleyTable = sniffValleyTablesThis{iSwp};
        whiskPeakTable = whiskPeakTablesThis{iSwp};
        whiskValleyTable = whiskValleyTablesThis{iSwp};
        
        % Plot peaks if they exist
        if ~isempty(sniffPeakTable) && ismember('peakIndex', sniffPeakTable.Properties.VariableNames)
            plot(sniffPeakTable.peakTime, sniffPeakTable.peakValue, ...
                markerSniffPeaksValleys, 'Color', colorSniffPeaksValleys, ...
                'MarkerSize', markerSizeForSample, 'LineWidth', lineWidthForSample);
        end
        if ~isempty(whiskPeakTable) && ismember('peakIndex', whiskPeakTable.Properties.VariableNames)
            plot(whiskPeakTable.peakTime, whiskPeakTable.peakValue, ...
                markerWhiskPeaksValleys, 'Color', colorWhiskPeaksValleys, ...
                'MarkerSize', markerSizeForSample, 'LineWidth', lineWidthForSample);
        end
     
        % Plot valleys if they exist
        if ~isempty(sniffValleyTable) && ismember('valleyIndex', sniffValleyTable.Properties.VariableNames)
            plot(sniffValleyTable.valleyTime, sniffValleyTable.valleyValue, ...
                markerSniffPeaksValleys, 'Color', colorSniffPeaksValleys, ...
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
    recNums = sniffStartWinTableForFile.sweepNumber;
    winStarts = sniffStartWinTableForFile.sniffStartWinStartTime;
    winEnds = sniffStartWinTableForFile.sniffStartWinEndTime;
    sniffStarts = sniffStartWinTableForFile.sniffStartTime;
    sniffPeakTimesForWin = sniffStartWinTableForFile.sniffPeakTimes;
    sniffPeakValuesForWin = sniffStartWinTableForFile.sniffPeakValues;
    sniffValleyTimesForWin = sniffStartWinTableForFile.sniffValleyTimes;
    sniffValleyValuesForWin = sniffStartWinTableForFile.sniffValleyValues;
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
    sniffVecsForWin = cell(nWindows, 1);
    for iWin = 1:nWindows
        % Get info for the current window
        swpNum = recNums(iWin);
        winStart = winStarts(iWin);
        winEnd = winEnds(iWin);
        
        % Get the full traces for the corresponding sweep
        tFull = tVecsThis(:, swpNum);
        whiskFull = whiskVecsThis(:, swpNum);
        sniffFull = sniffVecsThis(:, swpNum);
        
        % Find indices within the current sniff start window
        indicesInWin = find(tFull >= winStart & tFull <= winEnd);
        
        % Use extract_subvectors to get the data segments
        tVecsForWin{iWin} = extract_subvectors(tFull, 'Indices', indicesInWin);
        whiskVecsForWin{iWin} = extract_subvectors(whiskFull, 'Indices', indicesInWin);
        sniffVecsForWin{iWin} = extract_subvectors(sniffFull, 'Indices', indicesInWin);
    end

    % Set figure properties for the new plot
    if toSpeedUp
        set_figure_properties('AlwaysNew', true, 'Visible', 'off');
    else
        set_figure_properties('FigNumber', 3, 'AlwaysNew', false, 'ClearFigure', true);
    end
    
    % Plot all extracted windows in parallel subplots
    handles3 = plot_traces(tVecsForWin, whiskVecsForWin, 'ColorMap', colorWhisk, ...
        'DataToCompare', sniffVecsForWin, 'ColorMapToCompare', colorSniff, ...
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
        sniffPeakTimes = sniffPeakTimesForWin{iWin};
        sniffPeakValues = sniffPeakValuesForWin{iWin};
        sniffValleyTimes = sniffValleyTimesForWin{iWin};
        sniffValleyValues = sniffValleyValuesForWin{iWin};
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
        
        % Plot sniff peaks used for analysis
        if ~isempty(sniffPeakTimes)
            plot(sniffPeakTimes, sniffPeakValues, ...
                markerSniffPeaksValleys, 'Color', colorSniffPeaksValleys, ...
                'LineWidth', lineWidthForAnalysis, 'MarkerSize', markerSizeForAnalysis);
        end

        % Plot sniff valleys used for analysis
        if ~isempty(sniffValleyTimes)
            plot(sniffValleyTimes, sniffValleyValues, ...
                markerSniffPeaksValleys, 'Color', colorSniffPeaksValleys, ...
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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sniffPeakTables, valleyTables, sniffStartTimes, sniffEndTimes, sniffFreqFundamental] = ...
    parse_sniff_vecs(sniffVecsThisFile, tVecsThisFile, nSweeps, ...
                    amplitudeDefinition, fundFreqRange, fCutoffRelToFund, filterOrder, ...
                    promThresholdPerc, minPeakDistanceMs, ...
                    sniffIpiThresholdMs)
%% Parses Sniff vectors for each file

% Hard-coded parameters for this function's logic
minSniffPeaks = 3;      % Minimum number of peaks to detect transitions

% Convert to seconds
sniffIpiThresholdSec = sniffIpiThresholdMs / 1000;

% Initialize cell arrays to store results for each sweep in this file.
sniffPeakTables = cell(nSweeps, 1);
valleyTables = cell(nSweeps, 1);
sniffStartTimes = cell(nSweeps, 1);
sniffEndTimes = cell(nSweeps, 1);
sniffFreqFundamental = cell(nSweeps, 1);

% If sniff vector is empty, return
if isempty(sniffVecsThisFile)
    return
end

% Loop through each sweep within the current file.
for iSwp = 1:nSweeps
    % Extract the sniff vector for the current sweep.
    sniffVec = sniffVecsThisFile(:, iSwp);
    % Extract the time vector for the current sweep.
    timeVec = tVecsThisFile(:, iSwp);

    % Find all peak times and preceding valley times for the current sniff vector.
    [sniffPeakTable, sniffValleyTable, otherResults] = ...
        parse_oscillation(sniffVec, ...
                        'TimeVec', timeVec, 'TimeUnits', 's', ...
                        'AmpMode', amplitudeDefinition, ...
                        'FundFreqRange', fundFreqRange, ...
                        'FilterCutoffsRelToFund', fCutoffRelToFund, ...
                        'FilterOrder', filterOrder, ...
                        'PromThresholdPerc', promThresholdPerc, ...
                        'MinPeakDistanceMs', minPeakDistanceMs);

    % Store the tables
    sniffPeakTables{iSwp} = sniffPeakTable;
    valleyTables{iSwp} = sniffValleyTable;

    % Store the fundamental frequency if it was computed
    if isfield(otherResults, 'freqFundamental')
        sniffFreqFundamental{iSwp} = otherResults.freqFundamental;
    else
        sniffFreqFundamental{iSwp} = NaN;
    end

    % If fewer than the minimum required peaks are found, we can't compute IPIs, so we skip.
    if height(sniffPeakTable) < minSniffPeaks
        sniffStartTimes{iSwp} = []; % Store an empty result for this sweep.
        sniffEndTimes{iSwp} = []; % Store an empty result for this sweep.
        continue; % Move to the next sweep.
    end

    % Extract the peak times from the resulting table.
    peakTimes = sniffPeakTable.peakTime;

    % Extract the times of the valleys that precede each peak.
    preValleyTimes = sniffPeakTable.preValleyTime;

    % Compute the inter-peak intervals (IPIs) between all sniff peaks.
    interPeakIntervals = diff(peakTimes);

    % Remove first and last peak from peaks to test
    preValleyTimesToTest = preValleyTimes(2:end-1);

    % Find the preceding inter-peak intervals for each peak 2:end
    preIPIs = interPeakIntervals(1:end-1);

    % Find the succeeding inter-peak intervals for each peak 2:end
    postIPIs = interPeakIntervals(2:end);

    % Detect basal-respiration-to-sniffing transition times (start of sniffing).
    % This is defined as a long IPI (basal) followed by a short IPI (sniffing).
    isSniffStart = preIPIs > sniffIpiThresholdSec & postIPIs <= sniffIpiThresholdSec;
    sniffStartTimes{iSwp} = preValleyTimesToTest(isSniffStart);

    % Detect sniffing-to-basal-respiration transition times (end of sniffing).
    % This is defined as a short IPI (sniffing) followed by a long IPI (basal).
    isSniffEnd = preIPIs <= sniffIpiThresholdSec & postIPIs > sniffIpiThresholdSec;
    sniffEndTimes{iSwp} = preValleyTimesToTest(isSniffEnd);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [whiskPeakTables, whiskValleyTables, whiskFreqFundamental, sniffStartWinTable] = ...
    parse_whisk_vecs(whiskVecsThisFile, tVecsThisFile, nSweeps, ...
                    sniffStartTimesThisFile, sniffEndTimesThisFile, ...
                    amplitudeDefinition, fundFreqRange, fCutoffRelToFund, filterOrder, ...
                    promThresholdPerc, minPeakDistanceMs, nWhisksToAnalyze)
%% Parses whisk vectors for a single file and generates an sniff start window table

% Initialize cell arrays to store results for each sweep in this file.
whiskPeakTables = cell(nSweeps, 1);
whiskValleyTables = cell(nSweeps, 1);
whiskFreqFundamental = cell(nSweeps, 1);

% Define empty table structure for when no windows are found
emptyTable = table(zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), {}, {}, {}, {}, {}, {}, {}, {}, ...
                'VariableNames', {'sweepNumber', 'sniffStartWinStartTime', 'sniffStartWinEndTime', 'sniffStartTime', ...
                                  'whiskPeakTimes', 'whiskPeakValues', 'whiskPeakAmplitudes', ...
                                  'whiskPreValleyTimes', 'whiskPostValleyTimes', ...
                                  'whiskValleyTimes', 'whiskValleyValues', 'whiskLogDecrements'});

% If whisk vector is empty, return empty results
if isempty(whiskVecsThisFile)
    sniffStartWinTable = emptyTable;
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
                        'FilterCutoffsRelToFund', fCutoffRelToFund, ...
                        'FilterOrder', filterOrder, ...
                        'PromThresholdPerc', promThresholdPerc, ...
                        'MinPeakDistanceMs', minPeakDistanceMs);
    
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

% Now, process records serially to build the sniff start window table
allWindowsCell = {};
for iSwp = 1:nSweeps
    % Obtain whisk peak and valley tables for this sweep
    whiskPeakTable = whiskPeakTables{iSwp};
    whiskValleyTable = whiskValleyTables{iSwp};

    % Obtain whisk peak times and sniff period start and end times
    peakTimes = whiskPeakTable.peakTime;
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

    % Loop through all sniff periods
    for iSniff = 1:numel(sniffStartTimes)
        % Get current sniff period start and end times
        currentSniffStart = sniffStartTimes(iSniff);
        currentSniffEnd = sniffEndTimes(iSniff);

        % Find all peaks within the sniff period
        isInSniffPeriod = peakTimes >= currentSniffStart & ...
                            peakTimes <= currentSniffEnd;
        whiskPeaksInSniffPeriod = whiskPeakTable(isInSniffPeriod, :);

        % If there are at least nWhisksToAnalyze peaks
        %   within this sniff period, add to sniff start window
        %   with appropriate boundaries
        if height(whiskPeaksInSniffPeriod) >= nWhisksToAnalyze
            % Extract whisk peaks to analyze
            whiskPeaksToAnalyze = whiskPeaksInSniffPeriod(1:nWhisksToAnalyze, :);

            % Start of sniff start window is the earlier of the sniff start time 
            % and the pre-valley time of the 1st peak to analyze
            firstPeakToAnalyze = whiskPeaksToAnalyze(1, :);
            sniffStartWinStartTime = min([currentSniffStart, firstPeakToAnalyze.preValleyTime]);

            % End of sniff start window is the post-valley time of the 
            %   last peak to analyze
            lastPeakToAnalyze = whiskPeaksToAnalyze(end, :);
            sniffStartWinEndTime = lastPeakToAnalyze.postValleyTime;

            % Find all valleys within the sniff start window
            isValleyInWin = valleyTimes >= sniffStartWinStartTime & ...
                          valleyTimes <= sniffStartWinEndTime;
            whiskValleysToAnalyze = whiskValleyTable(isValleyInWin, :);

            % Compute the logarithmic decrements of successive whisk peak amplitudes
            peakAmplitudes = whiskPeaksToAnalyze.amplitude;
            if numel(peakAmplitudes) > 1
                logDecrements = log(peakAmplitudes(2:end) ./ peakAmplitudes(1:end-1));
            else
                logDecrements = nan(size(nWhisksToAnalyze - 1, 1));
            end
            
            % Create a one-row table for this window
            newWindow = table(iSwp, sniffStartWinStartTime, sniffStartWinEndTime, currentSniffStart, ...
                {whiskPeaksToAnalyze.peakTime}, {whiskPeaksToAnalyze.peakValue}, ...
                {peakAmplitudes}, {whiskPeaksToAnalyze.preValleyTime}, {whiskPeaksToAnalyze.postValleyTime}, ...
                {whiskValleysToAnalyze.valleyTime}, {whiskValleysToAnalyze.valleyValue}, {logDecrements}, ...
                'VariableNames', {'sweepNumber', 'sniffStartWinStartTime', 'sniffStartWinEndTime', 'sniffStartTime', ...
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
    sniffStartWinTable = emptyTable;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sniffStartWinTable = augment_sniffstart_windows(sniffStartWinTable, ...
                                        sniffPeakTablesThisFile, ...
                                        sniffValleyTablesThisFile)
%% Augments a single sniff start window table with sniff data for one file

% If no sniff start windows were found for this file, return the empty table
if isempty(sniffStartWinTable)
    % Define empty columns for the case where the input table is empty
    newCols = {'sniffPeakTimes', 'sniffPeakValues', 'sniffPeakAmplitudes', ...
                'sniffPreValleyTimes', 'sniffPostValleyTimes', ...
                'sniffValleyTimes', 'sniffValleyValues'};
    for i = 1:numel(newCols)
        sniffStartWinTable.(newCols{i}) = cell(0, 1);
    end

    return;
end

% Get the number of sniff start windows
nWindows = height(sniffStartWinTable);

% Initialize new columns as cell arrays
sniffPeakTimesInWin = cell(nWindows, 1);
sniffPeakValuesInWin = cell(nWindows, 1);
sniffPeakAmpsInWin = cell(nWindows, 1);
sniffPreValleyTimesInWin = cell(nWindows, 1);
sniffPostValleyTimesInWin = cell(nWindows, 1);
sniffValleyTimesInWin = cell(nWindows, 1);
sniffValleyValuesInWin = cell(nWindows, 1);

% Loop through each sniff start window in the table
for iWin = 1:nWindows
    % Get window start/end times and sweep number
    winStart = sniffStartWinTable.sniffStartWinStartTime(iWin);
    winEnd = sniffStartWinTable.sniffStartWinEndTime(iWin);
    swpNum = sniffStartWinTable.sweepNumber(iWin);

    % Get the sniff data for the corresponding sweep
    sniffPeakTable = sniffPeakTablesThisFile{swpNum};
    sniffValleyTable = sniffValleyTablesThisFile{swpNum};

    % Find sniff peaks within the window
    if ~isempty(sniffPeakTable)
        isPeakInWin = sniffPeakTable.peakTime >= winStart & ...
                      sniffPeakTable.peakTime <= winEnd;
        sniffPeakTimesInWin{iWin} = sniffPeakTable.peakTime(isPeakInWin);
        sniffPeakValuesInWin{iWin} = sniffPeakTable.peakValue(isPeakInWin);
        sniffPeakAmpsInWin{iWin} = sniffPeakTable.amplitude(isPeakInWin);
        sniffPreValleyTimesInWin{iWin} = sniffPeakTable.preValleyTime(isPeakInWin);
        sniffPostValleyTimesInWin{iWin} = sniffPeakTable.postValleyTime(isPeakInWin);
    end

    % Find sniff valleys within the window
    if ~isempty(sniffValleyTable)
        isValleyInWin = sniffValleyTable.valleyTime >= winStart & ...
                        sniffValleyTable.valleyTime <= winEnd;
        sniffValleyTimesInWin{iWin} = sniffValleyTable.valleyTime(isValleyInWin);
        sniffValleyValuesInWin{iWin} = sniffValleyTable.valleyValue(isValleyInWin);
    end
end

% Add the new data as columns to the table
sniffStartWinTable.sniffPeakTimes = sniffPeakTimesInWin;
sniffStartWinTable.sniffPeakValues = sniffPeakValuesInWin;
sniffStartWinTable.sniffPeakAmplitudes = sniffPeakAmpsInWin;
sniffStartWinTable.sniffPreValleyTimes = sniffPreValleyTimesInWin;
sniffStartWinTable.sniffPostValleyTimes = sniffPostValleyTimesInWin;
sniffStartWinTable.sniffValleyTimes = sniffValleyTimesInWin;
sniffStartWinTable.sniffValleyValues = sniffValleyValuesInWin;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [isUsedToAverage, fileNumsToAverage] = ...
                identify_files_for_averaging(trialNames, excludeStringsFromAverage)
%% Identifies files for averaging by excluding specific trial names.

% Initialize as all true
isUsedToAverage = true(size(trialNames));

% Iteratively apply exclusion criteria
for i = 1:numel(excludeStringsFromAverage)
    isUsedToAverage = isUsedToAverage & ~contains(trialNames, excludeStringsFromAverage{i});
end

% Find the file numbers (indices) to be used for averaging
fileNumsToAverage = find(isUsedToAverage);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combinedTable = combine_sniffstart_windows(sniffStartWinTablesAll, ...
                                        fileNumsToAverage, trialNames)
%% Combines sniff start window tables from selected files into one large table.

% Initialize a cell array to hold tables that will be concatenated
tablesToCombine = {};

% Loop through each file number selected for averaging
for i = 1:numel(fileNumsToAverage)
    fileNum = fileNumsToAverage(i);
    
    % Get the sniff start window table for the current file
    currentTable = sniffStartWinTablesAll{fileNum};
    
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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function average_log_decrement_jitter(sniffStartWinTableToAverage, pathOutDir, figTypes, figTitle4, figName4)
%% Averages and plots whisk logarithmic decrements as a grouped jitter plot,
%  with mean, 95% CI, and statistical test results overlaid.

fprintf('Plotting aggregated whisk log decrements jitter plot...\n');

% If there's no data to plot, exit early
if isempty(sniffStartWinTableToAverage)
    fprintf('No sniff start windows to average. Skipping jitter plot.\n');
    return;
end

% Extract the necessary columns from the combined table
whiskLogDecrementsCell = sniffStartWinTableToAverage.whiskLogDecrements;
fileNumbersEachWin = sniffStartWinTableToAverage.fileNumber;

% Check if there is anything to plot after unpacking
if isempty(whiskLogDecrementsCell) || isempty(whiskLogDecrementsCell{1})
    fprintf('No log decrement data found to plot.\n');
    return;
end

% Convert cell array to a matrix, with each column being a decrement number
%   and each row being an sniff start window
allLogDecrementsMatrix = transpose(horzcat(whiskLogDecrementsCell{:}));

% Count the number of sniff start windows
nAnalysisWin = size(allLogDecrementsMatrix, 1);

% Count the number of log decrements
nLogDecrements = size(allLogDecrementsMatrix, 2);

% Create matching decrement numbers and file numbers
allDecrementNumbersMatrix = repmat((1:nLogDecrements), nAnalysisWin, 1);
allFileNumbersMatrix = repmat(fileNumbersEachWin, 1, nLogDecrements);

%% Compute statistics
% Calculate mean and 95% confidence intervals
meanWhiskLogDecrements = compute_combined_trace(allLogDecrementsMatrix', 'mean');
lower95WhiskLogDecrements = compute_combined_trace(allLogDecrementsMatrix', 'lower95');
upper95WhiskLogDecrements = compute_combined_trace(allLogDecrementsMatrix', 'upper95');

% Perform significance tests for each decrement
stats = vecfun(@test_difference, allLogDecrementsMatrix);

% Extract results
pValues = extract_fields(stats, 'pValue');
testFunctions = extract_fields(stats, 'testFunction');
symbols = extract_fields(stats, 'symbol');

%% Prepare data for the jitter plot
allLogDecrements = allLogDecrementsMatrix(:);
allDecrementNumbers = allDecrementNumbersMatrix(:);
allFileNumbers = allFileNumbersMatrix(:);
xTickLabels = arrayfun(@(x) sprintf('ln(A%d/A%d)', x+1, x), 1:nLogDecrements, 'UniformOutput', false);

%% Plot the jitter plot and overlay stats
% Set up figure properties
fig4 = set_figure_properties('FigNumber', 4, 'AlwaysNew', false, 'ClearFigure', true);
figPath = fullfile(pathOutDir, figName4);

% 1. Generate the base jitter plot, ensuring it doesn't plot its own stats
plot_grouped_jitter(allLogDecrements, allFileNumbers, allDecrementNumbers, ...
    'XTickLabels', xTickLabels, ...
    'YLabel', 'Log Decrement (ln(A_{n+1}/A_{n}))', ...
    'LegendLocation', 'suppress', ...
    'PlotMeanValues', false, 'PlotErrorBars', false, ...
    'RunTTest', false, 'RunRankTest', false);

% 2. Hold the plot to overlay new elements
hold on;

% 3. Plot the mean and 95% CI error bars with a distinct style
xValues = 1:nLogDecrements;
errLower = meanWhiskLogDecrements - lower95WhiskLogDecrements;
errUpper = upper95WhiskLogDecrements - meanWhiskLogDecrements;

errorbar(xValues, meanWhiskLogDecrements, errLower, errUpper, '_', ...
         'Color', 'k', 'LineWidth', 2.5, 'CapSize', 20, 'Marker', 'none');
plot(xValues, meanWhiskLogDecrements, 'o', 'MarkerEdgeColor', 'k', ...
     'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 1.5);

% 4. Add a horizontal line at y=0 for the null hypothesis
yline(0, '--k', 'LineWidth', 1);

% 5. Annotate plot with symbol and p-values near the top of the plot area
plot_test_result(pValues, 'TestFunction', testFunctions, 'Symbol', symbols, ...        
                 'XLocText', xValues, 'YLocTextRel', 0.90, 'YLocStarRel', 0.95);

% 6. Finalize plot
title(figTitle4);
grid on;
hold off;

% 7. Save the figure
save_all_figtypes(fig4, figPath, figTypes);

fprintf('Finished plotting jitter plot with mean/CI overlay.\n\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:



%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%