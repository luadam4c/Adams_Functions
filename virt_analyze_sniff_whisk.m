%% Analyzes rat sniff whisk data 
% 
% This script must be run in a directory that contains the data directory
% specified below.
%
% Dataset: sniff_whisk_data
% Description from Jeff: 
%   - unilateral whisking and sniffing data, head restrained
% 	- whisking data are from tracking a single vibrissa in full-field video
%   - breathing data are from a thermocouple
% 	- some trials contain puffs of ammonia or control puffs of unscented air to the nose
%       time of application is indicated by the pulse in channel "piezo"
% 	- check sign
%
% Requires:
%       cd/all_files.m
%       cd/argfun.m
%       cd/array_fun.m
%       cd/check_dir.m
%       cd/count_vectors.m
%       cd/create_labels_from_numbers.m
%       cd/extract_fields.m
%       cd/extract_fileparts.m
%       cd/find_matching_files.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/plot_traces.m
%       cd/print_cellstr.m
%
% Used by:

% File History:
% 2025-08-31 Created by Adam Lu
% 

%% Hard-coded parameters
% Directory and file naming conventions
nameDataDir = 'data_sniff_whisk';       % Name of the directory containing the data files
nameOutDir = 'output_sniff_whisk';      % Name of the output directory
fileNameMetaData = 'metadata.csv';      % File name of the meta data file
extDataFile = 'mat';                    % Extension for the data files
prefixDataFile = 'swdata_';             % Prefix for the data files
suffixDataFile = '_angle';              % Suffix for the data files
suffixParamFile = '_parameters';        % Suffix for the parameter files

% Detection parameters
ammoniaPuffString = 'ammpuff';
airPuffString = 'airpuff';
minPulseAmplitude = 2;                  % Minimum pulse amplitude for piezo trace

% Plotting parameters
sampleFileNumsToPlot = (1:38)';         % The sample file number(s) to plot (max 38)
toSpeedUp = true;                       % Whether to use parpool and hide figures
whiskAngleLimits = [-75, 75];           % Whisk angle limits to be plotted
piezoLimits = [-1, 10];
whiskAngleLabel = 'Whisk Angle (degrees)';
sniffToWhiskRangeRatio = 0.8;           % Sniff amplitude to whisk amplitude ratio for plotting
timeLabel = 'Time (s)';
whiskLabel = 'Whisking';
sniffLabel = 'Breathing';
legendLocation1 = 'northeast';
legendLocation2 = 'suppress';
figTitlePrefix1 = 'Whisking (blue) and breathing (red) data';
figTitlePrefix2 = 'Stimulation (Puff) data';
figPrefix1 = 'sniff_whisk_all_traces_';
figPrefix2 = 'sniff_whisk_all_stims_';
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
nTrials = count_vectors(tVecs);

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

% Detect the first pulse start and end times
pulseParams = cellfun(@(x, y) parse_pulse(x, 'TimeVecs', y, ...
        'MinPulseAmplitude', minPulseAmplitude), pulseVecs, tVecs, 'UniformOutput', false);
pulseStartTimes = cellfun(@(x) x.timeAfterStartMs, pulseParams, 'UniformOutput', false);
pulseEndTimes = cellfun(@(x) x.timeBeforeEndMs, pulseParams, 'UniformOutput', false);

%% Plot all data
[handles1, handles2] = ...
    array_fun(@(a) plot_one_file (a, tVecs, whiskVecs, sniffVecs, ...
        pulseVecs, pulseStartTimes, pulseEndTimes, ...
        isAmmoniaPuff, isAirPuff, whiskLabel, sniffLabel, ...
        pathOutDir, timeLabel, whiskAngleLimits, whiskAngleLabel, piezoLimits, ...
        legendLocation1, legendLocation2, figTitlePrefix1, figTitlePrefix2, ...
        figPrefix1, figPrefix2, figTypes, toSpeedUp), ...
        sampleFileNumsToPlot, 'UseParpool', toSpeedUp);

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
    metaData.nRecords = nTrials;
    metaData.isAmmoniaPuff = isAmmoniaPuff;
    metaData.isAirPuff = isAirPuff;

    % Expand the cell arrays to a single string
    metaDataToPrint = metaData;
    [metaDataToPrint.roiswitch, metaDataToPrint.records, ...
        metaDataToPrint.anesrecords] = ...
        argfun(@(x) cellfun(@(a) print_cellstr(a, 'Delimiter', ' ', ...
                                'ToPrint', false, 'OmitQuotes', true, ...
                                'OmitBraces', true, 'OmitNewline', true), ...
                                x, 'UniformOutput', false), ...
            metaData.roiswitch, metaData.records, metaData.anesrecords);

    % Display the resulting table in the command window
    fprintf('Generated Metadata Table:\n');
    disp(metaDataToPrint);

    % Save metadata table
    pathMetaData = fullfile(pathOutDir, fileNameMetaData);
    writetable(metaDataToPrint, pathMetaData);
    fprintf('Metadata table saved to %s\n', pathMetaData);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [handles1, handles2] = ...
    plot_one_file (sampleFileNum, tVecs, whiskVecs, sniffVecs, ...
        pulseVecs, pulseStartTimes, pulseEndTimes, ...
        isAmmoniaPuff, isAirPuff, whiskLabel, sniffLabel, ...
        pathOutDir, timeLabel, whiskAngleLimits, whiskAngleLabel, piezoLimits, ...
        legendLocation1, legendLocation2, figTitlePrefix1, figTitlePrefix2, ...
        figPrefix1, figPrefix2, figTypes, toSpeedUp)

% Extract data to plot
tVecsThis = tVecs{sampleFileNum};
whiskVecsThis = whiskVecs{sampleFileNum};
sniffVecsThis = sniffVecs{sampleFileNum};
pulseVecsThis = pulseVecs{sampleFileNum};
pulseStartTimesThis = pulseStartTimes{sampleFileNum};
pulseEndTimesThis = pulseEndTimes{sampleFileNum};
isAmmoniaPuffThis = isAmmoniaPuff(sampleFileNum);
isAirPuffThis = isAirPuff(sampleFileNum);

% Count the number of plots
nPlots = size(tVecsThis, 2);
iPlots = (1:nPlots)';

% Create labels for legend
whiskLabels = create_labels_from_numbers(iPlots, 'Prefix', whiskLabel);
sniffLabels = create_labels_from_numbers(iPlots, 'Prefix', sniffLabel);

% Decide on puff window colors and titles
if isAmmoniaPuffThis
    figTitlePrefix1 = replace(figTitlePrefix1, 'data', 'data with Ammonia Puffs (green)');
    figTitlePrefix2 = replace(figTitlePrefix2, 'Puff', 'Ammonia Puff');
    puffWindowColor = [0.5625, 0.9297, 0.5625];  % rgb('LightGreen')
elseif isAirPuffThis
    figTitlePrefix1 = replace(figTitlePrefix1, 'data', 'data with Control Air Puffs (gray)');
    figTitlePrefix2 = replace(figTitlePrefix2, 'Puff', 'Control Air Puff');
    puffWindowColor = 0.8242 * [1, 1, 1]; % rgb('LightGray')
else
    puffWindowColor = 0.8242 * [1, 1, 1]; % rgb('LightGray')
end

% Create figure titles
figTitle1 = [figTitlePrefix1, ' for File # ', num2str(sampleFileNum)];
figTitle2 = [figTitlePrefix2, ' for File # ', num2str(sampleFileNum)];

% Create paths for saving
figName1 = [figPrefix1, num2str(sampleFileNum)];
figName2 = [figPrefix2, num2str(sampleFileNum)];
figPath1 = fullfile(pathOutDir, figName1);
figPath2 = fullfile(pathOutDir, figName2);

% Create sample traces plot
fprintf('Plotting sample traces for file %d ...\n', sampleFileNum);
if toSpeedUp
    set_figure_properties('AlwaysNew', true, 'Visible', 'off');
else
    set_figure_properties('FigNumber', 1, 'AlwaysNew', false, 'ClearFigure', true);
end
if ~isempty(sniffVecsThis)
    handles1 = plot_traces(tVecsThis, whiskVecsThis, 'ColorMap', 'b', ...
                'DataToCompare', sniffVecsThis, 'ColorMapToCompare', 'r', ...
                'XBoundaries', [pulseStartTimesThis, pulseEndTimesThis], ...
                'XBoundaryColor', puffWindowColor, ...
                'TraceLabels', whiskLabels, 'TraceLabelsToCompare', sniffLabels, ...
                'PlotMode', 'parallel', 'FigTitle', figTitle1, ...
                'XLabel', timeLabel, 'LegendLocation', legendLocation1, ...
                'YLimits', whiskAngleLimits, 'YLabel', whiskAngleLabel, ...
                'FigName', figPath1, 'FigTypes', figTypes);
else
    handles1 = plot_traces(tVecsThis, whiskVecsThis, 'ColorMap', 'b', ...
                'TraceLabels', whiskLabels, 'TraceLabelsToCompare', sniffLabels, ...
                'XBoundaries', [pulseStartTimesThis, pulseEndTimesThis], ...
                'XBoundaryColor', puffWindowColor, ...
                'PlotMode', 'parallel', 'FigTitle', figTitle1, ...
                'XLabel', timeLabel, 'LegendLocation', legendLocation1, ...
                'YLimits', whiskAngleLimits, 'YLabel', whiskAngleLabel, ...
                'FigName', figPath1, 'FigTypes', figTypes);
end

% Plot stim traces with detection
if ~isempty(pulseVecsThis)
    fprintf('Plotting stim traces with detection for file %d ...\n', sampleFileNum);
    if toSpeedUp
        set_figure_properties('AlwaysNew', true, 'Visible', 'off');
    else
        set_figure_properties('FigNumber', 2, 'AlwaysNew', false, 'ClearFigure', true);
    end
    handles2 = plot_traces(tVecsThis, pulseVecsThis, 'ColorMap', 'k', ...
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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:



%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%