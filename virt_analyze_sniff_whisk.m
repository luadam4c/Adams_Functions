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
%       cd/check_dir.m
%       cd/extract_fields.m
%       cd/extract_fileparts.m
%       cd/find_matching_files.m
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

% Plotting parameters
whiskAngleLimits = [-75, 75];           % Whisk angle limits to be plotted
whiskAngleLabel = 'Whisk Angle (degrees)';
thermLimits = [-7.55, 7.55] * 1e-5;     % Thermocouple limits to be plotted
sampleExpNum = 38;                      % The sample experiment number to plot (max 38)
timeLabel = 'Time (s)';
whiskLabel = 'Whisking';
sniffLabel = 'Breathing';
legendLocation1 = 'northeast';
figTitle1 = 'Sample whisking (blue) and breathing (red) data';
figName1 = ['sniff_whisk_all_traces_', num2str(sampleExpNum)];
figTypes = {'eps', 'png'};

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

%% Reformat data for ease of plotting
% Find the maximum and minimum whisk angles
maxWhiskAngle = apply_iteratively(@max, whiskVecs); % == 72.5
minWhiskAngle = apply_iteratively(@min, whiskVecs); % == -71.9

% Find the maximum and minimum thermocouple values
maxThermValue = apply_iteratively(@max, thermVecs); % == 7.54e-05
minThermValue = apply_iteratively(@min, thermVecs); % == -6.88e-05

% Create sniff/breath vectors scaled to the order of whisk angle vectors
whiskAngleLimRange = diff(whiskAngleLimits);
thermLimRange = diff(thermLimits);
sniffVecs = cellfun(@(x) whiskAngleLimits(1) + (x - thermLimits(1)) * ...
                    whiskAngleLimRange / thermLimRange, thermVecs, ...
                    'UniformOutput', false);

%% Plot all data
% Extract data to plot
tVecsThis = tVecs{sampleExpNum};
whiskVecsThis = whiskVecs{sampleExpNum};
sniffVecsThis = sniffVecs{sampleExpNum};

% Count the number of plots
nPlots = size(tVecsThis, 2);
iPlots = (1:nPlots)';

% Create labels for legend
whiskLabels = create_labels_from_numbers(iPlots, 'Prefix', whiskLabel);
sniffLabels = create_labels_from_numbers(iPlots, 'Prefix', sniffLabel);

% Create paths for saving
figPath1 = fullfile(pathOutDir, figName1);

% Create sample plot
figure(1); clf;
if ~isempty(sniffVecsThis)
    handles = plot_traces(tVecsThis, whiskVecsThis, 'ColorMap', 'b', ...
                'DataToCompare', sniffVecsThis, 'ColorMapToCompare', 'r', ...
                'TraceLabels', whiskLabels, 'TraceLabelsToCompare', sniffLabels, ...
                'PlotMode', 'parallel', 'FigTitle', figTitle1, ...
                'XLabel', timeLabel, 'LegendLocation', legendLocation1, ...
                'YLimits', whiskAngleLimits, 'YLabel', whiskAngleLabel, ...
                'FigName', figPath1, 'FigTypes', figTypes);
else
    handles = plot_traces(tVecsThis, whiskVecsThis, 'ColorMap', 'b', ...
                'TraceLabels', whiskLabels, 'TraceLabelsToCompare', sniffLabels, ...
                'PlotMode', 'parallel', 'FigTitle', figTitle1, ...
                'XLabel', timeLabel, 'LegendLocation', legendLocation1, ...
                'YLimits', whiskAngleLimits, 'YLabel', whiskAngleLabel, ...
                'FigName', figPath1, 'FigTypes', figTypes);
end

figure(2); clf;
plot_traces(tVecs{end}, pulseVecs{end}, 'PlotMode', 'parallel');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:



%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%