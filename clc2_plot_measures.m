% clc2_plot_measures.m
%% Plots all measures of interest across slices

% Requires:
%       cd/combine_variables_across_tables.m
%       cd/extract_fileparts.m

% File History:
% 2019-03-15 Created by Adam Lu
% 

%% Hard-coded parameters
% Protocol parameters
sweepLengthSec = 60;

% File patterns
sliceFilePattern = '.*slice.*';
outFolder = pwd;
timeLabel = 'Time';

% Must be consistent with parse_multiunit.m
varsToPlot = {'oscIndex1'; 'oscIndex2'; 'oscIndex3'; 'oscIndex4'; ...
                'oscPeriod1Ms'; 'oscPeriod2Ms'; ...
                'oscDurationSec'; 'nSpikesTotal'};
varLabels = {'Oscillatory Index 1'; 'Oscillatory Index 2'; ...
                'Oscillatory Index 3'; 'Oscillatory Index 4'; ...
                'Oscillation Period 1 (ms)'; 'Oscillation Period 2 (ms)'; ...
                'Oscillation Duration (s)'; 'Total Spike Count'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Find all files with the pattern *slice*_params in the file name
[~, sliceParamSheets] = all_files('Keyword', 'slice', 'Suffix', 'params');

% Extract the common prefix
prefix = extract_fileparts(sliceParamSheets, 'commonprefix');

% Extract the distinct parts of the file names
fileLabels = extract_fileparts(sliceParamSheets, 'distinct');

% Read all slice parameter spreadsheets
sliceParamsTables = cellfun(@readtable, sliceParamSheets, ...
                            'UniformOutput', false);

% Create a time column (time in minutes since drug on)
sliceParamsTables = ...
    cellfun(@(x) create_time_rel_to_drugon(x, sweepLengthSec), ...
                sliceParamsTables, 'UniformOutput', false);

% Combine variables across tables
measureTables = combine_variables_across_tables(sliceParamsTables, ...
                'Keys', 'Time', 'VariableNames', varsToPlot, ...
                'InputNames', fileLabels, 'OmitVarName', false, ...
                'OutFolder', outFolder, 'Prefix', prefix, 'SaveFlag', true);

% Create table labels
tableLabels = strcat(prefix, {': '}, varLabels);

% Create figure names
figNames = fullfile(outFolder, strcat(prefix, '_', varsToPlot, '.png'));

%% Do the job
% Convert to timetables
measureTimeTables = cellfun(@table2timetable, ...
                            measureTables, 'UniformOutput', false);

% Plot all columns together
% TODO: Plot drug onset and offset
figs = cellfun(@(x, y, z, w) plot_table(x, 'PlotSeparately', false, ...
                                'ReadoutLabel', y, 'TableLabel', z, ...
                                'XLabel', timeLabel, 'FigName', w, ...
                                'RemoveOutliers', true), ...
                measureTimeTables, varLabels, tableLabels, figNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function myTable = create_time_rel_to_drugon(myTable, sweepLengthSec)
%% Creates a time column in minutes since drug onset

% Count the number of rows
nRows = height(myTable);

% Get the set numbers
setNum = myTable.setNumber;

% Get the first row that is drug on
%   Note: this is set #2
rowDrugOn = find(setNum == 2, 1, 'first');

% Compute time in sweeps relative to drug on
timeSwps = transpose(1:nRows) - rowDrugOn;

% Convert to minutes
timeMin = timeSwps / (sweepLengthSec / 60);

% Convert to a duration vector
Time = minutes(timeMin);

% Add a time column to the table
myTable = addvars(myTable, Time, 'Before', 1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%