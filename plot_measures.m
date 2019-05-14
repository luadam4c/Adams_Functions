% plot_measures.m
%% Plots all measures of interest across slices

% Requires:
%       cd/combine_variables_across_tables.m
%       cd/count_vectors.m
%       cd/create_indices.m
%       cd/extract_common_directory.m
%       cd/extract_fileparts.m
%       cd/force_matrix.m
%       cd/ismatch.m
%       cd/plot_table.m
% Used by:
%       cd/parse_all_multiunit.m

% File History:
% 2019-03-15 Created by Adam Lu
% 2019-03-25 Now colors by phase number
% 2019-04-08 Renamed as plot_measures.m
% TODO: Add normalizeToBaseline as an optional argument
% TODO: extract specific usage to clc2_analyze.m
% 

%% Hard-coded parameters
% Flags (to be made arguments later)
normalizeToBaseline = true; %false;
plotType = 'tuning';

% Protocol parameters
sweepLengthSec = 60;
timeLabel = 'Time';
phaseLabel = 'Phase';
phaseStrs = {'Baseline', 'Wash-on', 'Wash-out'};

% Analysis parameters
%   Note: must be consistent with plot_struct.m
nSweepsLastOfPhase = 10;
nSweepsToAverage = 5;
maxRange2Mean = 40;

% File patterns
sliceFilePattern = '.*slice.*';
outFolder = pwd;

% Must be consistent with parse_multiunit.m
varsToPlot = {'oscIndex1'; 'oscIndex2'; 'oscIndex3'; 'oscIndex4'; ...
                    'oscPeriod1Ms'; 'oscPeriod2Ms'; ...
                    'oscDurationSec'; ...
                    'nSpikesTotal'; 'nSpikesIn10s'; 'nSpikesInOsc'; ...
                    'nBurstsTotal'; 'nBurstsIn10s'; 'nBurstsInOsc'; ...
                    'nSpikesPerBurst'; 'nSpikesPerBurstIn10s'; ...
                    'nSpikesPerBurstInOsc'};
varLabels = {'Oscillatory Index 1'; 'Oscillatory Index 2'; ...
                'Oscillatory Index 3'; 'Oscillatory Index 4'; ...
                'Oscillation Period 1 (ms)'; 'Oscillation Period 2 (ms)'; ...
                'Oscillation Duration (s)'; ...
                'Total Spike Count'; 'Number of Spikes in First 10 s'; 
                'Number of Spikes in Oscillation'; ...
                'Total Number of Bursts'; 'Number of Bursts in First 10 s'; ...
                'Number of Bursts in Oscillation'; ...
                'Number of Spikes Per Burst'; ...
                'Number of Spikes Per Burst in First 10 s'; ...
                'Number of Spikes Per Burst in Oscillation'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Find all files with the pattern *slice*_params in the file name
[~, sliceParamSheets] = all_files('Keyword', 'slice', 'Suffix', 'params', ...
                                    'ForceCellOutput', true);
                                
% Extract the common prefix
prefix = extract_fileparts(sliceParamSheets, 'commonprefix');

% If no common prefix, use the directory name as the prefix
if isempty(prefix)
    prefix = extract_common_directory(sliceParamSheets, 'BaseNameOnly', true);
end

% Extract the distinct parts of the file names
fileLabels = extract_fileparts(sliceParamSheets, 'distinct');

% Create table labels
tableLabels = strcat(prefix, {': '}, varLabels);

% Create table names
tableNames = strcat(prefix, '_', varsToPlot);

% Create figure names
figNames = fullfile(outFolder, strcat(tableNames, '.png'));
figNamesByPhase = fullfile(outFolder, strcat(tableNames, '_byphase.png'));
figNamesAvgd = fullfile(outFolder, strcat(tableNames, '_averaged.png'));

% Create paths for averaged tables
avgdTablePaths = fullfile(outFolder, strcat(tableNames, '_averaged.csv'));

% Read all slice parameter spreadsheets
sliceParamsTables = cellfun(@readtable, sliceParamSheets, ...
                            'UniformOutput', false);

% Create a time column (time in minutes since drug on)
sliceParamsTables = ...
    cellfun(@(x) create_time_rel_to_drugon(x, sweepLengthSec), ...
                sliceParamsTables, 'UniformOutput', false);

% Combine with phase number information
varsToCombine = [varsToPlot, repmat({'phaseNumber'}, size(varsToPlot))];

% Create the phaseNumber variables for the combined tables
phaseVars = strcat('phaseNumber_', fileLabels);

% Combine variables across tables
measureTables = combine_variables_across_tables(sliceParamsTables, ...
                'Keys', 'Time', 'VariableNames', varsToCombine, ...
                'InputNames', fileLabels, 'OmitVarName', false, ...
                'OutFolder', outFolder, 'Prefix', prefix, 'SaveFlag', true);

% Average over the last nSweepsToAverage sweeps
chevronTables = ...
    cellfun(@(x, y) average_last_of_each_phase(x, nSweepsLastOfPhase, ...
                                nSweepsToAverage, maxRange2Mean, y), ...
                    measureTables, avgdTablePaths, 'UniformOutput', false);

% Generate figure titles for Chevron plots
figTitlesChevron = strcat(varLabels, [' avg over ', ...
                            num2str(nSweepsToAverage), ' of last ', ...
                            num2str(nSweepsLastOfPhase), ' sweeps']);
                
% Normalize to baseline if requested
if normalizeToBaseline
    chevronTables = ...
        cellfun(@(x) normalize_to_first_row(x), ...
                        chevronTables, 'UniformOutput', false);

    varLabelsChevron = repmat({'% of baseline'}, size(varLabels));
else
    varLabelsChevron = varLabels;
end

%% Do the job
% Convert to timetables
measureTimeTables = cellfun(@table2timetable, ...
                            measureTables, 'UniformOutput', false);

% Plot Chevron plots
figs = cellfun(@(x, y, z, w, v, u) plot_table(x, 'PlotSeparately', false, ...
                                'PlotType', plotType, ...
                                'VariableNames', strcat(y, '_', fileLabels), ...
                                'PTicks', [1, 2, 3], 'PTickLabels', phaseStrs, ...
                                'ReadoutLabel', z, 'TableLabel', w, ...
                                'PLabel', phaseLabel, ...
                                'FigTitle', v, 'FigName', u, ...
                                'LegendLocation', 'eastoutside', ...
                                'RemoveOutliers', false), ...
            chevronTables, varsToPlot, varLabelsChevron, ...
            tableLabels, figTitlesChevron, figNamesAvgd);

close all;

% Plot all columns together
figs = cellfun(@(x, y, z, w, v) plot_table(x, 'PlotSeparately', false, ...
                                'PlotType', plotType, ...
                                'VariableNames', strcat(y, '_', fileLabels), ...
                                'ReadoutLabel', z, 'TableLabel', w, ...
                                'PLabel', timeLabel, 'FigName', v, ...
                                'RemoveOutliers', false), ...
                measureTimeTables, varsToPlot, varLabels, tableLabels, figNames);

close all;

% Plot all columns together colored by phase
figs = cellfun(@(x, y, z, w, v) plot_table(x, 'PlotSeparately', false, ...
                                'PlotType', plotType, ...
                                'VariableNames', strcat(y, '_', fileLabels), ...
                                'PhaseVariables', phaseVars, ...
                                'ReadoutLabel', z, 'TableLabel', w, ...
                                'PLabel', timeLabel, 'FigName', v, ...
                                'RemoveOutliers', false), ...
        measureTimeTables, varsToPlot, varLabels, tableLabels, figNamesByPhase);

%                                 'PhaseLabels', phaseStrs, ...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function myTable = create_time_rel_to_drugon(myTable, sweepLengthSec)
%% Creates a time column in minutes since drug onset

% Count the number of rows
nRows = height(myTable);

% Get the phase numbers
phaseNum = myTable.phaseNumber;

% Get the first row that is drug on
%   Note: this is set #2
rowDrugOn = find(phaseNum == 2, 1, 'first');

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

function outTable = average_last_of_each_phase(inTable, nSweepsLastOfPhase, ...
                                nSweepsToAverage, maxRange2Mean, sheetPath)
%% Averages over the last nSweepsToAverage sweeps of each phase
% Note: all distinct identifiers must have a matching phase variable column
% TODO: Fix

% Get all variable names
allVarNames = inTable.Properties.VariableNames;

% Remove the Time column if exists
% TODO FOR UNDERGRAD: is_variable_of_table.m
if ismember('Time', allVarNames)
    inTable = removevars(inTable, 'Time');
end

% Get all variable names that are left
allVarNames = transpose(inTable.Properties.VariableNames);

% Find the phase variable names
[~, phaseVars] = find_in_strings('phase.*', allVarNames, 'SearchMode', 'reg');

% Collect the rest of the variable names
readoutVars = setdiff(allVarNames, phaseVars);

% Extract the unique identifiers
uniqueIds = extract_distinct_fileparts(readoutVars, 'Delimiter', '_');

% Find matching phase variables for each unique identifiers
[~, phaseVars] = ...
    cellfun(@(x) find_in_strings(x, phaseVars, 'SearchMode', 'substrings'), ...
            uniqueIds, 'UniformOutput', false);

% Find the row indices for the last nSweepsLastOfPhase sweeps for each phase,
%   for each unique phase ID
indLastOfPhase = cellfun(@(x) find_last_ind_each_phase(inTable{:, x}, ...
                                                nSweepsLastOfPhase), ...
                    phaseVars, 'UniformOutput', false);

% Count the number of phases for each unique phase ID
nPhases = count_vectors(indLastOfPhase);

% Compute the maximum number of phases
maxNPhases = max(nPhases);

% Creat a phase number column
phaseNumber = transpose(1:maxNPhases);

% Average the nSweepsToAverage sweeps within maxRange2Mean in
%    the last nSweepsLastOfPhase sweeps of each phase
readoutAvg = ...
    cellfun(@(x, y) cellfun(@(z) compute_phase_average(inTable{:, x}, ...
                                        'Indices', z, ...
                                        'NToAverage', nSweepsToAverage, ...
                                        'MaxRange2Mean', maxRange2Mean), ...
                            y), ...
            readoutVars, indLastOfPhase, 'UniformOutput', false);

% Force as a matrix, padding each vector with NaNs at the end if necessary
readoutAvg = force_matrix(readoutAvg, 'AlignMethod', 'leftAdjustPad');

% Create an averaged table
outTable = array2table(readoutAvg);
outTable = addvars(outTable, phaseNumber, 'Before', 1);
outTable.Properties.VariableNames = vertcat({'phaseNumber'}, readoutVars);

% Save the table
writetable(outTable, sheetPath);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indLastEachPhase = find_last_ind_each_phase(phaseVec, nLastIndices)
%% Find the last nLastIndices indices for each phase in the phaseVec

% Get the unique phases
% TODO: use unique_custom.m
phaseVecNoNaN = phaseVec(~isnan(phaseVec));
uniquePhases = unique(phaseVecNoNaN, 'stable');

% Find the last index for each phase
lastIndexEachPhase = arrayfun(@(x) find(ismatch(phaseVec, x), 1, 'last'), ...
                                uniquePhases);

% Compute the first index for each phase
firstIndexEachPhase = lastIndexEachPhase - nLastIndices + 1;

% Construct the last nLastIndices indices
indLastEachPhase = create_indices('IndexStart', firstIndexEachPhase, ...
                                    'IndexEnd', lastIndexEachPhase, ...
                                    'MaxNum', nLastIndices, ...
                                    'ForceCellOutput', true);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function normalizedTable = normalize_to_first_row(table)
%% Normalizes all values to the first row

% Extract the first row
firstRow = table{1, :};

% Count the number of rows
nRows = height(table);

% Divide each row by the first row and multiply by 100 %
% normalizedTable = rowfun(@(x) x ./ firstRow, table);
for iRow = 1:nRows
    table{iRow, :} = (table{iRow, :} ./ firstRow) .* 100;
end

normalizedTable = table;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Find the row indices for the last nSweepsToAverage sweeps for each phase,
%   for each unique phase ID
indToAvg = cellfun(@(x) find_last_ind_each_phase(inTable{:, x}, nSweepsToAverage), ...
                    phaseVars, 'UniformOutput', false);
% Average the last nSweepsToAverage sweeps for each phase
readoutAvg = cellfun(@(x, y) cellfun(@(z) nanmean(inTable{:, x}(z)), y), ...
                    readoutVars, indToAvg, 'UniformOutput', false);

% Count the number of unique phases
nPhases = numel(uniquePhases);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%