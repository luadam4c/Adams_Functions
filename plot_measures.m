% plot_measures.m
%% Plots all measures of interest across slices

% Requires:
%       cd/argfun.m
%       cd/combine_variables_across_tables.m
%       cd/compute_phase_average.m
%       cd/compute_stats.m
%       cd/count_vectors.m
%       cd/create_indices.m
%       cd/extract_common_directory.m
%       cd/extract_fileparts.m
%       cd/force_matrix.m
%       cd/ismatch.m
%       cd/plot_table.m
%       cd/renamevars.m
%       cd/unique_custom.m
%
% Used by:
%       cd/parse_all_multiunit.m

% File History:
% 2019-03-15 Created by Adam Lu
% 2019-03-25 Now colors by phase number
% 2019-04-08 Renamed as plot_measures.m
% 2019-06-11 Now plots both normalized and original versions of the Chevron Plot
% 2019-08-06 Now plots markers for the Chevron plots
% 2019-08-06 Added phaseNumbers
% 2019-08-06 Renamed '_averaged' -> '_chevron'
% TODO: extract specific usage to clc2_analyze.m
% 

%% Hard-coded parameters
% Flags (TODO: make optional arguments)
plotType = 'tuning';
% TODO: plotAllFlag
plotChevronFlag = true;
plotNormalizedFlag = true;
plotByFileFlag = true;
plotByPhaseFlag = true;
plotPopAverageFlag = true;
computeChevronFlag = false;
computeNormalizedFlag = false;
computeTimeTablesFlag = false;
computePopAverageFlag = false;

% Protocol parameters
sweepLengthSec = 60;
timeLabel = 'Time';
phaseLabel = 'Phase';
phaseStrs = {'Baseline', 'Wash-on', 'Wash-out'};

% Analysis parameters
phaseNumbers = [1, 2];
% phaseNumbers = [];

%   Note: must be consistent with parse_multiunit.m
nSweepsLastOfPhase = 10;
nSweepsToAverage = 5;
maxRange2Mean = 40;
% maxRange2Mean = 20;
% maxRange2Mean = 200;

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
% Deal with flags
if plotPopAverageFlag
    computePopAverageFlag = true;
end
if plotByFileFlag || plotByPhaseFlag || computePopAverageFlag
    computeTimeTablesFlag = true;
end
if plotNormalizedFlag
    computeNormalizedFlag = true;
end
if plotChevronFlag || computeNormalizedFlag
    computeChevronFlag = true;
end

fprintf('Finding all spreadsheets ...\n');
% Find all files with the pattern *slice*_params in the file name
[~, sliceParamSheets] = all_files('Keyword', 'slice', 'Suffix', 'params', ...
                                    'ForceCellOutput', true);
                                
% Extract the common prefix
prefix = extract_fileparts(sliceParamSheets, 'commonprefix');

% If no common prefix, use the directory name as the prefix
if isempty(prefix)
    prefix = extract_common_directory(sliceParamSheets, 'BaseNameOnly', true);
end

% Modify prefix if necessary
if ~isempty(phaseNumbers)
    % Concatenate phase numbers into a string
    phaseNumbersString = num2str(phaseNumbers, '%d');

    % Append the phase numbers to the prefix
    prefix = [prefix, '_phase', phaseNumbersString];
end

% Extract the distinct parts of the file names
fileLabels = extract_fileparts(sliceParamSheets, 'distinct');

% Create table labels
tableLabels = strcat(prefix, {': '}, varLabels);

% Create table names
tableNames = strcat(prefix, '_', varsToPlot);

% Create figure names
[figNames, figNamesByPhase, figNamesAvgd, figNamesNormAvgd, figNamesPopAvg] = ...
    argfun(@(x) fullfile(outFolder, strcat(tableNames, '_', x, '.png')), ...
            'byfile', 'byphase', 'chevron', 'normChevron', 'popAverage');

% Create paths for chevron tables
% TODO: Use argfun
chevronTablePaths = fullfile(outFolder, strcat(tableNames, '_chevron.csv'));
normChevronTablePaths = fullfile(outFolder, strcat(tableNames, '_normChevron.csv'));
popAvgTablePaths = fullfile(outFolder, strcat(tableNames, '_popAverage.csv'));

% Create paths for mat files
% TODO: Use argfun
measureTablesMatPath = fullfile(outFolder, [prefix, '_', 'measureTables.mat']);
chevronTablesMatPath = fullfile(outFolder, [prefix, '_', 'chevronTables.mat']);
normalizedChevronTablesMatPath = fullfile(outFolder, [prefix, '_', 'normalizedChevronTables.mat']);
measureTimeTablesMatPath = fullfile(outFolder, [prefix, '_', 'measureTimeTables.mat']);
popAvgTablesMatPath = fullfile(outFolder, [prefix, '_', 'popAvgTables.mat']);

%% Read in data and preprocess
% Read all slice parameter spreadsheets
fprintf('Reading measure spreadsheets ...\n');
sliceParamsTables = cellfun(@readtable, sliceParamSheets, ...
                            'UniformOutput', false);

% Create a time column (time in minutes since drug on)
fprintf('Creating time columns ...\n');
sliceParamsTables = ...
    cellfun(@(x) create_time_rel_to_drugon(x, sweepLengthSec), ...
                sliceParamsTables, 'UniformOutput', false);

% Combine with phase number information
%   Note: each row results in a different table
varsToCombine = [varsToPlot, repmat({'phaseNumber'}, size(varsToPlot))];

% Create the phaseNumber variables for the combined tables
phaseVars = strcat('phaseNumber_', fileLabels);

% Restrict to certain phases if requested
if ~isempty(phaseNumbers)
    fprintf('Restricting to phases %s ...\n', num2str(phaseNumbers, '%d, '));
    sliceParamsTables = ...
        cellfun(@(x) x(ismember(x.phaseNumber, phaseNumbers), :), ...
                sliceParamsTables, 'UniformOutput', false);
end

% Generate phase tick locations
if isempty(phaseNumbers)
    phaseTickLocs = 1:3;
else
    phaseTickLocs = 1:numel(phaseNumbers);
    phaseStrs = phaseStrs(phaseNumbers);
end

% Combine variables across tables
fprintf('Combining variables across tables ...\n');
measureTables = combine_variables_across_tables(sliceParamsTables, ...
                'Keys', 'Time', 'VariableNames', varsToCombine, ...
                'InputNames', fileLabels, 'OmitVarName', false, ...
                'OutFolder', outFolder, 'Prefix', prefix, 'SaveFlag', true);
save(measureTablesMatPath, 'measureTables', '-mat');

%% Compute population data
%% Average over each phase
if computeChevronFlag
    % Average over the last nSweepsToAverage sweeps
    fprintf('Creating Chevron tables ...\n');
    chevronTables = ...
        cellfun(@(x, y) average_last_of_each_phase(x, nSweepsLastOfPhase, ...
                                    nSweepsToAverage, maxRange2Mean, y), ...
                        measureTables, chevronTablePaths, 'UniformOutput', false);

    % Generate figure titles for Chevron plots
    figTitlesChevron = strcat(varLabels, [' avg over ', ...
                                num2str(nSweepsToAverage), ' of last ', ...
                                num2str(nSweepsLastOfPhase), ' sweeps']);


    % Generate variable labels
    varLabelsChevron = varLabels;

    save(chevronTablesMatPath, 'chevronTables', '-mat');
end

%% Normalize to baseline if requested
if computeNormalizedFlag
    fprintf('Normalizing to baseline ...\n');
    normalizedChevronTables = ...
        cellfun(@(x, y) normalize_to_first_row(x, y), ...
                    chevronTables, normChevronTablePaths, 'UniformOutput', false);

    % Generate variable labels
    varLabelsNormAvgd = repmat({'% of baseline'}, size(varLabels));

    save(normalizedChevronTablesMatPath, 'normalizedChevronTables', '-mat');
end

%% Convert to timetables
if computeTimeTablesFlag
    fprintf('Converting measure tables to time tables ...\n');
    measureTimeTables = cellfun(@table2timetable, ...
                                measureTables, 'UniformOutput', false);

    save(measureTimeTablesMatPath, 'measureTimeTables', '-mat');
end

%% Average over each file
if computePopAverageFlag
    fprintf('Computing population averages ...\n');
    popAvgTables = cellfun(@(x, y, z) compute_population_average(x, y, z), ...
                            measureTimeTables, varsToPlot, popAvgTablePaths, ...
                            'UniformOutput', false);

    save(popAvgTablesMatPath, 'popAvgTables', '-mat');
end

%% Do the job
if plotChevronFlag
    close all;

    % Plot Chevron plots
    fprintf('Plotting Chevron plots ...\n');
    figs = cellfun(@(x, y, z, w, v, u) plot_table(x, 'PlotSeparately', false, ...
                                    'PlotType', plotType, ...
                                    'VariableNames', strcat(y, '_', fileLabels), ...
                                    'PTicks', phaseTickLocs, ...
                                    'PTickLabels', phaseStrs, ...
                                    'ReadoutLabel', z, 'TableLabel', w, ...
                                    'PLabel', phaseLabel, ...
                                    'FigTitle', v, 'FigName', u, ...
                                    'LegendLocation', 'eastoutside', ...
                                    'RemoveOutliers', false, ...
                                    'Marker', 'o', 'MarkerFaceColor', 'auto'), ...
                chevronTables, varsToPlot, varLabelsChevron, ...
                tableLabels, figTitlesChevron, figNamesAvgd);
end

%%
if plotNormalizedFlag
    close all;

    % Plot Normalized Chevron plots
    fprintf('Plotting Normalized Chevron plots ...\n');
    figs = cellfun(@(x, y, z, w, v, u) plot_table(x, 'PlotSeparately', false, ...
                                    'PlotType', plotType, ...
                                    'VariableNames', strcat(y, '_', fileLabels), ...
                                    'PTicks', phaseTickLocs, ...
                                    'PTickLabels', phaseStrs, ...
                                    'ReadoutLabel', z, 'TableLabel', w, ...
                                    'PLabel', phaseLabel, ...
                                    'FigTitle', v, 'FigName', u, ...
                                    'LegendLocation', 'eastoutside', ...
                                    'RemoveOutliers', false, ...
                                    'Marker', 'o', 'MarkerFaceColor', 'auto'), ...
                normalizedChevronTables, varsToPlot, varLabelsNormAvgd, ...
                tableLabels, figTitlesChevron, figNamesNormAvgd);
end

%%
if plotByFileFlag
    close all;

    % Plot each column separately
    fprintf('Plotting each column with a different color ...\n');
    figs = cellfun(@(x, y, z, w, v) plot_table(x, 'PlotSeparately', false, ...
                                    'PlotType', plotType, ...
                                    'VariableNames', strcat(y, '_', fileLabels), ...
                                    'ReadoutLabel', z, 'TableLabel', w, ...
                                    'PLabel', timeLabel, 'FigName', v, ...
                                    'RemoveOutliers', false), ...
                    measureTimeTables, varsToPlot, varLabels, tableLabels, figNames);
end

%%
if plotByPhaseFlag
    close all;

    % Plot all columns together colored by phase
    fprintf('Plotting each phase with a different color ...\n');
    figs = cellfun(@(x, y, z, w, v) plot_table(x, 'PlotSeparately', false, ...
                                    'PlotType', plotType, ...
                                    'VariableNames', strcat(y, '_', fileLabels), ...
                                    'PhaseVariables', phaseVars, ...
                                    'PhaseLabels', phaseStrs, ...
                                    'ReadoutLabel', z, 'TableLabel', w, ...
                                    'PLabel', timeLabel, 'FigName', v, ...
                                    'RemoveOutliers', false), ...
            measureTimeTables, varsToPlot, varLabels, tableLabels, figNamesByPhase);

    %                                 'PhaseLabels', phaseStrs, ...
end

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
uniquePhases = unique_custom(phaseVec, 'stable', 'IgnoreNaN', true);

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

function normalizedTable = normalize_to_first_row(table, sheetPath)
%% Normalizes all values to the first row

% Extract the first row
firstRow = table{1, :};
% firstRow = table{3, :};

% Count the number of rows
nRows = height(table);

% Divide each row by the first row and multiply by 100 %
% normalizedTable = rowfun(@(x) x ./ firstRow, table);
for iRow = 1:nRows
    table{iRow, :} = (table{iRow, :} ./ firstRow) .* 100;
end

normalizedTable = table;

% Save the table
writetable(normalizedTable, sheetPath);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outTimeTable = compute_population_average (inTimeTable, varName, sheetPath)
%% Computes the population mean and confidence intervals

% Get all the variable names
columnNames = inTimeTable.Properties.VariableNames;

% Select the columns that contains given variable name
columnsToAverage = contains(columnNames, varName);

% Extract the data
popData = inTimeTable{:, columnsToAverage};

% Compute the means and bounds of 95% confidence intervals
[means, lower95s, upper95s] = ...
    argfun(@(x) compute_stats(popData, x, 2), ...
            'mean', 'lower95', 'upper95');

% Save the row times from inTable
rowTimes = inTimeTable.Properties.RowTimes;

% Create output table
outTimeTable = timetable(rowTimes, means, lower95s, upper95s, ...
                    'VariableNames', {'Mean', 'Lower95', 'Upper95'});


% Convert to a time table
outTable = timetable2table(outTimeTable);

% Rename 'rowTimes' -> 'Time'
outTable = renamevars(outTable, 'rowTimes', 'Time');

% Save the table
writetable(outTable, sheetPath);

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

phaseStrs = {'baseline', 'washon', 'washoff'};

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
