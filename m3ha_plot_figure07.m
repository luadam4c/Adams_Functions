% m3ha_plot_figure07.m
%% Plots Figure 07 for the GAT Blocker paper
%
% Requires:
%       cd/apply_over_cells.m
%       cd/archive_dependent_scripts.m
%       cd/create_label_from_sequence.m
%       cd/find_matching_files.m
%       cd/force_column_cell.m
%       cd/m3ha_network_plot_gabab.m
%       cd/m3ha_network_plot_essential.m
%       cd/m3ha_plot_violin.m
%       cd/match_positions.m
%       cd/match_row_count.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/update_figure_for_corel.m

% File History:
% 2020-01-30 Modified from m3ha_plot_figure05.m

%% Hard-coded parameters
% Flags
plotIpscComparison = true;
plot2CellExamples = true;

combine2CellPopulation = false; %true;
plot2CellViolins = false; %true;
plot2CellBars = false; %true;

archiveScriptsFlag = true;

% Directories
parentDirectory = fullfile('/media', 'adamX', 'm3ha');
figure07Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure07');
networkDirectory = fullfile(parentDirectory, 'network_model');
% exampleIterName = '20200131T1345_using_bestparams_20200126_singleneuronfitting101';  % 20200131
exampleIterName = 'TODO';
popIterName = '20200204T1042_using_bestparams_20200203_manual_singleneuronfitting0-102_vtraub_-65';
candCellSheetName = 'candidate_cells.csv';
oscParamsSuffix = 'oscillation_params';

% Well-fitted, good 2-cell network response
rankNumsToUse = [2, 4, 5, 7, 9, 10, 12, 13, 16, 20, 21, 23, 25, 29];

% Files

% Analysis settings
% Should be consistent with m3ha_plot_figure03.m & m3ha_plot_figure07.m
exampleCellNames = {'D101310'; 'G101310'};

gIncr = 200;            % Original dynamic clamp gIncr value
measuresOfInterest = {'oscDurationSec'; 'oscPeriod2Ms'; ...
                        'oscIndex4'; 'hasOscillation'};
measureTitles = {'Oscillation Duration (sec)'; 'Oscillation Period (ms)'; ...
                    'Oscillatory Index'; 'Has Oscillation'};

% Plot settings
ipscFigWidth = 8.5;
ipscFigHeight = 4;
exampleFigWidth = 8.5;
exampleFigHeight = 1.5 * 6;
pharmLabelsShort = {'{\it s}-Con', '{\it s}-GAT1', ...
                    '{\it s}-GAT3', '{\it s}-Dual'};

figTypes = {'png', 'epsc2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Find the directory for this iteration
exampleIterDir = fullfile(networkDirectory, exampleIterName);
popIterDir = fullfile(networkDirectory, popIterName);

% Construct the full path to the candidate cell spreadsheet
candCellSheetPath = fullfile(networkDirectory, candCellSheetName);

% Create a rank string
rankStr = ['rank', create_label_from_sequence(rankNumsToUse)];

% Create a condition label
conditionLabel2D = [popIterName, '_', rankStr, '_gIncr', num2str(gIncr)];

% Create a population data spreadsheet name
popDataSheetName = [popIterName, '_', rankStr, '_', oscParamsSuffix, '.csv'];

% Contruct the full path to the population data spreadsheet
popDataPath = fullfile(figure07Dir, popDataSheetName);

%% Find example files and directories
if plotIpscComparison || plot2CellExamples
    % Find example network directories
    [~, exampleDirs] = ...
        cellfun(@(x) all_subdirs('Directory', popIterDir, 'Keyword', x), ...
                exampleCellNames, 'UniformOutput', false);
end

%% Plots figures for comparing dynamic clamp ipsc
if plotIpscComparison
    cellfun(@(x, y) plot_ipsc_comparison(x, exampleIterName, gIncr, y, ...
                                        figure07Dir, figTypes, ...
                                        ipscFigWidth, ipscFigHeight), ...
            exampleCellNames, exampleDirs);
end

%% Plots example 2-cell networks
if plot2CellExamples
    arrayfun(@(z) ...
        cellfun(@(x, y) plot_2cell_examples(x, exampleIterName, gIncr, z, y, ...
                                        figure07Dir, figTypes, ...
                                        exampleFigWidth, exampleFigHeight), ...
                exampleCellNames, exampleDirs), ...
        pharm);
end

%% Combines quantification over all 2-cell networks
if combine2CellPopulation
    combine_2cell_data(popIterDir, candCellSheetPath, rankNumsToUse, popDataPath);
end

%% Plots oscillation duration and period over pharm condition 
%       across all 2-cell networks
if plot2CellViolins
    % Construct stats table path
    stats2dPath = fullfile(figure07Dir, strcat(conditionLabel2D, '_stats.mat'));

    % Compute statistics if not done already
    if ~isfile(stats2dPath)
        % Compute statistics for all features
        disp('Computing statistics for violin plots ...');
        statsTable = m3ha_network_compute_statistics(popDataPath, gIncr, ...
                                            measuresOfInterest, measureTitles);

        % Generate labels
        conditionLabel = conditionLabel2D;
        pharmLabels = pharmLabelsShort;

        % Save stats table
        save(stats2dPath, 'statsTable', 'pharmLabels', ...
                            'conditionLabel', '-v7.3');
    end

    % Plot violin plots
    m3ha_plot_violin(stats2dPath, 'RowsToPlot', measuresOfInterest, ...
                    'OutFolder', figure07Dir);
end

%% Plots percentage of 2-cell networks having oscillations over pharm condition
if plot2CellBars
end

%% Archive all scripts for this run
if archiveScriptsFlag
    archive_dependent_scripts(mfilename, 'OutFolder', figure07Dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_ipsc_comparison (cellName, popIterName, gIncr, ...
                                inFolder, outFolder, ...
                                figTypes, figWidth, figHeight)
% Plot an IPSC comparison plot

% Create a gIncr string
gIncrStr = ['gIncr', num2str(gIncr)];

% Create figure names
figPathBase = fullfile(outFolder, [cellName, '_', popIterName, '_', gIncrStr, ...
                                    '_gabab_ipsc_comparison']);
figPathBaseOrig = [figPathBase, '_orig'];

% Create the figure
fig = set_figure_properties('AlwaysNew', true);

% Plot comparison
m3ha_network_plot_gabab('SaveNewFlag', false, 'InFolder', inFolder, ...
                            'FigTitle', 'suppress', ...
                            'AmpScaleFactor', gIncr);

% Save original figure
save_all_figtypes(fig, figPathBaseOrig, 'png');

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Width', figWidth, 'Height', figHeight, ...
                        'AlignSubplots', true);

% Save the figure
save_all_figtypes(fig, figPathBase, figTypes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_2cell_examples (cellName, popIterName, gIncr, pharm, ...
                                inFolder, outFolder, ...
                                figTypes, figWidth, figHeight)
% Plot 2-cell network examples

% Create a gIncr string
gIncrStr = ['gIncr', num2str(gIncr)];
pharmStr = ['pharm', num2str(pharm)];

% Create figure names
figPathBase = fullfile(outFolder, [cellName, '_', popIterName, '_', gIncrStr, ...
                                    '_', pharmStr, '_example']);
figPathBaseOrig = [figPathBase, '_orig'];

% Create the figure
fig = set_figure_properties('AlwaysNew', true);

% Plot example
m3ha_network_plot_essential('SaveNewFlag', false, 'InFolder', inFolder, ...
                            'FigTitle', 'suppress', ...
                            'AmpScaleFactor', gIncr, 'PharmCondition', pharm);

% Save original figure
save_all_figtypes(fig, figPathBaseOrig, 'png');

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Width', figWidth, 'Height', figHeight, ...
                        'AlignSubplots', true);

% Save the figure
save_all_figtypes(fig, figPathBase, figTypes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combine_2cell_data (popIterDir, candCellSheetPath, ...
                            rankNumsToUse, popDataPath);

%% Hard-coded parameters
rankNumStr = 'rankNum';
cellNameStr = 'cellName';
oscParamsSuffix = 'oscillation_params';

%% Do the job
% Read the candidate cell table
candCellTable = readtable(candCellSheetPath, 'ReadRowNames', true);

% Find the cell names to use from the table
% TODO: table_lookup.m
% TODO: cellNamesToUse = table_lookup(candCellTable, cellNameStr, rankNumStr, rankNumsToUse)
rankNumbersAll = candCellTable.(rankNumStr);
cellNamesAll = candCellTable.(cellNameStr);
cellNamesToUse = match_positions(cellNamesAll, rankNumbersAll, rankNumsToUse);

% Locate corresponding oscillation parameter paths
[~, oscParamPaths] = ...
    find_matching_files(cellNamesToUse, 'Directory', popIterDir, ...
                        'Suffix', oscParamsSuffix, 'Extension', 'csv', ...
                        'Recursive', true, 'ForceCellOutput', true);

% Read the oscillation parameter tables
oscParamTables = cellfun(@readtable, oscParamPaths, 'UniformOutput', false);

% Add the cell name to the tables
oscParamTables = ...
    cellfun(@(x, y) addvars_custom(x, {y}, 'NewVariableNames', cellNameStr, ...
                                    'Before', 1), ...
            oscParamTables, cellNamesToUse, 'UniformOutput', false);

% Vertically concatenate the tables
oscPopTable = apply_over_cells(@vertcat, oscParamTables);

% Join the candidate cell info to the table
oscPopTable = join(oscPopTable, candCellTable, 'Keys', cellNameStr);

% Save the table
writetable(oscPopTable, popDataPath);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function myTable = addvars_custom (myTable, value, varargin)
%% Adds a column to a table, matching rows if necessary 
%% TODO: Pull out as its own function

% Count the number of rows
nRows = height(myTable);

% Make sure value has the same number of rows
value = match_row_count(value, nRows);

% Add the column to the table
myTable = addvars(myTable, value, varargin{:});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function statsTable = m3ha_network_compute_statistics (popDataPath, gIncr, ...
                                                    measureStr, measureTitle)
%% Computes all statistics for the 2-cell network

%% Hard-coded parameters
cellNameStr = 'cellName';
gIncrStr = 'gIncr';
pharmStr = 'pCond';
dclamp2NetworkAmpRatio = 12;

%% Do the job
% Read the data table
popDataTable = readtable(popDataPath);

% Restrict to the rows with given gIncr
rowsToUse = round(popDataTable.(gIncrStr) * dclamp2NetworkAmpRatio) == gIncr;

% Locate the columns of interest
colsOfInterest = [{cellNameStr}; {pharmStr}; measureStr];

% Extract the table of interest
popTableOfInterest = popDataTable(rowsToUse, colsOfInterest);

% Compute statistics for each measure of interest
[allValues, pharmCondition] = ...
    cellfun(@(x) m3ha_network_stats_helper(popTableOfInterest, pharmStr, x), ...
                    measureStr, 'UniformOutput', false);

% Create the statistics table
statsTable = table(measureTitle, measureStr, pharmCondition, allValues, ...
                    'RowNames', measureStr);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [allValues, pharmCondition] = ...
                m3ha_network_stats_helper (popDataTable, pharmStr, measureStr)
%% Plots one violin plot for the 2-cell network

%% Do the job
% Get all possible pharmacological conditions
pharmAll = popDataTable.(pharmStr);
pharmCondition = force_column_cell(num2cell(unique(pharmAll, 'sorted')));

% Find corresponding row numbers
rowsEachPharm = cellfun(@(x) pharmAll == x, ...
                        pharmCondition, 'UniformOutput', false);

% Get all values for this measure under each pharm condition
allValues = cellfun(@(x) popDataTable{x, measureStr}, ...
                    rowsEachPharm, 'UniformOutput', false);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

rowToUse = find_in_list(rankNumsToUse, rankNumbersAll);
cellNamesToUse = cellNamesAll(rowToUse);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
