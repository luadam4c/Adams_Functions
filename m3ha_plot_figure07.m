% m3ha_plot_figure07.m
%% Plots Figure 07 for the GAT Blocker paper
%
% Requires:
%       cd/apply_over_cells.m
%       cd/archive_dependent_scripts.m
%       cd/argfun.m
%       cd/create_label_from_sequence.m
%       cd/decide_on_colormap.m
%       cd/find_matching_files.m
%       cd/force_column_cell.m
%       cd/m3ha_network_plot_gabab.m
%       cd/m3ha_network_plot_essential.m
%       cd/m3ha_plot_violin.m
%       cd/match_positions.m
%       cd/match_row_count.m
%       cd/plot_scale_bar.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/update_figure_for_corel.m

% File History:
% 2020-01-30 Modified from m3ha_plot_figure05.m

%% Hard-coded parameters
% Flags
plotIpscComparison = false; %true;
plot2CellExamples = true;

combine2CellPopulation = false; %true;
plot2CellViolins = false; %true;
plot2CellBars = false; %true;

plot200CellExamples = false; %true;

combine200CellPopulation = false; %true;
plot200CellViolins = false; %true;

archiveScriptsFlag = false; %true;

% Directories
parentDirectory = fullfile('/media', 'adamX', 'm3ha');
figure07Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure07');
figure08Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure08');
networkDirectory = fullfile(parentDirectory, 'network_model');
% exampleIterName2Cell = '20200131T1345_using_bestparams_20200126_singleneuronfitting101';  % 20200131
exampleIterName2Cell = '20200205T1353_using_bestparams_20200203_manual_singleneuronfitting0-102_2cell_examples';
exampleIterName200Cell = '20200204T1239_using_bestparams_20200203_manual_singleneuronfitting0-102_200cell_spikes';
popIterName2Cell = '20200204T1042_using_bestparams_20200203_manual_singleneuronfitting0-102_vtraub_-65_2cell_spikes';
popIterName200Cell = exampleIterName200Cell;
candCellSheetName = 'candidate_cells.csv';
oscParamsSuffix = 'oscillation_params';

% Well-fitted, good 2-cell network response
rankNumsToUse = [2, 4, 5, 7, 9, 10, 12, 13, 16, 20, 21, 23, 25, 29];

% Files

% Analysis settings
% Should be consistent with m3ha_plot_figure03.m & m3ha_plot_figure07.m
exampleCellNames = {'D101310'; 'G101310'};

gIncr = 200;                % Original dynamic clamp gIncr value
pharmConditions = (1:4)';   % Pharmacological conditions
                            %   1 - Control
                            %   2 - GAT 1 Block
                            %   3 - GAT 3 Block
                            %   4 - Dual Block
measuresOfInterest = {'oscDurationSec'; 'oscPeriod2Ms'; ...
                        'oscIndex4'; 'hasOscillation'};
measureTitles = {'Oscillation Duration (sec)'; 'Oscillation Period (ms)'; ...
                    'Oscillatory Index'; 'Has Oscillation'};

% Plot settings
ipscFigWidth = 8.5;
ipscFigHeight = 4;
xLimits2Cell = [2800, 4800];
% yLimits2Cell = {[], [], [], [], [], []};
yLimits2Cell = {[-100, 100], [-100, 100], [-1, 12], [-15, 5], ...
                    [1e-10, 1], [1e-10, 1]};
example2CellFigWidth = 8.5;
example2CellFigHeight = 1.5 * 6;
example200CellFigWidth = 8.5;
example200CellFigHeight = 3;
pharmLabelsShort = {'{\it s}-Con', '{\it s}-GAT1', ...
                    '{\it s}-GAT3', '{\it s}-Dual'};

figTypes = {'png', 'epsc'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Find the directory for this iteration
exampleIterDir2Cell = fullfile(networkDirectory, exampleIterName2Cell);
exampleIterDir200Cell = fullfile(networkDirectory, exampleIterName200Cell);
popIterDir2Cell = fullfile(networkDirectory, popIterName2Cell);
popIterDir200Cell = fullfile(networkDirectory, popIterName200Cell);

% Construct the full path to the candidate cell spreadsheet
candCellSheetPath = fullfile(networkDirectory, candCellSheetName);

% Create a rank string
rankStr = ['rank', create_label_from_sequence(rankNumsToUse)];

% Create a condition label
[conditionLabel2Cell, conditionLabel200Cell] = ...
    argfun(@(x) [x, '_', rankStr, '_gIncr', num2str(gIncr)], ...
            popIterName2Cell, popIterName200Cell);

% Create a population data spreadsheet name
popDataSheetName2Cell = [popIterName2Cell, '_', rankStr, '_', ...
                            oscParamsSuffix, '.csv'];
popDataSheetName200Cell = [popIterName200Cell, '_', rankStr, '_', ...
                            oscParamsSuffix, '.csv'];

% Contruct the full path to the population data spreadsheet
popDataPath2Cell = fullfile(figure07Dir, popDataSheetName2Cell);
popDataPath200Cell = fullfile(figure08Dir, popDataSheetName2Cell);

% Create color maps
colorMapPharm = decide_on_colormap([], 4);
colorMapPharmCell = arrayfun(@(x) colorMapPharm(x, :), ...
                            transpose(1:4), 'UniformOutput', false);

%% Find example files and directories
if plotIpscComparison || plot2CellExamples
    % Find example network directories
    [~, exampleDirs2Cell] = ...
        cellfun(@(x) all_subdirs('Directory', exampleIterDir2Cell, ...
                                'Keyword', x), ...
                exampleCellNames, 'UniformOutput', false);
end
if plot200CellExamples
    % Find example network directories
    [~, exampleDirs200Cell] = ...
        cellfun(@(x) all_subdirs('Directory', exampleIterDir200Cell, ...
                                'Keyword', x), ...
                exampleCellNames, 'UniformOutput', false);
end

%% Plots figures for comparing dynamic clamp ipsc
if plotIpscComparison
    cellfun(@(x, y) plot_ipsc_comparison(x, exampleIterName2Cell, gIncr, y, ...
                                        figure07Dir, figTypes, ...
                                        ipscFigWidth, ipscFigHeight), ...
            exampleCellNames, exampleDirs2Cell);
end

%% Plots example 2-cell networks
if plot2CellExamples
    cellfun(@(a, b) ...
        cellfun(@(x, y) plot_2cell_examples(x, exampleIterName2Cell, ...
                            gIncr, a, y, figure07Dir, figTypes, ...
                            example2CellFigWidth, example2CellFigHeight, ...
                            xLimits2Cell, yLimits2Cell, b), ...
                exampleCellNames, exampleDirs2Cell), ...
        num2cell(pharmConditions), colorMapPharmCell);
end

%% Plots example 200-cell networks
if plot200CellExamples
    arrayfun(@(z) ...
        cellfun(@(x, y) plot_200cell_examples(x, exampleIterName200Cell, ...
                            gIncr, z, y, figure08Dir, figTypes, ...
                        example200CellFigWidth, example200CellFigHeight), ...
                exampleCellNames, exampleDirs200Cell), ...
        pharmConditions);
end

%% Combines quantification over all 2-cell networks
if combine2CellPopulation
    combine_osc_params_data(popIterDir2Cell, candCellSheetPath, ...
                            rankNumsToUse, popDataPath2Cell);
end

%% Combines quantification over all 2-cell networks
if combine200CellPopulation
    combine_osc_params_data(popIterDir200Cell, candCellSheetPath, ...
                            rankNumsToUse, popDataPath200Cell);
end

%% Plots oscillation measures over pharm condition 
%       across all 2-cell networks
if plot2CellViolins
    % Construct stats table path
    stats2dPath2Cell = ...
        fullfile(figure07Dir, strcat(conditionLabel2Cell, '_stats.mat'));

    % Compute statistics if not done already
    if ~isfile(stats2dPath2Cell)
        % Compute statistics for all features
        disp('Computing statistics for violin plots ...');
        statsTable = m3ha_network_compute_statistics(popDataPath2Cell, ...
                                    gIncr, measuresOfInterest, measureTitles);

        % Generate labels
        conditionLabel = conditionLabel2Cell;
        pharmLabels = pharmLabelsShort;

        % Save stats table
        save(stats2dPath2Cell, 'statsTable', 'pharmLabels', ...
                            'conditionLabel', '-v7.3');
    end

    % Plot violin plots
    m3ha_plot_violin(stats2dPath2Cell, 'RowsToPlot', measuresOfInterest, ...
                    'OutFolder', figure07Dir);
end

%% Plots oscillation measures over pharm condition 
%       across all 200-cell networks
if plot200CellViolins
    % Construct stats table path
    stats2dPath200Cell = ...
        fullfile(figure08Dir, strcat(conditionLabel200Cell, '_stats.mat'));

    % Compute statistics if not done already
    if ~isfile(stats2dPath200Cell)
        % Compute statistics for all features
        disp('Computing statistics for violin plots ...');
        statsTable = m3ha_network_compute_statistics(popDataPath200Cell, ...
                                    gIncr, measuresOfInterest, measureTitles);

        % Generate labels
        conditionLabel = conditionLabel200Cell;
        pharmLabels = pharmLabelsShort;

        % Save stats table
        save(stats2dPath200Cell, 'statsTable', 'pharmLabels', ...
                            'conditionLabel', '-v7.3');
    end

    % Plot violin plots
    m3ha_plot_violin(stats2dPath200Cell, 'RowsToPlot', measuresOfInterest, ...
                    'OutFolder', figure08Dir);
end

%% Plots percentage of 2-cell networks having oscillations over pharm condition
if plot2CellBars
end

%% Archive all scripts for this run
if archiveScriptsFlag
    archive_dependent_scripts(mfilename, 'OutFolder', figure07Dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_ipsc_comparison (cellName, popIterName2Cell, gIncr, ...
                                inFolder, outFolder, ...
                                figTypes, figWidth, figHeight)
% Plot an IPSC comparison plot

% Create a gIncr string
gIncrStr = ['gIncr', num2str(gIncr)];

% Create figure names
figPathBase = fullfile(outFolder, [cellName, '_', popIterName2Cell, ...
                        '_', gIncrStr, '_gabab_ipsc_comparison']);
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

function plot_2cell_examples (cellName, iterName, gIncr, pharm, ...
                            inFolder, outFolder, figTypes, ...
                            figWidth, figHeight, xLimits, yLimits, colorMap)
% Plot 2-cell network examples

% Create a gIncr string
gIncrStr = ['gIncr', num2str(gIncr)];
pharmStr = ['pharm', num2str(pharm)];

% Create figure names
figPathBase = fullfile(outFolder, [cellName, '_', iterName, ...
                            '_', gIncrStr, '_', pharmStr, '_2cell_example']);
figPathBaseOrig = [figPathBase, '_orig'];

% Create the figure
fig = set_figure_properties('AlwaysNew', true);

% Plot example
handles = ...
    m3ha_network_plot_essential('SaveNewFlag', false, 'InFolder', inFolder, ...
                        'XLimits', xLimits, 'YLimits', yLimits, ...
                        'FigTitle', 'suppress', ...
                        'AmpScaleFactor', gIncr, 'PharmCondition', pharm, ...
                        'Color', colorMap);

% Save original figure
drawnow;
save_all_figtypes(fig, figPathBaseOrig, 'png');


% Plot a scale bar in the first subplot
subPlots = handles.subPlots;
subplot(subPlots(1));
plot_scale_bar('x', 'XBarUnits', 'ms', 'XBarLength', 200, ...
                'XPosNormalized', 0.9, 'YPosNormalized', 0.9);

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Width', figWidth, 'Height', figHeight, ...
                        'RemoveXRulers', true, 'AlignSubplots', true);

% Save the figure
drawnow;
save_all_figtypes(fig, figPathBase, figTypes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_200cell_examples (cellName, iterName, gIncrDclamp, pharm, ...
                            inFolder, outFolder, figTypes, figWidth, figHeight)
% Plot 200-cell network examples

% Get the gIncr value for the network
gIncr = gIncrDclamp / 12;

% Create strings
gIncrStr = ['gIncr', num2str(gIncrDclamp)];
pharmStr = ['pharm', num2str(pharm)];

% Find the appropriate simulation number
simNumber = m3ha_network_find_sim_number(inFolder, pharm, gIncr);

% Create figure names
figPathBase = fullfile(outFolder, [cellName, '_', iterName, ...
                        '_', gIncrStr, '_', pharmStr, '_200cell_example']);
figPathBaseOrig = [figPathBase, '_orig'];
figPathBaseNoRasters = [figPathBase, '_noRasters'];
figPathBaseRasterOnly = [figPathBase, '_rastersOnly'];

%% Full figure
% Create the figure
fig = set_figure_properties('AlwaysNew', true);

% Plot spike raster plot
m3ha_network_raster_plot(inFolder, 'OutFolder', outFolder, ...
                        'SingleTrialNum', simNumber, ...
                        'PlotSpikes', true, 'PlotTuning', false, ...
                        'PlotOnly', true);

% Save original figure
drawnow;
save_all_figtypes(fig, figPathBaseOrig, 'png');

%% Figure with no rasters
% Create the figure with no rasters
fig = set_figure_properties('AlwaysNew', true, 'Units', 'centimeters', ...
                            'Width', figWidth, 'Height', figHeight);

% Plot spike raster plot again with no rasters
m3ha_network_raster_plot(inFolder, 'OutFolder', outFolder, ...
                        'SingleTrialNum', simNumber, ...
                        'PlotSpikes', true, 'PlotTuning', false, ...
                        'PlotOnly', true, 'NoRasters', true);

% Plot a scale bar only for the Dual Blockade condition
if pharm == 4
    plot_scale_bar('x', 'XBarUnits', 'sec', 'XBarLength', 2, ...
                    'XPosNormalized', 0.6, 'YPosNormalized', 0.2);
end

% Update figure for CorelDraw
update_figure_for_corel(fig, 'RemoveXLabels', true, 'RemoveTitles', true, ...
                        'BoxOn', true, 'RemoveXRulers', true);

% Save figure with no rasters
drawnow;
save_all_figtypes(fig, figPathBaseNoRasters, figTypes);

%% Figure with rasters only
% Create the figure with rasters only
fig = set_figure_properties('AlwaysNew', true, 'Units', 'centimeters', ...
                            'Width', figWidth, 'Height', figHeight);

% Plot spike raster plot again with rasters only
m3ha_network_raster_plot(inFolder, 'OutFolder', outFolder, ...
                        'SingleTrialNum', simNumber, ...
                        'PlotSpikes', true, 'PlotTuning', false, ...
                        'PlotOnly', true, 'RastersOnly', true);

% Update figure for CorelDraw
update_figure_for_corel(fig, 'RemoveXLabels', true, 'RemoveTitles', true, ...
                        'BoxOn', true, 'RemoveXRulers', true);

% Save the figure with rasters only
drawnow;
save_all_figtypes(fig, figPathBaseRasterOnly, figTypes);

% Close all figures
close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function simNumber = m3ha_network_find_sim_number (inFolder, pharm, gIncr)
%% TODO: Move this to m3ha_network_raster_plot.m

% Create strings
paramsPrefix = 'sim_params';
pharmStr = ['pCond_', num2str(pharm)];
gIncrStr = ['gIncr_', num2str(gIncr)];

% Create the keyword
keyword = [pharmStr, '_', gIncrStr];

% Find the sim params file
[~, paramPath] = all_files('Directory', inFolder, 'Prefix', paramsPrefix, ...
                            'Keyword', keyword, 'MaxNum', 1);

% Find the corresponding simNumber
paramTable = readtable(paramPath, 'ReadRowNames', true);
simNumber = paramTable{'simNumber', 'Value'};

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combine_osc_params_data (popIterDir2Cell, candCellSheetPath, ...
                                    rankNumsToUse, popDataPath2Cell);

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
    find_matching_files(cellNamesToUse, 'Directory', popIterDir2Cell, ...
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
writetable(oscPopTable, popDataPath2Cell);

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

function statsTable = m3ha_network_compute_statistics (popDataPath2Cell, ...
                                                gIncr, measureStr, measureTitle)
%% Computes all statistics for the 2-cell network

%% Hard-coded parameters
cellNameStr = 'cellName';
gIncrStr = 'gIncr';
pharmStr = 'pCond';
dclamp2NetworkAmpRatio = 12;

%% Do the job
% Read the data table
popDataTable = readtable(popDataPath2Cell);

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
