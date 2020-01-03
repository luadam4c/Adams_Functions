% function [output1] = m3ha_rank_neurons (reqarg1, varargin)
%% Ranks neurons according to how well they are fitted
% Usage: [output1] = m3ha_rank_neurons (reqarg1, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
%
% Arguments:
%       reqarg1     - TODO: Description of reqarg1
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/argfun.m
%       cd/create_labels_from_numbers.m
%       cd/create_label_from_sequence.m
%       cd/create_subplots.m
%       cd/create_time_stamp.m
%       cd/check_dir.m
%       cd/compile_mod_files.m
%       cd/convert_to_char.m
%       cd/copy_into.m
%       cd/create_error_for_nargin.m
%       cd/force_matrix.m
%       cd/m3ha_decide_on_plot_vars.m
%       cd/m3ha_locate_homedir.m
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_select_cells.m
%       cd/m3ha_select_raw_traces.m
%       cd/plot_bar.m
%       cd/plot_histogram.m
%       cd/plot_table_parallel.m
%       cd/save_all_figtypes.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-12-03 Created by Adam Lu
% 2019-12-18 Added bar plots and error history
% 2019-12-30 Added error param comparison across cells
% 

%% Hard-coded parameters
% Flags
chooseBestParamsFlag = false; %true;
plotErrorHistoryFlag = true;
plotErrorComparisonFlag = true;
plotParamHistoryFlag = true;
plotIndividualFlag = true;
rankNeuronsFlag = false; %true;
plotHistogramsFlag = false; %true;
plotBarPlotFlag = false; %true;
plotParamViolinsFlag = false; %true;
plotErrorParamComparisonFlag = true;

% Fitting parameters 
%   Note: Must be consistent with singleneuronfitting91.m
useHH = false;
buildMode = 'active';
simMode = 'active';
dataMode = 2; %0; %3;                   % data mode:
                                    %   0 - all data
                                    %   1 - all of g incr = 100%, 200%, 400% 
                                    %   2 - same g incr but exclude 
                                    %       cell-pharm-g_incr sets 
                                    %       containing problematic sweeps
                                    %   3 - all data but exclude 
                                    %       cell-pharm-g_incr sets 
                                    %       containing problematic sweeps
columnMode = 1;                     % optimization mode:
                                    %   1 - across trials
                                    %   2 - across cells TODO
rowModeAcrossTrials = 1;            % row mode when fitting across trials:
                                    %   1 - each row is a pharm condition
                                    %   2 - each row is a pharm, g incr pair
attemptNumberAcrossTrials = 4;      % attempt number for across trials:
                                    %   1 - Use 4 traces @ 200% gIncr 
                                    %           for this data mode
                                    %   2 - Use all traces @ 200% gIncr 
                                    %           for this data mode
                                    %   3 - Use all traces for this data mode
                                    %   4 - Use 1 trace for each pharm x gIncr 
                                    %           for this data mode
                                    %   5 - Use 4 traces @ 400% gIncr 
                                    %       for this data mode

% Directory names
parentDirectoryTemp = '/media/adamX/m3ha';
fitDirName = 'optimizer4gabab';
dataDirName = fullfile('data_dclamp', 'take4');
matFilesDirName = 'matfiles';
specialCasesDirName = 'special_cases';
defaultOutFolderStr = 'ranked';

% File info
%   Note: Must be consistent with m3ha_neuron_choose_best_params.m
errorSheetSuffix = 'error_param_table';
errorSheetExtension = 'csv';

rankSuffix = 'ranked';
rankSheetExtension = 'csv';

% Default parameters used in computing errors
%   Note: Should be consistent with singleneuronfitting78.m
%       & compute_lts_errors.m & compute_single_neuron_errors.m
ltsFeatureStrings = {'peak amp', 'peak time', 'max slope value'};
sweepWeights = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];
%sweepWeights = [];
errorWeights = [1; 3; 1; 1; 1];
ltsFeatureWeights = errorWeights(3:5);  
                                % weights for LTS feature errors
missedLtsError = 1.5;           % how much error (dimensionless) to 
                                %   penalize a sweep that failed to match 
                                %   a recorded LTS existence
falseLtsError = 0.5;            % how much error (dimensionless) to 
                                %   penalize a sweep that produced an LTS 
                                %   that is not recorded
match2FeatureErrorRatio = errorWeights(2) / sum(ltsFeatureWeights);
                                % ratio of LTS match error to LTS feature error
lts2SweepErrorRatio = sum(errorWeights(2:5)) / errorWeights(1);
                                % ratio of LTS error to sweep error
normalize2InitErrFlag = 0;      % whether to normalize errors to initial values

%   Note: The following must be consistent with compute_single_neuron_errors.m
idxSweep = 1;
idxMatch = 2;
idxAmp = 3;
idxTime = 4;
idxSlope = 5;
totalErrorStr = 'totalError';
avgSwpErrorStr = 'avgSwpError';
ltsMatchErrorStr = 'ltsMatchError';
avgLtsAmpErrorStr = 'avgLtsAmpError';
avgLtsDelayErrorStr = 'avgLtsDelayError';
avgLtsSlopeErrorStr = 'avgLtsSlopeError';
errorWeightsStr = 'errorWeights';
cellNameStr = 'cellName';

% TODO: Make optional argument
% outFolder = '20191227_ranked_singleneuronfitting0-90';
% rankNumsToPlot = [1, 2, 5, 6, 8, 9, 10, 11, 23, 34];
% iterSetStr = 'singleneuronfitting0-90';

% outFolder = '20191229_ranked_singleneuronfitting0-91';
% rankNumsToPlot = [1, 2, 5, 7, 8, 9, 10, 13, 17, 34];

outFolder = '';
figTypes = {'png', 'epsc2'};
rankNumsToPlot = 1:10;
errorParamXTicks = 6:6:36;
iterSetStr = 'singleneuronfitting0-94';
paramDirNames = fullfile('best_params', ...
                        {'bestparams_20191112_singleneuronfitting0', ...
                        'bestparams_20191112_singleneuronfitting1', ...
                        'bestparams_20191112_singleneuronfitting57', ...
                        'bestparams_20191120_singleneuronfitting60', ...
                        'bestparams_20191122_singleneuronfitting61', ...
                        'bestparams_20191123_singleneuronfitting62', ...
                        'bestparams_20191125_singleneuronfitting63', ...
                        'bestparams_20191129_singleneuronfitting72', ...
                        'bestparams_20191201_singleneuronfitting73', ...
                        'bestparams_20191203_singleneuronfitting74', ...
                        'bestparams_20191205_singleneuronfitting75', ...
                        'bestparams_20191211_singleneuronfitting76', ...
                        'bestparams_20191218_singleneuronfitting78', ...
                        'bestparams_20191219_singleneuronfitting85', ...
                        'bestparams_20191221_singleneuronfitting86', ...
                        'bestparams_20191225_singleneuronfitting90', ...
                        'bestparams_20191227_singleneuronfitting91', ...
                        'bestparams_20191230_singleneuronfitting92', ...
                        'bestparams_20191231_singleneuronfitting93', ...
                        'bestparams_20191231_singleneuronfitting94'});
outFolder = '20200102_ranked_singleneuronfitting0-94';
% rankNumsToPlot = [1:6, 21, 36];
% rankNumsToPlot = [1:6, 8, 10, 11, 13:18, 21, 22, 24, 27, 28, 36];
rankNumsToPlot = [1:6, 8:12, 13:18, 21, 22, 24, 26:32, 34, 36];

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% Deal with arguments
% Check number of required arguments
if nargin < 1    % TODO: 1 might need to be changed
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'reqarg1');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, reqarg1, varargin{:});
param1 = iP.Results.param1;

% Keep unmatched arguments for the TODO() function
otherArguments = iP.Unmatched;
%}

%% Preparation
% Locate the home directory
% parentDirectory = m3ha_locate_homedir;
parentDirectory = parentDirectoryTemp;

% Locate the fit directory
fitDirectory = fullfile(parentDirectory, fitDirName);

% Decide on output folder
if isempty(outFolder)
    % Create output folder name
    outFolderName = strcat(create_time_stamp('FormatOut', 'yyyymmdd'), ...
                            '_', defaultOutFolderStr, '_', iterSetStr);

    % Create full path to output folder
    outFolder = fullfile(fitDirectory, outFolderName);
end

% Check if output folder exists
check_dir(outFolder);

% Create a path for the combined spreadsheet
rankBase = [iterSetStr, '_', rankSuffix];
rankPathBase = fullfile(outFolder, rankBase);
rankSheetPath = [rankPathBase, '.', rankSheetExtension];

paramViolinPathBase = [rankPathBase, '_params_violin'];

%% Choose the best parameters for each cell
if chooseBestParamsFlag
    %% Preparation
    % Locate the data directory
    dataDir = fullfile(parentDirectory, dataDirName);

    % Construct full paths to other directories used 
    %   and previously analyzed results under dataDir
    [matFilesDir, specialCasesDir] = ...
        argfun(@(x) fullfile(dataDir, x), matFilesDirName, specialCasesDirName);

    % Decide on candidate parameters directory(ies)
    paramDirs = fullfile(fitDirectory, paramDirNames);

    %% Backup and compile
    % Display message
    fprintf('Backing up parameters ... \n');

    % Copy the candidate params directory(ies) over for backup
    copy_into(paramDirs, outFolder);

    % Display message
    fprintf('Compiling all .mod files ... \n');

    % Compile or re-compile .mod files in the fitting directory
    compile_mod_files(fitDirectory);

    %% Select recorded data
    % Display message
    fprintf('Selecting sweeps to fit for all cells ... \n');

    % Select cells with sweeps to fit for all pharm-gIncr pairs
    [cellIdsToFit, cellInfo, swpInfo] = ...
        m3ha_select_cells('DataMode', dataMode, 'CasesDir', specialCasesDir);

    % Get all cell names to fit
    cellNamesToFit = cellInfo{cellIdsToFit, 'cellName'};

    % Count the number of cells that were fitted 
    nCellsToFit = numel(cellNamesToFit);

    % Select the raw traces to import for each cell to fit
    [fileNamesToFit, rowConditionsToFit] = ...
        arrayfun(@(x) m3ha_select_raw_traces(x, 'ColumnMode', columnMode, ...
                    'RowMode', rowModeAcrossTrials, ...
                    'AttemptNumber', attemptNumberAcrossTrials, ...
                    'SwpInfo', swpInfo, 'CellInfo', cellInfo), ...
                cellIdsToFit, 'UniformOutput', false);

    %% Find the best parameters for each cell
    % Find all the possible initial parameters files for each cell to fit
    [~, customInitPathsToFit] = ...
        cellfun(@(x) all_files('Directory', paramDirs, 'Keyword', x), ...
                                cellNamesToFit, 'UniformOutput', false);

    % Find the best parameters for each cell
     for iCellToFit = nCellsToFit:-1:1
        % Extract stuff for this cell
        cellNameThis = cellNamesToFit{iCellToFit};
        fileNamesThis = fileNamesToFit{iCellToFit};
        rowConditionsThis = rowConditionsToFit{iCellToFit};
        customInitPathsThis = customInitPathsToFit{iCellToFit};

        % Display message
        fprintf('Choosing initial parameters for cell %s ... \n', cellNameThis);

        % Choose the best initial parameters for each cell among all the
        %   custom files
        [previousBestParamsTable, chosenTableLabel] = ...
            m3ha_neuron_choose_best_params(customInitPathsThis, ...
                'PlotErrorHistoryFlag', plotErrorHistoryFlag, ...
                'PlotErrorComparisonFlag', plotErrorComparisonFlag, ...
                'PlotParamHistoryFlag', plotParamHistoryFlag, ...
                'PlotIndividualFlag', plotIndividualFlag, ...
                'OutFolder', outFolder, 'FileNames', fileNamesThis, ...
                'BuildMode', buildMode, 'SimMode', simMode, ...
                'UseHH', useHH, 'RowConditionsIpscr', rowConditionsThis, ...
                'SweepWeightsIpscr', sweepWeights, ...
                'LtsFeatureWeights', ltsFeatureWeights, ...
                'MissedLtsError', missedLtsError, ...
                'FalseLtsError', falseLtsError, ...
                'Lts2SweepErrorRatio', lts2SweepErrorRatio, ...
                'Match2FeatureErrorRatio', match2FeatureErrorRatio, ...
                'Normalize2InitErrFlag', normalize2InitErrFlag, ...
                'SaveParamsFlag', false, 'SaveSimCmdsFlag', false, ...
                'SaveSimOutFlag', false, 'SaveStdOutFlag', false, ...
                'SaveLtsInfoFlag', false, 'SaveLtsStatsFlag', false, ...
                'PlotConductanceFlag', false, 'PlotCurrentFlag', false, ...
                'PlotResidualsFlag', false, 'PlotOverlappedFlag', false, ...
                'PlotIpeakFlag', false, 'PlotLtsFlag', false, ...
                'PlotStatisticsFlag', false, 'PlotSwpWeightsFlag', false);
    end
end

%% Combine the errors and rank neurons
if rankNeuronsFlag
    % Display message
    fprintf('Ranking all cells ... \n');

    % Locate all individual error spreadsheets
    [~, errorSheetPaths] = all_files('Directory', outFolder, ...
                                'Suffix', errorSheetSuffix, ...
                                'Extension', errorSheetExtension);

    % Read all error tables
    errorTables = cellfun(@readtable, errorSheetPaths, 'UniformOutput', false);

    % Restrict to the row with the minimal error
    minimalErrorTables = ...
        cellfun(@(x, y) process_table_for_ranking(x, totalErrorStr), ...
                errorTables, 'UniformOutput', false);

    % Combine into a single table
    minimalErrorTableCombined = vertcat(minimalErrorTables{:});

    % Sort by total error, with minimal error first
    rankTable = sortrows(minimalErrorTableCombined, 'totalError', 'ascend');

    % Set ranking
    ranking = transpose(1:height(rankTable));

    % Add a column for ranking
    rankTable = addvars(rankTable, ranking, 'Before', 1);

    % Save the rank table
    writetable(rankTable, rankSheetPath);

    % Copy and rename .png files according to ranking
    copy_and_rename_png_files(rankTable, outFolder, outFolder);
end

%% Plot histograms
if plotHistogramsFlag
    % Display message
    fprintf('Plotting histograms ... \n');

    % Read the rank table
    rankTable = readtable(rankSheetPath);

    % Plot histograms for error
    % TODO: Pass in 'NBins', 10, instead
    figure;
    plot_histogram(rankTable.(totalErrorStr), 'Edges', 0:0.5:6, ...
                    'XLabel', 'Total Error', 'YLabel', 'Cell Count', ...
                    'FigTitle', 'Total Error Distribution');
    figure;
    plot_histogram(rankTable.(avgSwpErrorStr), 'Edges', 0:1:10, ...
                    'XLabel', 'Average Sweep Error', 'YLabel', 'Cell Count', ...
                    'FigTitle', 'Average Sweep Error Distribution');
    figure;
    plot_histogram(rankTable.(ltsMatchErrorStr), 'Edges', 0:0.5:6, ...
                    'XLabel', 'LTS Mismatch Error', 'YLabel', 'Cell Count', ...
                    'FigTitle', 'LTS Mismatch Error Distribution');
end

%% Plot a stacked horizontal bar plot comparing errors
if plotBarPlotFlag
    % Display message
    fprintf('Plotting bar plot ... \n');

    % Read the rank table
    rankTable = readtable(rankSheetPath);

    % Count the number of cells (number of rows)
    nCells = height(rankTable);

    % Extract fields
    [totalErrors, avgSwpErrors, ltsMatchErrors, ...
            avgLtsAmpErrors, avgLtsDelayErrors, avgLtsSlopeErrors, ...
            cellNames] = ...
        argfun(@(x) rankTable.(x), ...
                totalErrorStr, avgSwpErrorStr, ltsMatchErrorStr, ...
                avgLtsAmpErrorStr, avgLtsDelayErrorStr, avgLtsSlopeErrorStr, ...
                cellNameStr);

    % Create error weight labels
    errorWeightsStrs = create_labels_from_numbers(1:5, ...
                                            'Prefix', [errorWeightsStr, '_']);

    % Extract errorWeights
    errorWeights = cellfun(@(x) rankTable.(x), errorWeightsStrs, ...
                            'UniformOutput', false);
    errorWeights = force_matrix(errorWeights);

    % Make sure error weights are normalized
    errorWeights = errorWeights ./ repmat(sum(errorWeights, 2), 1, ...
                                            size(errorWeights, 2));

    % Compute components of total error
    %   Note: must match groupLabels
    componentErrors = [ltsMatchErrors, avgSwpErrors, avgLtsAmpErrors, ...
                        avgLtsDelayErrors, avgLtsSlopeErrors] .* ...
            errorWeights(:, [idxMatch, idxSweep, idxAmp, idxTime, idxSlope]);
   
    % Decide on group labels
    %   Note: must match componentErrors
    groupLabels = {'LTS Match Error', 'Sweep Error', 'LTS Amp Error', ...
                    'LTS Time Error', 'LTS Slope Error'};

    % Decide on ticks and tick labels
    pTicks = transpose(1:nCells);
    pTickLabels = cellNames;

    % Create a figure with two subplots
    [fig, ax] = create_subplots(1, 2, 'AlwaysNew', true, ...
                                'FigExpansion', [2.4, 1.2]);

    % Plot components of total error stacked
    subplot(ax(1));
    plot_bar(componentErrors, 'BarDirection', 'horizontal', ...
            'ReverseOrder', true, 'GroupStyle', 'stacked', ...
            'PLabel', 'suppress', 'ReadoutLabel', 'Error (dimensionless)', ...
            'PTicks', pTicks, 'PTickLabels', pTickLabels, ...
            'ColumnLabels', groupLabels, ...
            'FigTitle', 'Total error separated by components');

    % Plot total error for verification
    subplot(ax(2));
    plot_bar(totalErrors, 'BarDirection', 'horizontal', ...
            'ReverseOrder', true, ...
            'PLabel', 'suppress', 'ReadoutLabel', 'Error (dimensionless)', ...
            'PTicks', pTicks, ...
            'ColumnLabels', groupLabels, ...
            'FigTitle', 'Total error');

    % Save figure
    save_all_figtypes(fig, rankPathBase, figTypes);
end

%% Plot a parameters as violin plots
if plotParamViolinsFlag
    % Display message
    fprintf('Plotting parameters as violin plots ... \n');

    % Plot violin plot for geometric parameters

    % Plot violin plot for conductances
    % Create figure
    fig = set_figure_properties('AlwaysNew', true);

    % Plot violin plot
    % TODO: Same plot? Standardize y axis?
    % violins = plot_violin(paramValues, 'XTickLabels', paramNames, ...
    %                         'YLabel', );

    % Save figure
    save_all_figtypes(fig, paramViolinPathBase, figTypes);
end

%% Plot an error and parameter comparison plot
if plotErrorParamComparisonFlag
    % Display message
    fprintf('Plotting error parameter comparison ... \n');

    % Read the rank table
    rankTable = readtable(rankSheetPath);

    % Decide on the errors and parameters to plot
    [errorParamToPlot, errorParamLabels, ...
            errorParamYLimits, errorParamIsLog] = m3ha_decide_on_plot_vars;

    % Create figure title and file name
    errorParamFigTitle = ['Error & Parameter Comparison for all cells'];
    errorParamFigName = strcat(rankPathBase, '_param_comparison_all');

    % Plot error & parameter comparison
    plot_table_parallel(rankTable, 'VarsToPlot', errorParamToPlot, ...
            'RowsToPlot', 'all', 'VarIsLog', errorParamIsLog, ...
            'YLimits', errorParamYLimits, 'XTicks', errorParamXTicks, ...
            'XLabel', 'Rank Number', 'YLabel', errorParamLabels, ...
            'FigTitle', errorParamFigTitle, 'FigTypes', figTypes, ...
            'FigName', errorParamFigName);

    % Create a rank string
    rankStr = ['rank', create_label_from_sequence(rankNumsToPlot)];

    % Create figure title and file name
    selectedFigTitle = ['Error & Parameter Comparison for cells with ranks ', rankStr];
    selectedFigName = strcat(rankPathBase, '_param_comparison_', rankStr);

    % Decide on x ticks
    selectedXTicks = 1:numel(rankNumsToPlot);

    % Plot error & parameter comparison
    plot_table_parallel(rankTable, 'VarsToPlot', errorParamToPlot, ...
            'RowsToPlot', rankNumsToPlot, 'VarIsLog', errorParamIsLog, ...
            'YLimits', errorParamYLimits, 'XTicks', selectedXTicks, ...
            'XLabel', 'Rank Number', 'YLabel', errorParamLabels, ...
            'FigTitle', selectedFigTitle, 'FigTypes', figTypes, ...
            'FigName', selectedFigName);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errorTable = process_table_for_ranking(errorTable, totalErrorStr)
%% Restricts table to row of minimal error

% Extract the total errors
totalError = errorTable.(totalErrorStr);

% Extract the row number with the minimal error
[~, rowMinimalError] = min(totalError);

% Restrict table to that row
errorTable = errorTable(rowMinimalError, :);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function copy_and_rename_png_files (rankTable, inFolder, outFolder)
%% Copy and rename png files associated with candidate labels

% Extract candidate labels
candLabels = rankTable.candLabel;

% Extract rankings
ranking = rankTable.ranking;

% Convert rankings to strings
% rankingStrs = convert_to_char(ranking, 'FormatSpec', '%2.f');
rankingStrs = convert_to_char(ranking, 'FormatSpec', '%d');

% Construct old paths
oldPngPaths = fullfile(inFolder, strcat(candLabels, '_individual.png'));

% Construct new paths
newPngPaths = fullfile(outFolder, strcat('rank_', rankingStrs, '_', ...
                                        candLabels, '_individual.png'));

% Copy files 
cellfun(@(x, y) copyfile(x, y), oldPngPaths, newPngPaths, ...
        'UniformOutput', false);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Force as a cell array
cellName = force_column_cell(cellName);

% Add a column for cell name
errorTable = addvars(errorTable, cellName, 'Before', 1);

% All cell names should be the same
if numel(unique(cellNames)) ~= 1
    error('Cell name not all the same in %s!', errorSheetPaths{iTable});
else
    cellName = cellNames{1};
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
