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
%       cd/apply_over_cells.m
%       cd/archive_dependent_scripts.m
%       cd/argfun.m
%       cd/create_label_from_sequence.m
%       cd/create_subplots.m
%       cd/create_time_stamp.m
%       cd/check_dir.m
%       cd/compile_mod_files.m
%       cd/convert_to_char.m
%       cd/copy_into.m
%       cd/create_error_for_nargin.m
%       cd/isemptycell.m
%       cd/m3ha_decide_on_plot_vars.m
%       cd/m3ha_decide_on_sweep_weights.m
%       cd/m3ha_extract_component_errors.m
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
% 2020-01-03 Restricted fitting window
% 2020-01-24 Now uses m3ha_decide_on_sweep_weights.m
% 2020-04-28 Added useCvode
% 2020-04-28 Added updateScriptsFlag
% 2019-04-28 Changed timeToStabilize from 2000 to 3000

%% Hard-coded parameters
% Flags
updateScriptsFlag = false; %true;
chooseBestParamsFlag = false; %true;
plotIndividualFlag = false; %true;
rankNeuronsFlag = false; %true;
plotHistogramsFlag = false; %true;
plotBarPlotFlag = false; %true;
plotErrorParamComparisonAllFlag = false; %true;
plotErrorParamComparisonSelectedFlag = true;
archiveScriptsFlag = true;

plotErrorHistoryFlag = false; %true;
plotErrorComparisonFlag = false; %true;
plotParamHistoryFlag = false; %true;

% Not implemented
plotParamViolinsFlag = false; %true;

% Fitting parameters 
%   Note: Must be consistent with singleneuronfitting91.m
useHH = false;
buildMode = 'active';
simMode = 'active';
columnMode = 1;                     % optimization mode:
                                    %   1 - across trials
                                    %   2 - across cells TODO
rowModeAcrossTrials = 1;            % row mode when fitting across trials:
                                    %   1 - each row is a pharm condition
                                    %   2 - each row is a pharm, g incr pair

% Directory names
parentDirectory = '/media/adamX/m3ha';
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
sweepWeights = [];              % uses m3ha_decide_on_sweep_weights.m
errorWeights = [1; 6; 5; 1; 1];
ltsFeatureWeights = errorWeights(3:5);  
                                % weights for LTS feature errors
missedLtsError = 18;            % how much error (dimensionless) to 
                                %   penalize a sweep that failed to match 
                                %   a recorded LTS existence
falseLtsError = 6;              % how much error (dimensionless) to 
                                %   penalize a sweep that produced an LTS 
                                %   that is not recorded
match2FeatureErrorRatio = errorWeights(2) / sum(ltsFeatureWeights);
                                % ratio of LTS match error to LTS feature error
lts2SweepErrorRatio = sum(errorWeights(2:5)) / errorWeights(1);
                                % ratio of LTS error to sweep error
normalize2InitErrFlag = 0;      % whether to normalize errors to initial values

% The following must be consistent with singleneuron4compgabab.hoc
timeToStabilize = 3000;         % padded time (ms) to make sure initial value 
                                %   of simulations are stabilized

% Fitting window
ipscrWindow = timeToStabilize + [0, 2800];     % only simulate up to that time
fitWindowIpscr = timeToStabilize + [1000, 2800];  
                                % the time window (ms) where all 
                                %   recorded LTS would lie

%   Note: The following must be consistent with compute_single_neuron_errors.m
totalErrorStr = 'totalError';
avgSwpErrorStr = 'avgSwpError';
ltsMatchErrorStr = 'ltsMatchError';
avgLtsAmpErrorStr = 'avgLtsAmpError';
avgLtsDelayErrorStr = 'avgLtsDelayError';
avgLtsSlopeErrorStr = 'avgLtsSlopeError';
cellNameStr = 'cellName';

% TODO: Make optional argument
% outFolderName = '20191227_ranked_singleneuronfitting0-90';
% iterSetStr = 'singleneuronfitting0-90';
% rankNumsToPlot = [1, 2, 5, 6, 8, 9, 10, 11, 23, 34];
% iterSetStr = 'singleneuronfitting0-90';
% dataMode = 2;
% attemptNumberAcrossTrials = 4;
% errorWeights = [1; 3; 1; 1; 1];
% sweepWeights = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];

% outFolderName = '20191229_ranked_singleneuronfitting0-91';
% iterSetStr = 'singleneuronfitting0-91';
% rankNumsToPlot = [1, 2, 5, 7, 8, 9, 10, 13, 17, 34];
% dataMode = 2;
% attemptNumberAcrossTrials = 4;
% errorWeights = [1; 3; 1; 1; 1];
% sweepWeights = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];

% outFolderName = '20200102_ranked_singleneuronfitting0-94';
% iterSetStr = 'singleneuronfitting0-94';
% rankNumsToPlot = [1:6, 21, 36];
% rankNumsToPlot = [1:6, 8, 10, 11, 13:18, 21, 22, 24, 27, 28, 36];
% rankNumsToPlot = [1:6, 8:12, 13:18, 21, 22, 24, 26:32, 34, 36];
% dataMode = 2;
% attemptNumberAcrossTrials = 4;
% errorWeights = [1; 3; 1; 1; 1];
% sweepWeights = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];

% outFolderName = '20200103_ranked_singleneuronfitting0-94';
% iterSetStr = 'singleneuronfitting0-94';
% rankNumsToPlot = 1:11;
% dataMode = 2;
% attemptNumberAcrossTrials = 4;
% errorWeights = [1; 6; 5; 1; 1];
% sweepWeights = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];

% outFolderName = '20200106_ranked_singleneuronfitting0-95';
% iterSetStr = 'singleneuronfitting0-95';
% rankNumsToPlot = 1:11;
% dataMode = 2;
% errorWeights = [1; 6; 5; 1; 1];
% sweepWeights = [];
% attemptNumberAcrossTrials = 3;

% outFolderName = '20200108_ranked_singleneuronfitting0-95';
% iterSetStr = 'singleneuronfitting0-95';
% rankNumsToPlot = 1:11;
% dataMode = 3;
% attemptNumberAcrossTrials = 3;
% errorWeights = [1; 6; 5; 1; 1];
% sweepWeights = [];

% outFolderName = '20200123_ranked_singleneuronfitting0-97';
% iterSetStr = 'singleneuronfitting0-97';
% rankNumsToPlot = 1:11;
% dataMode = 2;
% attemptNumberAcrossTrials = 6;
% errorWeights = [1; 6; 5; 1; 1];
% sweepWeights = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];

% outFolderName = '20200129_ranked_singleneuronfitting101';
% rankNumsToPlot = 1:11;
% rankNumsToPlot = [8, 18];
% rankNumsToPlot = [8, 18, 20, 23, 24, 26, 27, 30, 31, 33, 35, 36];
% rankNumsToPlot = [1, 4, 7, 9, 15, 16];
% rankNumsToPlot = [1, 4, 7, 9, 15, 16, 30];
% rankNumsToPlot = [1, 4, 7, 9, 15, 16, 8, 18, 20, 23, 24, 26, 27, 30, 31, 33, 35, 36];
% iterSetStr = 'singleneuronfitting101';
% paramDirNames = fullfile('best_params', ...
%                         {'bestparams_20200126_singleneuronfitting101'});

% outFolderName = '20200131_ranked_singleneuronfitting0-102';
% iterSetStr = 'singleneuronfitting0-102';
% rankNumsToPlot = 1:11;
% dataMode = 3;
% attemptNumberAcrossTrials = 3;
% sweepWeights = [];              % uses m3ha_decide_on_sweep_weights.m
% errorWeights = [1; 6; 5; 1; 1];

% outFolderName = '20200202_ranked_singleneuronfitting0-102';
% iterSetStr = 'singleneuronfitting0-102';
% rankNumsToPlot = 1:11;
% dataMode = 2;
% attemptNumberAcrossTrials = 3;
% sweepWeights = [];              % uses m3ha_decide_on_sweep_weights.m
% errorWeights = [1; 6; 5; 1; 1];

% outFolderName = '20200203_ranked_manual_singleneuronfitting0-102';
% iterSetStr = 'manual_singleneuronfitting0-102';
% rankNumsToPlot = 1:11;
% dataMode = 2;
% attemptNumberAcrossTrials = 3;
% sweepWeights = [];              % uses m3ha_decide_on_sweep_weights.m
% errorWeights = [1; 6; 5; 1; 1];

% outFolderName = '20200203_ranked_manual_singleneuronfitting0-102';
% iterSetStr = 'manual_singleneuronfitting0-102';
% rankNumsToPlot = [1, 2, 4, 7, 10];            % old best-fitted
% rankNumsToPlot = [1, 2, 5:10, 12:25, 29, 33]; % mistake
% rankNumsToPlot = [1, 2, 4:10, 12:25, 29, 33]; % old well-fitted
% rankNumsToPlot = [2, 7, 10, 12, 20];      % old well-fitted, good response, not from D101310
% rankNumsToPlot = [7, 10, 22, 33];         % old well-fitted, from curve-fitted geometry
% rankNumsToPlot = [1:2, 4:10, 12:25];      % new well-fitted
% dataMode = 2;
% attemptNumberAcrossTrials = 3;
% sweepWeights = [];              % uses m3ha_decide_on_sweep_weights.m
% errorWeights = [1; 6; 5; 1; 1];
% paramDirNames = fullfile('best_params', ...
%                         'bestparams_20200203_manual_singleneuronfitting0-102');
% errorParamXTicks = 6:6:36;
% rankYTickLocs = [1, 8, 15, 22, 29, 36];
% selectedXTicks = [1, 25];

% outFolderName = '20200207_ranked_manual_singleneuronfitting0-102';
% rankNumsToPlot = 1:23;                  % new well-fitted
% rankNumsToPlot = [1:3, 6, 9];               % best-fitted
% iterSetStr = 'manual_singleneuronfitting0-102';
% dataMode = 2;
% attemptNumberAcrossTrials = 3;
% useCvode = true;

% outFolderName = '20200429_ranked_manual_singleneuronfitting0-102';
% rankNumsToPlot = 1:23;
% iterSetStr = 'manual_singleneuronfitting0-102';
% dataMode = 2;
% attemptNumberAcrossTrials = 3;
% useCvode = false;

outFolderName = '20200207_ranked_manual_singleneuronfitting0-102';
rankNumsToPlot = 1:31;
iterSetStr = 'manual_singleneuronfitting0-102';
dataMode = 2;                       % data mode:
                                    %   0 - all data
                                    %   1 - all of g incr = 100%, 200%, 400% 
                                    %   2 - same g incr but exclude 
                                    %       cell-pharm-g_incr sets 
                                    %       containing problematic sweeps
                                    %   3 - all data but exclude 
                                    %       cell-pharm-g_incr sets 
                                    %       containing problematic sweeps
attemptNumberAcrossTrials = 3;      % attempt number for across trials:
                                    %   1 - Use 4 traces @ 200% gIncr 
                                    %           for this data mode
                                    %   2 - Use all traces @ 200% gIncr 
                                    %           for this data mode
                                    %   3 - Use all traces for this data mode
                                    %   4 - Use 1 trace for each pharm x gIncr 
                                    %           for this data mode
                                    %   5 - Use 4 traces @ 400% gIncr 
                                    %       for this data mode
                                    %   6 - Same as 4 but prioritize least vHold
                                    %   7 - Same as 1 but prioritize least vHold
                                    %   8 - Same as 5 but prioritize least vHold
useCvode = true;

figTypes = {'png', 'epsc2'};
% paramDirNames = fullfile('best_params', ...
%                         {'bestparams_20191112_singleneuronfitting0', ...
%                         'bestparams_20191112_singleneuronfitting1', ...
%                         'bestparams_20191112_singleneuronfitting57', ...
%                         'bestparams_20191120_singleneuronfitting60', ...
%                         'bestparams_20191122_singleneuronfitting61', ...
%                         'bestparams_20191123_singleneuronfitting62', ...
%                         'bestparams_20191125_singleneuronfitting63', ...
%                         'bestparams_20191129_singleneuronfitting72', ...
%                         'bestparams_20191201_singleneuronfitting73', ...
%                         'bestparams_20191203_singleneuronfitting74', ...
%                         'bestparams_20191205_singleneuronfitting75', ...
%                         'bestparams_20191211_singleneuronfitting76', ...
%                         'bestparams_20191218_singleneuronfitting78', ...
%                         'bestparams_20191219_singleneuronfitting85', ...
%                         'bestparams_20191221_singleneuronfitting86', ...
%                         'bestparams_20191225_singleneuronfitting90', ...
%                         'bestparams_20191227_singleneuronfitting91', ...
%                         'bestparams_20191230_singleneuronfitting92', ...
%                         'bestparams_20191231_singleneuronfitting94', ...
%                         'bestparams_20200103_singleneuronfitting95', ...
%                         'bestparams_20200120_singleneuronfitting97', ...
%                         'bestparams_20200124_singleneuronfitting99', ...
%                         'bestparams_20200126_singleneuronfitting101', ...
%                         'bestparams_20200129_singleneuronfitting102'});
paramDirNames = fullfile('best_params', ...
                        'bestparams_20200207_manual_singleneuronfitting0-102');

% For plots
rankYTickLocs = 1:4:33;
decisionCutoff = 23.5;
rankFigXLimits = [0, 3.5];
rankFigWidth = 7.5;
rankFigHeight = 8;
selectedFigWidth = 10;
selectedFigHeight = 6;

errorParamXTicks = 1:4:33;
selectedYLimits = {[8, 300]; [5, 400]; [1, 100]; ...
                        [1e-6, 1e-4]; [-80, -60]; ...
                    [1e-9, 1e-1]; [1e-9, 1e-2]; [1e-9, 1e-2]; ...
                        [1e-9, 1e-1]; [1e-9, 1e-2]; ...
                    [1e-9, 1e-1]; [1e-9, 1e-2]; [1e-9, 1e-2]; ...
                        [1e-9, 1e-1]; [1e-9, 1e-2]; ...
                    [1e-9, 1e-1]; [1e-9, 1e-2]; [1e-9, 1e-2]; ...
                        [1e-9, 1e-1]; [1e-9, 1e-2]};
selectedYTicks = {[8, 150, 300]; [5, 200, 400]; [1, 50, 100]; ...
                        [1e-6, 1e-5, 1e-4]; [-70]; ...
                    [1e-9, 1e-1]; [1e-9, 1e-2]; [1e-9, 1e-2]; ...
                        [1e-9, 1e-1]; [1e-9, 1e-2]; ...
                    [1e-9, 1e-1]; [1e-9, 1e-2]; [1e-9, 1e-2]; ...
                        [1e-9, 1e-1]; [1e-9, 1e-2]; ...
                    [1e-9, 1e-1]; [1e-9, 1e-2]; [1e-9, 1e-2]; ...
                        [1e-9, 1e-1]; [1e-9, 1e-2]};

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1
outFolder = '';

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
if isempty(parentDirectory)
    parentDirectory = m3ha_locate_homedir;
end

% Locate the fit directory
fitDirectory = fullfile(parentDirectory, fitDirName);

% Decide on output folder name
if isempty(outFolderName)
    % Create output folder name
    outFolderName = strcat(create_time_stamp('FormatOut', 'yyyymmdd'), ...
                            '_', defaultOutFolderStr, '_', iterSetStr);
end

% Decide on output folder
if isempty(outFolder)
    % Create full path to output folder
    outFolder = fullfile(fitDirectory, outFolderName);
end

% Check if output folder exists
check_dir(outFolder);

% Create a path for the combined spreadsheet
rankBase = [iterSetStr, '_', rankSuffix];
rankPathBase = fullfile(outFolder, rankBase);
rankSheetPath = [rankPathBase, '.', rankSheetExtension];
rankPathBaseOrig = [rankPathBase, '_orig'];

paramViolinPathBase = [rankPathBase, '_params_violin'];
paramViolinPathBaseOrig = [paramViolinPathBase, '_orig'];

%% Make sure NEURON scripts are up to date in outFolder
if updateScriptsFlag
    % Display message
    fprintf('Updating NEURON scripts ... \n');
    update_neuron_scripts(fitDirectory, outFolder);
end

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

    % Change to outFolder
    cd(outFolder);

    %% Backup and compile
    % Display message
    fprintf('Backing up parameters ... \n');

    % Copy the candidate params directory(ies) over for backup
    copy_into(paramDirs, outFolder);

    %% Select recorded data
    % Display message
    fprintf('Selecting sweeps to fit for all cells ... \n');

    % Select cells with sweeps to fit for all pharm-gIncr pairs
    [cellIdsToFit, cellInfo, swpInfo] = ...
        m3ha_select_cells('DataMode', dataMode, 'CasesDir', specialCasesDir);

    % Get all cell names to fit
    cellNamesToFit = cellInfo{cellIdsToFit, 'cellName'};

    % Find all the possible parameters files for each cell to fit
    [~, customInitPathsToFit] = ...
        cellfun(@(x) all_files('Directory', paramDirs, 'Keyword', x), ...
                                cellNamesToFit, 'UniformOutput', false);

    % Remove cells that don't have parameter files
    toRank = ~isemptycell(customInitPathsToFit);
    cellNamesToRank = cellNamesToFit(toRank);
    cellIdsToRank = cellIdsToFit(toRank);
    customInitPathsToRank = customInitPathsToFit(toRank);

    % Select the raw traces to import for each cell to fit
    [fileNamesToRank, rowConditionsToRank] = ...
        arrayfun(@(x) m3ha_select_raw_traces(x, 'ColumnMode', columnMode, ...
                    'RowMode', rowModeAcrossTrials, ...
                    'AttemptNumber', attemptNumberAcrossTrials, ...
                    'SwpInfo', swpInfo, 'CellInfo', cellInfo), ...
                cellIdsToRank, 'UniformOutput', false);

    %% Find the best parameters for each cell
    % Count the number of cells that will be ranked
    nCellsToRank = numel(cellNamesToRank);

    % Find the best parameters for each cell
    for iCellToRank = nCellsToRank:-1:1
        % Extract stuff for this cell
        cellNameThis = cellNamesToRank{iCellToRank};
        fileNamesThis = fileNamesToRank{iCellToRank};
        rowConditionsThis = rowConditionsToRank{iCellToRank};
        customInitPathsThis = customInitPathsToRank{iCellToRank};

        % Display message
        fprintf('Choosing initial parameters for cell %s ... \n', cellNameThis);

        % Set default sweep weights
        sweepWeights = ...
            m3ha_decide_on_sweep_weights(sweepWeights, fileNamesThis);

        % Choose the best initial parameters for each cell among all the
        %   custom files
        m3ha_neuron_choose_best_params(customInitPathsThis, ...
                'PlotErrorHistoryFlag', plotErrorHistoryFlag, ...
                'PlotErrorComparisonFlag', plotErrorComparisonFlag, ...
                'PlotParamHistoryFlag', plotParamHistoryFlag, ...
                'PlotIndividualFlag', plotIndividualFlag, ...
                'OutFolder', outFolder, 'FileNames', fileNamesThis, ...
                'BuildMode', buildMode, 'SimMode', simMode, ...
                'UseHH', useHH, 'UseCvode', useCvode, ...
                'RowConditionsIpscr', rowConditionsThis, ...
                'IpscrWindow', ipscrWindow, ...
                'FitWindowIpscr', fitWindowIpscr, ...
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
    minimalErrorTableCombined = ...
        apply_over_cells(@outerjoin, minimalErrorTables, 'MergeKeys', true);

    % Sort by total error, with minimal error first
    rankTable = sortrows(minimalErrorTableCombined, 'totalError', 'ascend');

    % Set ranking
    ranking = transpose(1:height(rankTable));

    % Add a column for ranking
    rankTable = addvars(rankTable, ranking, 'Before', 1);

    % Save the rank table
    writetable(rankTable, rankSheetPath);

    % Copy and rename .png files according to ranking
    if plotIndividualFlag
        copy_and_rename_png_files(rankTable, outFolder, outFolder);
    end
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
    pLimits = 0.5 + [0, nCells];

    % Extract component errors
    [componentErrors, groupLabels] = m3ha_extract_component_errors(rankTable);

    % Extract cell names and total errors
    [totalErrors, cellNames] = ...
        argfun(@(x) rankTable.(x), totalErrorStr, cellNameStr);

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
            'BarWidth', 1, 'PLimits', pLimits, ...
            'FigTitle', 'Total error separated by components', ...
            'PBoundaries', decisionCutoff);

    % Plot total error for verification
    subplot(ax(2));
    plot_bar(totalErrors, 'BarDirection', 'horizontal', ...
            'ReverseOrder', true, ...
            'PLabel', 'suppress', 'ReadoutLabel', 'Error (dimensionless)', ...
            'PTicks', pTicks, ...
            'ColumnLabels', groupLabels, ...
            'BarWidth', 1, 'PLimits', pLimits, ...
            'FigTitle', 'Total error', ...
            'PBoundaries', decisionCutoff);

    % Save figure
    save_all_figtypes(fig, rankPathBaseOrig, 'png');

    % Create new figure
    fig = set_figure_properties('Units', 'centimeters', ...
                        'Width', rankFigWidth, 'Height', rankFigHeight, ...
                        'AlwaysNew', true);

    plot_bar(componentErrors, 'BarDirection', 'horizontal', ...
            'ReverseOrder', true, 'GroupStyle', 'stacked', ...
            'PLabel', 'Rank', 'ReadoutLabel', 'Error (dimensionless)', ...
            'PTicks', pTicks, 'ColumnLabels', groupLabels, ...
            'BarWidth', 1, 'PLimits', pLimits, ...
            'ReadoutLimits', rankFigXLimits, ...
            'FigTitle', 'Errors of Single Neuron Fits', ...
            'PBoundaries', decisionCutoff);

    % Update figure for CorelDraw
    update_figure_for_corel(fig, 'YTickLocs', rankYTickLocs);

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
if plotErrorParamComparisonAllFlag || plotErrorParamComparisonSelectedFlag
    % Display message
    fprintf('Plotting error parameter comparison ... \n');

    % Read the rank table
    rankTable = readtable(rankSheetPath);

    % Decide on the errors and parameters to plot
    [errorParamToPlot, errorParamLabels, ...
            errorParamYLimits, errorParamIsLog] = m3ha_decide_on_plot_vars;

    % Decide on the errors and parameters to plot for CorelDraw
    [errorParamToPlotForCorel, errorParamIsLogForCorel, ...       
            errorParamYLimitsForCorel, errorParamLabelsForCorel] = ...
        argfun(@(x) x(6:end), errorParamToPlot, errorParamIsLog, ...
                errorParamYLimits, errorParamLabels);

    if plotErrorParamComparisonAllFlag
        % Create figure title and file name
        errorParamFigTitle = ['Error & Parameter Comparison for all cells'];
        errorParamFigName = strcat(rankPathBase, '_param_comparison_all');

        % Plot error & parameter comparison
        plot_table_parallel(rankTable, 'VarsToPlot', errorParamToPlot, ...
                'RowsToPlot', 'all', 'VarIsLog', errorParamIsLog, ...
                'YLimits', errorParamYLimits, 'XTicks', errorParamXTicks, ...
                'SubplotDimensions', [5, 5], ...
                'XLabel', 'Rank Number', 'YLabel', errorParamLabels, ...
                'FigTitle', errorParamFigTitle, 'FigTypes', figTypes, ...
                'FigName', errorParamFigName);
    end

    if plotErrorParamComparisonSelectedFlag
        % Create a rank string
        rankStr = ['rank', create_label_from_sequence(rankNumsToPlot)];

        % Create figure title and file name
        selectedFigTitle = ['Error & Parameter Comparison ', ...
                            'for cells with ranks ', rankStr];
        selectedFigName = strcat(rankPathBase, '_param_comparison_', rankStr);
        selectedFigNameOrig = strcat(selectedFigName, '_orig');

        % Plot error & parameter comparison
        handles = ...
            plot_table_parallel(rankTable, 'VarsToPlot', errorParamToPlot, ...
                'RowsToPlot', rankNumsToPlot, 'VarIsLog', errorParamIsLog, ...
                'YLimits', errorParamYLimits, ...
                'SubplotDimensions', [5, 5], ...
                'XLabel', 'Rank Number', 'YLabel', errorParamLabels, ...
                'FigTitle', selectedFigTitle, 'FigTypes', 'png', ...
                'FigName', selectedFigNameOrig);

        % Create figure
        fig = set_figure_properties('AlwaysNew', true);

        % Plot error & parameter comparison
        handles = ...
            plot_table_parallel(rankTable, ...
                'VarsToPlot', errorParamToPlotForCorel, ...
                'RowsToPlot', rankNumsToPlot, ...
                'VarIsLog', errorParamIsLogForCorel, ...
                'YLimits', selectedYLimits, ...
                'AxTitles', errorParamLabelsForCorel, ...
                'SubplotDimensions', [4, 5], ...
                'XLabel', 'suppress', 'YLabel', 'suppress');
        
        % Change marker sizes
        lines = findobj(fig, 'Type', 'Line');
        set(lines, 'MarkerSize', 2);

        % Decide on x ticks
        if numel(rankNumsToPlot) <= 5
            selectedXTicks = 1:numel(rankNumsToPlot);
        else
            selectedXTicks = linspace(1, numel(rankNumsToPlot), 3);
        end

        % Update figure for CorelDraw
        update_figure_for_corel(fig, 'Units', 'centimeters', ...
            'Width', selectedFigWidth, 'Height', selectedFigHeight, ...
            'XTickLocs', selectedXTicks, 'YTickLocs', selectedYTicks);

        % Save figure
        save_all_figtypes(fig, selectedFigName, figTypes);
    end
end

% Archive all scripts for this run
if archiveScriptsFlag
    archive_dependent_scripts(mfilename, 'OutFolder', outFolder);
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

                        'bestparams_20191231_singleneuronfitting93', ...

% Update ticks
ax = handles.ax;
for i = 6:20
    subplot(ax(i));
    yticks([1e-7, 1e-2]);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
