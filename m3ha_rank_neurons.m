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
%       cd/check_dir.m
%       cd/compile_mod_files.m
%       cd/convert_to_char.m
%       cd/copy_into.m
%       cd/create_error_for_nargin.m
%       cd/m3ha_locate_homedir.m
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_select_cells.m
%       cd/m3ha_select_raw_traces.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-12-03 Created by Adam Lu
% 

%% Hard-coded parameters
% Flags
chooseBestParamsFlag = true;
plotErrorHistoryFlag = true;
rankNeuronsFlag = true;
plotHistogramsFlag = true;

% Fitting parameters 
%   Note: Must be consistent with singleneuronfitting75.m
simMode = 'active';
dataMode = 0; %3;                   % data mode:
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
paramDirNames = fullfile('best_params', ...
                        {'bestparams_20191120_singleneuronfitting60', ...
                        'bestparams_20191122_singleneuronfitting61', ...
                        'bestparams_20191123_singleneuronfitting62', ...
                        'bestparams_20191125_singleneuronfitting63', ...
                        'bestparams_20191129_singleneuronfitting72', ...
                        'bestparams_20191201_singleneuronfitting73', ...
                        'bestparams_20191203_singleneuronfitting74', ...
                        'bestparams_20191205_singleneuronfitting75', ...
                        'bestparams_20191211_singleneuronfitting76'});
dataDirName = fullfile('data_dclamp', 'take4');
matFilesDirName = 'matfiles';
specialCasesDirName = 'special_cases';
defaultOutFolderName = 'ranked';

% File info
cellNamePattern = '[A-Z][0-9]{6}';

%   Note: Must be consistent with m3ha_neuron_choose_best_params.m
errorSheetSuffix = 'error_comparison';
errorSheetExtension = 'csv';

rankPrefix = 'singleneuronfitting60-76';
rankSheetSuffix = 'ranked';
rankSheetExtension = 'csv';

% Default parameters used in computing errors
%   Note: Should be consistent with singleneuronfitting75.m
%       & compute_lts_errors.m & compute_single_neuron_errors.m
ltsFeatureStrings = {'peak amp', 'peak time', 'max slope value'};
ltsFeatureWeights = [2; 2; 2];  % default weights for optimizing LTS statistics
missedLtsError = 1.5;           % how much error (dimensionless) to 
                                %   penalize a sweep that failed to match 
                                %   a recorded LTS existence
falseLtsError = 0.5;            % how much error (dimensionless) to 
                                %   penalize a sweep that produced an LTS 
                                %   that is not recorded
lts2SweepErrorRatio = 3;        % ratio of LTS error to sweep error
normalize2InitErrFlag = 0;      % whether to normalize errors to initial values
% sweepWeights = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];
sweepWeights = [];

% TODO: Make optional argument
outFolder = '';

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
    outFolder = fullfile(fitDirectory, defaultOutFolderName);
end

% Check if output folder exists
check_dir(outFolder);

% Create a path for the combined spreadsheet
rankName = [rankPrefix, '_', rankSheetSuffix, '.', rankSheetExtension];
rankPath = fullfile(outFolder, rankName);

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
        arrayfun(@(x) m3ha_select_raw_traces(rowModeAcrossTrials, ...
                        columnMode, attemptNumberAcrossTrials, ...
                        x, swpInfo, cellInfo), ...
                cellIdsToFit, 'UniformOutput', false);

    %% Find the best parameters for each cell
    % Find all the possible initial parameters files for each cell to fit
    [~, customInitPathsToFit] = ...
        cellfun(@(x) all_files('Directory', paramDirs, 'Keyword', x), ...
                                cellNamesToFit, 'UniformOutput', false);

    % Find the best parameters for each cell
    for iCellToFit = 1:nCellsToFit
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
                'PlotIndividualFlag', true, ...
                'OutFolder', outFolder, 'FileNames', fileNamesThis, ...
                'SimMode', simMode, 'RowConditionsIpscr', rowConditionsThis, ...
                'SweepWeightsIpscr', sweepWeights, ...
                'LtsFeatureWeights', ltsFeatureWeights, ...
                'MissedLtsError', missedLtsError, ...
                'FalseLtsError', falseLtsError, ...
                'Lts2SweepErrorRatio', lts2SweepErrorRatio, ...
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

%% Plot the error history for each cell
if plotErrorHistoryFlag
end

%% Combine the errors and rank neurons
if rankNeuronsFlag
    % Display message
    fprintf('Ranking all cells ... \n');

    % Locate all individual error spreadsheets
    [~, sheetPaths] = all_files('Directory', outFolder, ...
                                'Suffix', errorSheetSuffix, ...
                                'Extension', errorSheetExtension);

    % Read all error tables
    errorTables = cellfun(@readtable, sheetPaths, 'UniformOutput', false);

    % Restrict to the row with the minimal error
    minimalErrorTables = cellfun(@(x, y) process_table_for_ranking(x), ...
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
    writetable(rankTable, rankPath);

    % Copy and rename .png files according to ranking
    copy_and_rename_png_files(rankTable, outFolder, outFolder);
end

%% Plot histograms
if plotHistogramsFlag
    % Read the rank table
    rankTable = readtable(rankPath);

    % Plot histograms for error
    % TODO: Pass in 'NBins', 10, instead
    figure;
    plot_histogram(rankTable.('totalError'), 'Edges', 0:0.5:6, ...
                    'XLabel', 'Total Error', 'YLabel', 'Cell Count', ...
                    'FigTitle', 'Total Error Distribution');
    figure;
    plot_histogram(rankTable.('ltsMisMatchError'), 'Edges', 0:0.5:6, ...
                    'XLabel', 'LTS Mismatch Error', 'YLabel', 'Cell Count', ...
                    'FigTitle', 'LTS Mismatch Error Distribution');
    figure;
    plot_histogram(rankTable.('avgSwpError'), 'Edges', 0:1:10, ...
                    'XLabel', 'Average Sweep Error', 'YLabel', 'Cell Count', ...
                    'FigTitle', 'Average Sweep Error Distribution');
end

%% Plot a stacked horizontal bar plot comparing errors
if plotBarPlotFlag
    % Read the rank table
    rankTable = readtable(rankPath);

end

%% Output results
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errorTable = process_table_for_ranking(errorTable)
%% Restricts table to row of minimal error

% Extract the total errors
totalError = errorTable.totalError;

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

% Convert rankings to strings so that 1 becomes '01'
rankingStrs = convert_to_char(ranking, 'FormatSpec', '%2.f');

% Construct old paths
oldPngPaths = fullfile(inFolder, strcat(candLabels, '_individual.png'));

% Construct new paths
newPngPaths = fullfile(inFolder, strcat('rank_', rankingStrs, '_', ...
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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%