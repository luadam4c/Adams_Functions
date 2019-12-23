% function [output1] = m3ha_simulate_population (reqarg1, varargin)
%% Generates simulated IPSC responses that can be compared with recorded data
% Usage: [output1] = m3ha_simulate_population (reqarg1, varargin)
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
% TODO
%       cd/all_subdirs.m
%       cd/all_files.m
%       cd/check_dir.m
%       cd/compile_mod_files.m
%       cd/copy_into.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/create_time_stamp.m
%       cd/find_matching_files.m
%       cd/m3ha_extract_cell_name.m
%       cd/m3ha_extract_iteration_string.m
%       cd/m3ha_load_sweep_info.m
%       cd/m3ha_locate_homedir.m
%       cd/m3ha_select_sweeps.m
%       cd/print_cellstr.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-12-11 Created by Adam Lu
% 

%% Hard-coded parameters
% Flags
chooseBestNeuronsFlag = true;
simulateFlag = true;
computeStatsFlag = false;

% Selection parameters
nCellsToSim = 10;

% Simulation parameters 
buildMode = 'active';
simMode = 'active';
dataMode = 0;                       % data mode:
                                    %   0 - all data
                                    %   1 - all of g incr = 100%, 200%, 400% 
                                    %   2 - same g incr but exclude 
                                    %       cell-pharm-g_incr sets 
                                    %       containing problematic sweeps

% Directory names
parentDirectoryTemp = '/media/adamX/m3ha';
fitDirName = 'optimizer4gabab';
rankDirName = 'ranked';
dataDirName = fullfile('data_dclamp', 'take4');
matFilesDirName = 'matfiles';
specialCasesDirName = 'special_cases';
defaultOutFolderSuffix = '_population';

% File info
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
lts2SweepErrorRatio = 3;        % default ratio of LTS error to sweep error
normalize2InitErrFlag = 0;      % whether to normalize errors to initial values
sweepWeights = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];

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

% Locate the ranked directory
rankDirectory = fullfile(fitDirectory, rankDirName);

% Load sweep info
swpInfo = m3ha_load_sweep_info;

% Decide on output folder
if isempty(outFolder)
    % Create output folder name
    oufFolderName = strcat(create_time_stamp('FormatOut', 'yyyymmdd'), ...
                            defaultOutFolderSuffix);

    % Create full path to output folder
    outFolder = fullfile(fitDirectory, oufFolderName);
end

% Check if output folder exists
check_dir(outFolder);

%% Choose the best cells and the best parameters for each cell
if chooseBestNeuronsFlag
    % Create rank number prefixes
    rankPrefixes = create_labels_from_numbers(1:nCellsToSim, ...
                                        'Prefix', 'rank_', 'Suffix', '_');

    % Find png files matching the rank prefixes
    [~, pngPaths] = find_matching_files(rankPrefixes, 'PartType', 'Prefix', ...
                            'Directory', rankDirectory, 'Extension', 'png', ...
                            'ExtractDistinct', false);

    % Extract the cell names
    cellNames = m3ha_extract_cell_name(pngPaths);

    % Extract the iteration numbers
    iterStrs = m3ha_extract_iteration_string(pngPaths);

    % Find the parameter file directories
    [~, paramDirs] = cellfun(@(x) all_subdirs('Directory', rankDirectory, ...
                                        'RegExp', x, 'MaxNum', 1), ...
                        iterStrs, 'UniformOutput', false);

    % Find the parameter files for each cell
    [~, paramPaths] = cellfun(@(x, y) all_files('Directory', x, ...
                            'Keyword', y, 'Suffix', 'params', 'MaxNum', 1), ...
                            paramDirs, cellNames, 'UniformOutput', false);

    % Copy the parameter files into 
    copy_into(paramPaths, outFolder);
end

%% Simulate
if simulateFlag
    %% Decide on the candidate parameter files and cells
    % Decide on candidate parameters files
    [~, paramPaths] = all_files('Directory', outFolder, 'Suffix', 'params');

    % Extract the cell names
    cellNames = m3ha_extract_cell_name(paramPaths);

    % Display message
    fprintf('All sweeps from the following cells will be simulated: \n');
    print_cellstr(cellNames, 'OmitBraces', true, 'Delimiter', '\n');

    %% Compile
    % Display message
    fprintf('Compiling all .mod files ... \n');

    % Compile or re-compile .mod files in the fitting directory
    compile_mod_files(fitDirectory);

    %% Select recorded data
    % Display message
    fprintf('Selecting recorded sweeps for all cells ... \n');

    % Locate the data directory
    dataDir = fullfile(parentDirectory, dataDirName);

    % Construct full paths to other directories used 
    %   and previously analyzed results under dataDir
    [matFilesDir, specialCasesDir] = ...
        argfun(@(x) fullfile(dataDir, x), matFilesDirName, specialCasesDirName);

    % Select the sweep indices that will be simulated
    swpInfo = m3ha_select_sweeps('SwpInfo', swpInfo, 'DataMode', dataMode, ...
                                    'CasesDir', specialCasesDir);

    % Select the raw traces to import for each cell to fit
    [fileNamesToFit, rowConditionsToFit] = ...
        arrayfun(@(x) m3ha_select_raw_traces(rowmodeAcrossTrials, ...
                        columnMode, attemptNumberAcrossTrials, ...
                        x, swpInfo, cellInfo), ...
                cellNamesToFit, 'UniformOutput', false);


end

if false
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
        arrayfun(@(x) m3ha_select_raw_traces(rowmodeAcrossTrials, ...
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
                'BuildMode', buildMode, 'SimMode', simMode, ...
                'RowConditionsIpscr', rowConditionsThis, ...
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

%% Combine the errors and rank neurons
if computeStatsFlag
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

    % Create a path for the combined spreadsheet
    rankName = [rankPrefix, '_', rankSheetSuffix, '.', rankSheetExtension];
    rankPath = fullfile(outFolder, rankName);

    % Save the rank table
    writetable(rankTable, rankPath);

    % Copy and rename .png files according to ranking
    copy_and_rename_png_files(rankTable, outFolder, outFolder);

    % Plot a histogram for total error
    figure;
    plot_histogram(rankTable.('totalError'), 'Edges', 0:0.5:6, ...
                        'XLabel', 'Total Error', 'YLabel', 'Cell Count', ...
                        'FigTitle', 'Total Error Distribution');
    figure;
    plot_histogram(rankTable.('ltsMatchError'), 'Edges', 0:0.5:6, ...
                        'XLabel', 'LTS Mismatch Error', 'YLabel', 'Cell Count', ...
                        'FigTitle', 'LTS Mismatch Error Distribution');
    figure;
    plot_histogram(rankTable.('avgSwpError'), 'Edges', 0:1:10, ...
                        'XLabel', 'Average Sweep Error', 'YLabel', 'Cell Count', ...
                        'FigTitle', 'Average Sweep Error Distribution');

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

%{
OLD CODE:

% Force as a cell array
cellName = force_column_cell(cellName);

% Add a column for cell name
errorTable = addvars(errorTable, cellName, 'Before', 1);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
