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
%       cd/all_files.m
%       cd/all_subdirs.m
%       cd/copy_into.m
%       cd/create_time_stamp.m
%       cd/create_labels_from_numbers.m
%       cd/find_matching_files.m
%       cd/m3ha_compute_statistics.m
%       cd/m3ha_extract_cell_name.m
%       cd/m3ha_extract_iteration_string.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_load_sweep_info.m
%       cd/print_cellstr.m
%       cd/renamevars.m
%       cd/vertcat_spreadsheets.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-12-11 Created by Adam Lu
% 

%% Hard-coded parameters
% Flags
chooseBestNeuronsFlag = false; %true;
simulateFlag = false; %true;
combineFeatureTablesFlag = false; %true;
computeStatsFlag = true;
plotStatsFlag = true;

% Selection parameters
nCellsToSim = 10;

% Simulation parameters 
buildMode = 'active';
simMode = 'active';
dataMode = 0;           % data mode:
                        %   0 - all data
                        %   1 - all of g incr = 100%, 200%, 400% 
                        %   2 - same g incr but exclude 
                        %       cell-pharm-g_incr sets 
                        %       containing problematic sweeps
attemptNumber = 3;      %   1 - Use 4 traces @ 200% gIncr for this data mode
                        %   2 - Use all traces @ 200% gIncr for this data mode
                        %   3 - Use all traces for this data mode
                        %   4 - Use 1 trace for each pharm x gIncr 
                        %           for this data mode
                        %   5 - Use 4 traces @ 400% gIncr for this data mode
% TODO:useHH = true;           % whether to use Hudgin-Huxley Na+ and K+ channels

% Directory names
parentDirectoryTemp = '/media/adamX/m3ha';
fitDirName = 'optimizer4gabab';
rankDirName = 'ranked';
defaultOutFolderSuffix = '_population';

% File names
simStr = 'sim';
ltsParamsSuffix = '_ltsParams';
simLtsParamsSuffix = strcat(simStr, '_ltsParams');
simSwpInfoSuffix = strcat(simStr, '_swpInfo');
stats3dSuffix = strcat(simStr, '_swpInfo');

% Note: The following must be consistent with m3ha_parse_dclamp_data.m
condVarStrs = {'cellidrow', 'prow', 'vrow', 'grow', 'swpnrow', ...
                'gabab_amp', 'gabab_Trise', 'gabab_TfallFast', ...
                'gabab_TfallSlow', 'gabab_w'};
pharmAll = [1; 2; 3; 4];          
pharmLabels = {'{\it d}-Control', '{\it d}-GAT1 Block', ...
                '{\it d}-GAT3 Block', '{\it d}-Dual Block'};
pharmLabelsShort = {'{\it d}-Con', '{\it d}-GAT1', ...
                    '{\it d}-GAT3', '{\it d}-Dual'};
gIncrAll = [25; 50; 100; 200; 400; 800];
gIncrLabels = {'25%', '50%', '100%', '200%', '400%', '800%'};
conditionLabel2D = 'pharm_1-4_gincr_200';
pCond2D = num2cell(pharmAll);
gCond2D = 200;
conditionLabel3D = 'pharm_1-4_gincr_all';
pCond3D = num2cell(pharmAll);
gCond3D = num2cell(gIncrAll);

% TODO: Make optional argument
outFolder = '';

% Plot settings
% Note: must be consistent with m3ha_compute_statistics.m
measuresOfInterest = {'ltsAmplitude'; 'ltsMaxSlope'; ...
                    'ltsConcavity'; 'ltsProminence'; ...
                    'ltsWidth'; 'ltsOnsetTime'; 'ltsTimeJitter'; ...
                    'ltsProbability'; 'spikesPerLts'; ...
                    'spikeMaxAmp'; 'spikeMinAmp'; ...
                    'spikeFrequency'; 'spikeAdaptation'
                    'burstOnsetTime'; 'burstTimeJitter'; ...
                    'burstProbability'; 'spikesPerBurst'};

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

% Extract output folder base name
outFolderName = extract_fileparts(outFolder, 'dirbase');

% Construct path to simulated LTS info
simSwpInfoPath = fullfile(outFolder, [outFolderName, '_', simSwpInfoSuffix, '.csv']);
stats3dPath = fullfile(outFolder, [outFolderName, '_', stats3dSuffix, '.mat']);

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
    % Decide on candidate parameters files
    [~, paramPaths] = all_files('Directory', outFolder, 'Suffix', 'params');

    % Extract the cell names
    cellNames = m3ha_extract_cell_name(paramPaths);

    % Display message
    fprintf('All sweeps from the following cells will be simulated: \n');
    print_cellstr(cellNames, 'OmitBraces', true, 'Delimiter', '\n');

    % Run simulations for each parameter file
    cellfun(@(x) m3ha_neuron_run_and_analyze (x, 'DataMode', dataMode, ...
                    'BuildMode', buildMode, 'SimMode', simMode, ...
                    'AttemptNumber', attemptNumber, ...
                    'SaveSimOutFlag', true, 'SaveLtsInfoFlag', true), ...
            paramPaths);

    % Find all simulated LTS stats spreadsheets
    [~, simLtsParamPaths] = ...
        all_files('Directory', outFolder, 'Recursive', true, ...
                    'Suffix', simLtsParamsSuffix, 'Extension', 'csv');

    % Copy over the spreadsheets
    copy_into(simLtsParamPaths, outFolder);
end

%% Combine LTS & burst feature tables
if combineFeatureTablesFlag
    % Display message
    fprintf('Combining LTS & burst statistics ... \n');

    % Find all simulated LTS stats spreadsheets
    [~, simLtsParamPaths] = ...
        all_files('Directory', outFolder, 'Suffix', simLtsParamsSuffix, ...
                    'Extension', 'csv');

    % Combine the spreadsheets
    simSwpInfo = vertcat_spreadsheets(simLtsParamPaths);

    % Rename variables
    simSwpInfo = renamevars(simSwpInfo, 'fileBase', 'simFileBase');

    % Extract the simulation file bases
    simFileBase = simSwpInfo.simFileBase;

    % Extract the original file bases
    fileBase = extractBefore(simFileBase, '_sim');

    % Make the original sweep name the row name
    simSwpInfo.Properties.RowNames = fileBase;

    % Load original sweep info
    origSwpInfo = m3ha_load_sweep_info;

    % Extract the condition info
    condInfo = origSwpInfo(:, condVarStrs);

    % Join the condition info
    simSwpInfo = join(simSwpInfo, condInfo, 'Keys', 'Row');

    % Save the simulated sweep info table
    writetable(simSwpInfo, simSwpInfoPath);
end

%% Compute statistics
if computeStatsFlag
    % Load sweep info
    simSwpInfo = readtable(simSwpInfoPath, 'ReadRowNames', true);

    % Compute statistics for all features
    disp('Computing statistics for 3D bar plots ...');
    statsTable = m3ha_compute_statistics('SwpInfo', simSwpInfo, ...
                                            'PharmConditions', pCond3D, ...
                                            'GIncrConditions', gCond3D, ...
                                            'DataMode', dataMode);

    % Save stats table
    save(stats3dPath, 'statsTable', '-v7.3');
end

%% Plot statistics
if plotStatsFlag
    % Load stats table
    disp('Loading statistics for 3D bar plots ...');
    load(stats3dPath, 'statsTable');

    % Restrict to measures of interest
    statsTable = statsTable(measuresOfInterest, :);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
