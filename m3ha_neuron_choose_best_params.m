function [bestParamsTable, bestParamsLabel, errorTable] = ...
                m3ha_neuron_choose_best_params (candParamsTablesOrFiles, varargin)
%% Chooses among candidates the NEURON parameters that fits a cell's data the best
% Usage: [bestParamsTable, bestParamsLabel] = ...
%               m3ha_neuron_choose_best_params (candParamsTablesOrFiles, varargin)
% Explanation:
%       Computes errors for more than one candidate sets of NEURON parameters
%            and choose the one with the least total error as the best 
%
% Example(s):
%       TODO
%
% Outputs:
%       bestParamsTable - the NEURON table for best parameters
%                       specified as a table
%       bestParamsLabel - file name or table name for the best parameters
%                       specified as a character vector
%
% Arguments:
%       candParamsTablesOrFiles  - candidate sets of NEURON parameter
%                                   tables or spreadsheet file names
%                   must be a cell array or string array
%       varargin    - 'UseHH': whether to use HH channels
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'SimMode': simulation mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'passive' - simulate a current pulse response
%                       'active'  - simulate an IPSC response
%                   default == 'active'
%                   - 'BuildMode': TC neuron build mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'passive' - insert leak channels only
%                       'active'  - insert both passive and active channels
%                       or a cell array of them TODO
%                   default == 'active'
%                   - 'PlotIndividualFlag': whether to plot fits
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotErrorHistoryFlag': whether to plot error history
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotErrorComparisonFlag': whether to plot 
%                                                   error comparison
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotParamHistoryFlag': whether to plot parameter history
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OutFolder': the directory where outputs will be placed
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'Prefix': prefix to prepend to file names
%                   must be a character array
%                   default == extract_common_prefix(fileBase)
%                   - 'FileNames': data file names to import
%                   must be a character vector, a string vector 
%                       or a cell array of character vectors
%                   default == none provided
%                   - 'IpscrWindow': IPSC response window in ms
%                   must be a numeric vector with 2 elements
%                   default == [0, 8000] + timeToStabilize
%                   - Any other parameter-value pair for 
%                           m3ha_neuron_run_and_analyze()
%
% Requires:
%       cd/argfun.m
%       cd/check_dir.m
%       cd/combine_strings.m
%       cd/create_error_for_nargin.m
%       cd/create_indices.m
%       cd/create_labels_from_numbers.m
%       cd/create_subplots.m
%       cd/extract_fields.m
%       cd/extract_fileparts.m
%       cd/extract_param_values.m
%       cd/force_matrix.m
%       cd/isemptycell.m
%       cd/istext.m
%       cd/read_params.m
%       cd/m3ha_extract_cell_name.m
%       cd/m3ha_extract_iteration_string.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/plot_bar.m
%       cd/plot_table_parallel.m
%       cd/save_all_figtypes.m
%       cd/save_params.m
%       cd/set_fields_zero.m
%
% Used by:
%       cd/m3ha_rank_neurons.m
%       /media/adamX/m3ha/optimizer4compgabab/singleneuronfitting63.m

% File History:
% 2019-11-23 Created by Adam Lu
% 2019-11-28 Now saves error table and plots individual plots for each set
%               of parameters
% 2019-12-19 Added 'PlotErrorHistoryFlag' as an optional argument
% 2019-12-27 Updated bounds
% 2019-04-28 Changed timeToStabilize from 2000 to 3000

%% Hard-coded parameters
validBuildModes = {'active', 'passive'};
validSimModes = {'active', 'passive'};

%   Note: The following must be consistent with compute_single_neuron_errors.m
idxSweep = 1;
idxMatch = 2;
idxAmp = 3;
idxTime = 4;
idxSlope = 5;

% The following must be consistent with both dclampDataExtractor.m & ...
%   singleneuron4compgabab.hoc
ipscrWinOrig = [0, 8000];       % IPSC response window (ms), original

% The following must be consistent with singleneuron4compgabab.hoc
timeToStabilize = 3000;         % padded time (ms) to make sure initial value 
                                %   of simulations are stabilized

% Spreadsheet settings
paramsToSave = { ...
    'diamSoma'; 'LDend'; 'diamDend'; 'gpas'; 'epas'; ...
    'pcabarITSoma'; 'pcabarITDend1'; 'pcabarITDend2'; ...
    'ghbarIhSoma'; 'ghbarIhDend1'; 'ghbarIhDend2'; ...
    'gkbarIKirSoma'; 'gkbarIKirDend1'; 'gkbarIKirDend2'; ...
    'gkbarIASoma'; 'gkbarIADend1'; 'gkbarIADend2'; ...
    'gnabarINaPSoma'; 'gnabarINaPDend1'; 'gnabarINaPDend2'; ...
    'cm'; 'Ra'; 'corrD'; ...
    'shiftmIT'; 'shifthIT'; 'slopemIT'; 'slopehIT'; ...
    'ehIh'; 'shiftmIh'; ...
    };

% Plot settings
errorsToPlot = {'totalError'; 'avgSwpError'; 'ltsMatchError'; ...
                'avgLtsAmpError'; 'avgLtsDelayError'; 'avgLtsSlopeError'};
errorLabelsToPlot = {'Total Error'; 'Sweep Error'; 'Match Error'; ...
            'Amp Error'; 'Time Error'; 'Slope Error'};
errorYLimits = [0, Inf];
errorXTicks = [];
errorFigNumber = 1105;

errorParamXTicks = [];
errorParamFigNumber = 1107;

% TODO: Make optional argument
bestParamsDirPrefix = 'bestparams';
errorParamSheetSuffix = 'error_param_table';
sheetExtension = 'csv';
figTypes = {'png'};

%% Default values for optional arguments
useHHDefault = false;           % don't use HH channels by default
buildModeDefault = 'active';    % insert active channels by default
simModeDefault = 'active';      % simulate active responses by default
plotIndividualFlagDefault = true;
plotErrorHistoryFlagDefault = false;
plotErrorComparisonFlagDefault = false;
plotParamHistoryFlagDefault = false;
outFolderDefault = pwd;         % use the present working directory for outputs
                                %   by default
prefixDefault = '';             % set later
fileNamesDefault = {};          % none provided by default
ipscrWindowDefault = ipscrWinOrig + timeToStabilize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'candParamsTablesOrFiles', ...
    @(x) validateattributes(x, {'cell', 'string', 'char'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'UseHH', useHHDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'BuildMode', buildModeDefault, ...
    @(x) any(validatestring(x, validBuildModes)));
addParameter(iP, 'SimMode', simModeDefault, ...
    @(x) any(validatestring(x, validSimModes)));
addParameter(iP, 'PlotIndividualFlag', plotIndividualFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotErrorHistoryFlag', plotErrorHistoryFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotErrorComparisonFlag', plotErrorComparisonFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotParamHistoryFlag', plotParamHistoryFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileNames', fileNamesDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['FileNames must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'IpscrWindow', ipscrWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));

% Read from the Input Parser
parse(iP, candParamsTablesOrFiles, varargin{:});
useHH = iP.Results.UseHH;
buildMode = validatestring(iP.Results.BuildMode, validBuildModes);
simMode = validatestring(iP.Results.SimMode, validSimModes);
plotIndividualFlag = iP.Results.PlotIndividualFlag;
plotErrorHistoryFlag = iP.Results.PlotErrorHistoryFlag;
plotErrorComparisonFlag = iP.Results.PlotErrorComparisonFlag;
plotParamHistoryFlag = iP.Results.PlotParamHistoryFlag;
outFolder = iP.Results.OutFolder;
prefix = iP.Results.Prefix;
fileNames = iP.Results.FileNames;
ipscrWindow = iP.Results.IpscrWindow;

% Keep unmatched arguments for the m3ha_neuron_run_and_analyze() function
otherArguments = iP.Unmatched;

%% Preparation
% Parse first argument
if istext(candParamsTablesOrFiles)
    candParamsFiles = candParamsTablesOrFiles;
    if ischar(candParamsFiles)
        candParamsFiles = {candParamsFiles};
    end
    candParamsTables = {};
else
    candParamsTables = candParamsTablesOrFiles;
    candParamsFiles = {};
end

% Load parameters if necessary
if isempty(candParamsTables)
    candParamsTables = cellfun(@read_params, candParamsFiles, ...
                                'UniformOutput', false);
end

% Count the number of tables
nTables = numel(candParamsTables);

% If there are no tables, return
if nTables == 0
    bestParamsTable = table.empty;
    bestParamsLabel = '';
    errorTable = table.empty;
    return
end

% Decide on iteration strings and cell names
%   Note: Make sure iterStr and cellName have nTables rows
if isempty(candParamsFiles)
    if ~isempty(fileNames)
        cellNameTemp = m3ha_extract_cell_name(fileNames{1});
    else
        cellNameTemp = 'some_cell';
    end
    cellName = repmat({cellNameTemp}, nTables, 1);

    iterStr = create_labels_from_numbers(1:nTables, 'Prefix', 'table');
else
    % Extract the chosen iteration string
    iterStr = m3ha_extract_iteration_string(candParamsFiles);
    if any(isemptycell(iterStr))
        error('Iteration string not found!');
    end

    % Extract the cell names
    cellName = m3ha_extract_cell_name(candParamsFiles);
end

% Get unique cell names
uniqueCellNames = unique(cellName);

% Check if all cell names are the same
if numel(uniqueCellNames) > 2
    error('Candidate parameters must all come from the same cell!');
else
    uniqueCellName = uniqueCellNames{1};
end

% Turn off all flags for stats and plots except plotIndividualFlag
% TODO: CHange this to use input parser for case-insensitivity
otherArguments = ...
    set_fields_zero(otherArguments, ...
        'saveLtsInfoFlag', 'saveLtsStatsFlag', ...
        'saveSimCmdsFlag', 'saveStdOutFlag', 'saveSimOutFlag', ...
        'plotConductanceFlag', 'plotCurrentFlag', ...
        'plotResidualsFlag', 'plotOverlappedFlag', ...
        'plotIpeakFlag', 'plotLtsFlag', 'plotStatisticsFlag', ...
        'plotSwpWeightsFlag');
otherArguments = ...
    set_fields_zero(otherArguments, ...
        'SaveLtsInfoFlag', 'SaveLtsStatsFlag', ...
        'SaveSimCmdsFlag', 'SaveStdOutFlag', 'SaveSimOutFlag', ...
        'PlotConductanceFlag', 'PlotCurrentFlag', ...
        'PlotResidualsFlag', 'PlotOverlappedFlag', ...
        'PlotIpeakFlag', 'PlotLtsFlag', 'PlotStatisticsFlag', ...
        'PlotSwpWeightsFlag');

% Extract name of output folder
outFolderName = extract_fileparts(outFolder, 'dirbase');

% Decide on prefix if not provided
if isempty(prefix)
    % Use the cell name
    prefix = uniqueCellName;

    % If no cell name, use the directory base
    if isempty(prefix)
        prefix = outFolderName;
    end
end

% Create full path to error sheet file
outPathPrefix = fullfile(outFolder, prefix);
sheetPathBase = [outPathPrefix, '_', errorParamSheetSuffix];
sheetPath = [sheetPathBase, '.', sheetExtension];

% Create best parameters directory name
bestParamsDirName = [bestParamsDirPrefix, '_', outFolderName];
bestParamsDir = fullfile(outFolder, bestParamsDirName);

% Check directory
check_dir(bestParamsDir);

%% Generate an error-params spreadsheet
if ~isfile(sheetPath)
    %% Run through all parameters
    % Display message
    fprintf('Choosing best parameters for %s ... \n', uniqueCellName);

    % Create candidate labels
    candLabel = combine_strings('Substrings', {prefix, 'from', iterStr});

    % Import data
    if ~isempty(fileNames)
        [realDataIpscr, sweepInfoIpscr] = ...
            m3ha_import_raw_traces(fileNames, 'ImportMode', 'active', ...
                        'Verbose', false, 'OutFolder', outFolder, ...
                        'ResponseWindow', ipscrWindow - timeToStabilize);

        nSweepsIpscr = height(sweepInfoIpscr);
        sweepInfoIpscr{:, 'holdCurrentNoise'} = zeros(nSweepsIpscr, 1);

        sweepInfoIpscr.Properties.VariableNames = ...
            strcat(sweepInfoIpscr.Properties.VariableNames, 'Ipscr');

        sweepInfoIpscrStruct = table2struct(sweepInfoIpscr, 'ToScalar', true);

        otherArguments = merge_structs(otherArguments, sweepInfoIpscrStruct);
    else
        realDataIpscr = [];
    end

    % Compute errors for all tables
    errorStructs = cellfun(@(x, y) m3ha_neuron_run_and_analyze(x, ...
                                'RealDataIpscr', realDataIpscr, ...
                                'FileNames', fileNames, ...
                                'IpscrWindow', ipscrWindow, ...
                                'SaveImportLogFlag', false, ...
                                'PlotIndividualFlag', plotIndividualFlag, ...
                                'BuildMode', buildMode, 'SimMode', simMode, ...
                                'UseHH', useHH, 'OutFolder', outFolder, ...
                                'Prefix', y, otherArguments), ...
                                candParamsTables, candLabel);

    % Extract scalar fields of interest
    %   Note: must be consistent with compute_single_neuron_errors.m
    [totalError, lts2SweepErrorRatio, match2FeatureErrorRatio, ...
            avgSwpError, avgLtsError, ltsMatchError, ...
            missedLtsError, falseLtsError, ...
            avgLtsAmpError, avgLtsDelayError, avgLtsSlopeError] = ...
        argfun(@(x) extract_fields(errorStructs, x, 'UniformOutput', true), ...
                'totalError', 'lts2SweepErrorRatio', 'match2FeatureErrorRatio', ...
                'avgSwpError', 'avgLtsError', 'ltsMatchError', ...
                'missedLtsError', 'falseLtsError', ...
                'avgLtsAmpError', 'avgLtsDelayError', 'avgLtsSlopeError');

    % Extract vector fields of interest
    %   Note: must be consistent with compute_single_neuron_errors.m
    [errorWeights, ltsFeatureWeights, swpErrors, ...
            ltsAmpErrors, ltsDelayErrors, ltsSlopeErrors] = ...
        argfun(@(x) extract_fields(errorStructs, x, 'UniformOutput', false), ...
                    'errorWeights', 'ltsFeatureWeights', 'swpErrors', ...
                    'ltsAmpErrors', 'ltsDelayErrors', 'ltsSlopeErrors');

    % Create iteration numbers
    iterNumber = transpose(1:nTables);

    % Add variables in the beginning
    errorTable = table(iterNumber, candLabel, cellName, iterStr, ...
                    totalError, avgSwpError, avgLtsError, ltsMatchError, ...
                    avgLtsAmpError, avgLtsDelayError, avgLtsSlopeError, ...
                    lts2SweepErrorRatio, missedLtsError, falseLtsError, ...
                    match2FeatureErrorRatio, ltsFeatureWeights, errorWeights, ...
                    swpErrors, ltsAmpErrors, ltsDelayErrors, ltsSlopeErrors, ...
                    'RowNames', iterStr);

    %% Join with parameters
    % Extract parameters of interest as a structure array 
    %   (each parameter is a field)
    paramValueStructs = extract_param_values(candParamsTables, ...
                                            'RowsToExtract', paramsToSave);

    % Convert to a table (each parameter is a variable)
    candParamValueTable = struct2table(paramValueStructs);

    % Add the iteration number column in the beginning
    candParamValueTable = addvars(candParamValueTable, iterNumber, 'Before', 1);

    % Join with the errorTable
    errorParamTable = join(errorTable, candParamValueTable);

    %% Save results
    % Save the errors-parameters table
    writetable(errorParamTable, sheetPath);
end

%% Plot error history
% TODO: FOR DEBUG but transfer to m3ha_plot_error_history.m
errorParamTable = readtable(sheetPath);
candLabel = errorParamTable.candLabel;
iterStr = errorParamTable.iterStr;
iterNumber = errorParamTable.iterNumber;

%% Make m3ha_plot_error_history.m
%   m3ha_plot_error_history(errorParamPath)
if plotErrorHistoryFlag
    % Display message
    fprintf('Plotting error history for %s ... \n', uniqueCellName);

    % Create figure title and file name
    errorFigTitle = ['Error History for ', uniqueCellName];
    errorFigName = strcat(outPathPrefix, '_error_history');

    % Plot error history
    plot_table_parallel(errorParamTable, 'VarsToPlot', errorsToPlot,...
                 'RowValues', iterNumber, 'RowLabel', 'Iteration Number', ...
                 'ReadoutLimits', errorYLimits, 'RowTicklocs', errorXTicks, ...
                 'ReadoutLabel', errorLabelsToPlot, 'FigTitle', errorFigTitle, ...
                 'FigNumber', errorFigNumber, 'FigName', errorFigName, ...
                 'FigTypes', figTypes);
end

%% Plot error comparison
if plotErrorComparisonFlag
    % Display message
    fprintf('Plotting error comparison for %s ... \n', uniqueCellName);

    % Create figure
    fig2 = create_subplots(1, 1, 'FigNumber', 1106, 'ClearFigure', true);

    % Extract component errors
    [componentErrors, groupLabels] = ...
        m3ha_extract_component_errors(errorParamTable);

    % Decide on ticks and tick labels
    pTicks = 1:numel(iterStr);
    pTickLabels = iterStr;

    % Plot components of total error stacked
    plot_bar(componentErrors, 'BarDirection', 'horizontal', ...
            'ReverseOrder', true, 'GroupStyle', 'stacked', ...
            'PLabel', 'suppress', 'ReadoutLabel', 'Error (dimensionless)', ...
            'PTicks', pTicks, 'PTickLabels', pTickLabels, ...
            'ColumnLabels', groupLabels, ...
            'FigTitle', ['Error Comparison for ', uniqueCellName]);

    % Save figure
    save_all_figtypes(fig2, strcat(outPathPrefix, '_error_comparison'), figTypes);
end

%% Plot parameter history
if plotParamHistoryFlag
    % Display message
    fprintf('Plotting parameter history for %s ... \n', uniqueCellName);

    % Create figure title and file name
    errorParamFigTitle = ['Error & Parameter History for ', uniqueCellName];
    errorParamFigName = strcat(outPathPrefix, '_param_history');

    % Decide on the errors and parameters to plot
    [errorParamToPlot, errorParamLabels, ...
            errorParamYLimits, errorParamIsLog] = m3ha_decide_on_plot_vars;

    % Decide on x tick locations
    if isempty(errorParamXTicks)
        % Count the number of iterations
        nIters = numel(iterNumber);

        % Set indices
        errorParamXTicks = create_indices([1, nIters], 'MaxNum', 18);
    end

    % Plot error & parameter history
    plot_table_parallel(errorParamTable, 'VarsToPlot', errorParamToPlot,...
             'XValues', iterNumber, 'VarIsLog', errorParamIsLog, ...
             'YLimits', errorParamYLimits, 'XTicks', errorParamXTicks, ...
             'XLabel', 'Iteration Number', 'YLabel', errorParamLabels, ...
             'FigTitle', errorParamFigTitle, 'FigTypes', figTypes, ...
             'FigNumber', errorParamFigNumber, 'FigName', errorParamFigName);
end

%% Choose best parameters
% Extract total errors
totalError = errorParamTable.totalError;

% Find the index of the table with the least error
[totalErrorBest, iTableBest] = min(totalError);

% Return the table with the least error
bestParamsTable = candParamsTables{iTableBest};
bestParamsLabel = candLabel{iTableBest};

% Save or copy the table with least error to the bestParamsDir
if isempty(candParamsFiles)
    % Create a new path to params table file
    bestParamsPath = fullfile(bestParamsDir, ...
                                [bestParamsDirName, '_', uniqueCellName]);

    % Save the file
    save_params(bestParamsTable, 'FileName', bestParamsPath);
else
    % Copy the file
    bestParamsFile = candParamsFiles{iTableBest};
    copy_into(bestParamsFile, bestParamsDir);
end

% Display result
fprintf('%s has the least error: %g!\n', bestParamsLabel, totalErrorBest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
