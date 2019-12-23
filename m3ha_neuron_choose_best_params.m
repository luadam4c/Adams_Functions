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
%       varargin    - 'SimMode': simulation mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'passive' - simulate a current pulse response
%                       'active'  - simulate an IPSC response
%                   default == 'active'
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
%                   - Any other parameter-value pair for 
%                           m3ha_neuron_run_and_analyze()
%
% Requires:
%       cd/argfun.m
%       cd/combine_strings.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/create_subplots.m
%       cd/extract_fields.m
%       cd/extract_param_values.m
%       cd/extract_vars.m
%       cd/find_in_strings.m
%       cd/force_column_cell.m
%       cd/force_matrix.m
%       cd/isemptycell.m
%       cd/istext.m
%       cd/load_params.m
%       cd/match_format_vector_sets.m
%       cd/m3ha_extract_cell_name.m
%       cd/m3ha_extract_iteration_string.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/plot_bar.m
%       cd/plot_tuning_curve.m
%       cd/save_all_figtypes.m
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
% 

%% Hard-coded parameters
validBuildModes = {'active', 'passive'};
validSimModes = {'active', 'passive'};

%   Note: The following must be consistent with compute_single_neuron_errors.m
idxSweep = 1;
idxMatch = 2;
idxAmp = 3;
idxTime = 4;
idxSlope = 5;

%   Note: The following must be consistent with 
%           m3ha_neuron_create_initial_params.m
% Names for each parameter
neuronParamNames = { ...
    'diamSoma', 'LDend', 'diamDend', ...
    'cm', 'Ra', 'corrD', 'gpas', 'epas', ...
    'pcabarITSoma', 'pcabarITDend1', 'pcabarITDend2', ...
    'shiftmIT', 'shifthIT', 'slopemIT', 'slopehIT', ...
    'ghbarIhSoma', 'ghbarIhDend1', 'ghbarIhDend2', 'ehIh', 'shiftmIh', ...
    'gkbarIKirSoma', 'gkbarIKirDend1', 'gkbarIKirDend2', ...
    'gkbarIASoma', 'gkbarIADend1', 'gkbarIADend2', ...
    'gnabarINaPSoma', 'gnabarINaPDend1', 'gnabarINaPDend2', ...
    };
cmInit = 0.88;          % specific membrane capacitance [uF/cm^2]
RaInit = 173;           % axial resistivity [Ohm-cm]
corrDInit = 1;          % dendritic surface area correction factor
% Lower bounds for each parameter
neuronParamsLowerBound = [ ...
    8, 5, 3, ...
    cmInit, RaInit, corrDInit, 1.0e-8, -95, ...
    1.0e-9, 1.0e-9, 1.0e-9, ...
    -30, -30, 0.1, 0.1, ...
    1.0e-9, 1.0e-9, 1.0e-9, -32, -30, ...
    1.0e-9, 1.0e-9, 1.0e-9, ...
    1.0e-9, 1.0e-9, 1.0e-9, ...
    1.0e-9, 1.0e-9, 1.0e-9, ...
    ];
% Upper bounds for each parameter
neuronParamsUpperBound = [ ...
    250, 150, 30, ...
    cmInit, RaInit, corrDInit, 1.0, -45, ...
    1.0e-2, 1.0e-2, 1.0e-2, ...
    30, 30, 10, 10, ...
    1.0e-2, 1.0e-2, 1.0e-2, -24, 30, ...
    1.0e-2, 1.0e-2, 1.0e-2, ...
    1.0e-2, 1.0e-2, 1.0e-2, ...
    1.0e-2, 1.0e-2, 1.0e-2, ...
    ];
% Whether each parameter will be varied log-scaled
neuronParamsIsLog = logical([ ...
    0, 0, 0, ...
    1, 1, 0, 1, 0, ...
    1, 1, 1, ...
    0, 0, 1, 1, ...
    1, 1, 1, 0, 0, ...
    1, 1, 1, ...
    1, 1, 1, ...
    1, 1, 1, ...
    ]);
neuronParamsYLimits = [neuronParamsLowerBound; neuronParamsUpperBound];

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
paramsToPlot = paramsToSave(1:15);
paramLabelsToPlot = paramsToPlot;

errorYLimits = [0, Inf];
errorXTicks = 'auto';
errorFigNumber = 1105;

errorParamXTicks = 'auto';
errorParamFigNumber = 1107;

% TODO: Make optional argument
errorParamSheetSuffix = 'error_param_table';
sheetExtension = 'csv';
figTypes = {'png'};

%% Default values for optional arguments
buildModeDefault = 'active';    % insert active channels by default
simModeDefault = 'active';      % simulate active responses by default
plotErrorHistoryFlagDefault = false;
plotErrorComparisonFlagDefault = false;
plotParamHistoryFlagDefault = false;
outFolderDefault = pwd;         % use the present working directory for outputs
                                %   by default
prefixDefault = '';             % set later

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
    @(x) validateattributes(x, {'cell', 'string'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BuildMode', buildModeDefault, ...
    @(x) any(validatestring(x, validBuildModes)));
addParameter(iP, 'SimMode', simModeDefault, ...
    @(x) any(validatestring(x, validSimModes)));
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

% Read from the Input Parser
parse(iP, candParamsTablesOrFiles, varargin{:});
buildMode = validatestring(iP.Results.BuildMode, validBuildModes);
simMode = validatestring(iP.Results.SimMode, validSimModes);
plotErrorHistoryFlag = iP.Results.PlotErrorHistoryFlag;
plotErrorComparisonFlag = iP.Results.PlotErrorComparisonFlag;
plotParamHistoryFlag = iP.Results.PlotParamHistoryFlag;
outFolder = iP.Results.OutFolder;
prefix = iP.Results.Prefix;

% Keep unmatched arguments for the m3ha_neuron_run_and_analyze() function
otherArguments = iP.Unmatched;

%% Preparation
% Parse first argument
if istext(candParamsTablesOrFiles)
    candParamsFiles = candParamsTablesOrFiles;
    candParamsTables = {};
else
    candParamsTables = candParamsTablesOrFiles;
    candParamsFiles = {};
end

% Load parameters if necessary
if isempty(candParamsTables)
    candParamsTables = cellfun(@load_params, candParamsTablesOrFiles, ...
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
    iterStr = create_labels_from_numbers(1:nTables, 'Prefix', 'table');
    cellName = repmat({'some_cell'}, nTables, 1);
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

% Decide on prefix if not provided
if isempty(prefix)
    % Use the cell name
    prefix = uniqueCellName;

    % If no cell name, use the directory base
    if isempty(prefix)
        prefix = extract_fileparts(outFolder, 'dirbase');
    end
end

% Create candidate labels
candLabel = combine_strings('Substrings', {prefix, 'from', iterStr});

%% Do the job
% Display message
fprintf('Choosing best parameters for %s ... \n', uniqueCellName);

% Compute errors for all tables
errorStructs = cellfun(@(x, y) m3ha_neuron_run_and_analyze(x, ...
                            'SaveImportLogFlag', false, ...
                            'PlotIndividualFlag', true, ...
                            'BuildMode', buildMode, 'SimMode', simMode, ...
                            'OutFolder', outFolder, ...
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

% Find the index of the table with the least error
[totalErrorBest, iTableBest] = min(totalError);

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
% Create full path to error sheet file
sheetBase = [prefix, '_', errorParamSheetSuffix];
sheetPathBase = fullfile(outFolder, sheetBase);
sheetPath = [sheetPathBase, '.', sheetExtension];

% Save the errors-parameters table
writetable(errorParamTable, sheetPath);

%% Plot error history
if plotErrorHistoryFlag
    % Display message
    fprintf('Plotting error history for %s ... \n', uniqueCellName);

    % Create figure title and file name
    errorFigTitle = ['Error History for ', uniqueCellName];
    errorFigName = strcat(sheetPathBase, '_error_history');

    % Plot error history
    plot_history_table(errorParamTable, errorsToPlot,...
                         iterNumber, [], [], errorYLimits, errorXTicks, ...
                         [], errorLabelsToPlot, errorFigTitle, ...
                         errorFigNumber, errorFigName, figTypes);
end

%% Plot error comparison
if plotErrorComparisonFlag
    % Display message
    fprintf('Plotting error comparison for %s ... \n', uniqueCellName);

    % Create figure
    fig2 = create_subplots(1, 1, 'FigNumber', 1106, 'ClearFigure', true);

    % Force as a matrix with each column corresponding to a type of error
    errorWeights = transpose(force_matrix(errorWeights));

    % Compute components of total error
    %   Note: must match groupLabels
    componentErrors = [ltsMatchError, avgSwpError, avgLtsAmpError, ...
                        avgLtsDelayError, avgLtsSlopeError] .* ...
            errorWeights(:, [idxMatch, idxSweep, idxAmp, idxTime, idxSlope]);
    
    % Decide on group labels
    %   Note: must match componentErrors
    groupLabels = {'LTS Match Error', 'Sweep Error', 'LTS Amp Error', ...
                    'LTS Time Error', 'LTS Slope Error'};

    % Decide on tick labels
    pTickLabels = iterStr;

    % Plot components of total error stacked
    plot_bar(componentErrors, 'BarDirection', 'horizontal', ...
            'ReverseOrder', true, 'GroupStyle', 'stacked', ...
            'PLabel', 'suppress', 'ReadoutLabel', 'Error (dimensionless)', ...
            'PTickLabels', pTickLabels, 'ColumnLabels', groupLabels, ...
            'FigTitle', ['Error Comparison for ', uniqueCellName]);

    % Save figure
    save_all_figtypes(fig2, strcat(sheetPathBase, '_error_comparison'), figTypes);
end

%% Plot parameter history
if plotParamHistoryFlag
    % Display message
    fprintf('Plotting parameter history for %s ... \n', uniqueCellName);

    % Create figure title and file name
    errorParamFigTitle = ['Error & Parameter History for ', uniqueCellName];
    errorParamFigName = strcat(sheetPathBase, '_param_history');

    % Decide on the errors and parameters to plot
    errorParamToPlot = vertcat(errorsToPlot(2:6), paramsToPlot);
    errorParamLabels = vertcat(errorLabelsToPlot(2:6), paramLabelsToPlot);

    % Get the original index for each parameter to plot
    indParamsToPlot = cellfun(@(x) find_in_strings(x, neuronParamNames, ...
                            'SearchMode', 'exact', 'MaxNum', 1), paramsToPlot);

    % Extract the y limits for each parameter
    paramYLimits = force_column_cell(neuronParamsYLimits(:, indParamsToPlot));

    % Determine whether each parameter should be plotted in a log scale
    paramIsLog = neuronParamsIsLog(indParamsToPlot);

    % Construct all error and parameter y limits
    errorParamYLimits = [repmat({errorYLimits}, 5, 1); paramYLimits];

    % Plot error & parameter history
    plot_history_table(errorParamTable, errorParamToPlot,...
                         iterNumber, paramIsLog, [], ...
                         errorParamYLimits, errorParamXTicks, ...
                         [], errorParamLabels, errorParamFigTitle, ...
                         errorParamFigNumber, errorParamFigName, figTypes);
end

%% Output results
% Return the table with the least error
bestParamsTable = candParamsTables{iTableBest};
bestParamsLabel = candLabel{iTableBest};

% Display result
fprintf('%s has the least error: %g!\n', bestParamsLabel, totalErrorBest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = plot_history_table(historyTable, errorParamToPlot, ...
                        iterNumbers, paramIsLog, xLimits, yLimits, xTicks, ...
                        colorMap, yLabels, figTitle, ...
                        figNumber, figName, figTypes)
%% Plots the history of variables given a table where each row is an iteration
% TODO: Pull out as its own function
% TODO: Make everything except historyTable optional arguments

%% Preparation
% Count the number of iterations
nIters = height(historyTable);

% Create iteration numbers if not provided
if isempty(iterNumbers)
    iterNumbers = transpose(1:nIters);
end

% TODO FOR SHINSHIN: count_strings.m
% Force character arrays as cell arrays
if ischar(errorParamToPlot)
    errorParamToPlot = {errorParamToPlot};
end

% Count the number of variables
nVarsToPlot = numel(errorParamToPlot);

% Decide on the number of rows for subplots
nSubplotRows = ceil(sqrt(nVarsToPlot));

% Compute the number of columns
nSubplotColumns = ceil(nVarsToPlot/nSubplotRows);

% Decide on axis limits
if isempty(xLimits)
    xLimits = [0, nIters + 1];
end

% Decide on color map
if isempty(colorMap)
    colorMap = {@lines};
end

% Decide on tick locations
if ischar(xTicks) && strcmp(xTicks, 'auto')
    xTicks = iterNumbers;
end

% Decide on whether to plot on a log scale
if isempty(paramIsLog)
    paramIsLog = repmat({false}, nVarsToPlot, 1);
elseif isnumeric(paramIsLog)
    paramIsLog = num2cell(paramIsLog);
end

%% Do the job
% Extract variables from table
dataToPlot = extract_vars(historyTable, errorParamToPlot);

% Match the number of items with dataToPlot
[paramIsLog, xLimits, yLimits, xTicks, colorMap, yLabels] = ...
    argfun(@(x) match_format_vector_sets(x, dataToPlot), ...
            paramIsLog, xLimits, yLimits, xTicks, colorMap, yLabels);

% Create figure
[fig, ax] = create_subplots(nSubplotRows, nSubplotColumns, ...
                'FigNumber', figNumber, 'ClearFigure', true, ...
                'FigExpansion', [nSubplotColumns / 2, nSubplotRows / 3]);

% Plot each variable on a separate subplot
dots = cellfun(@(a, b, c, d, e, f, g, h) ...
                    update_subplot(a, iterNumbers, b, c, d, e, f, g, h), ...
                num2cell(ax), dataToPlot, paramIsLog, xLimits, yLimits, ...
                xTicks, colorMap, yLabels, 'UniformOutput', false);

% Create an overarching title
suptitle(figTitle);

% Save figure
save_all_figtypes(fig, figName, figTypes);

%% Outputs
handles.fig = fig;
handles.ax = ax;
handles.dots = dots;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dots = update_subplot(axHandle, iterNumber, vecToPlot, ...
                        paramIsLog, xLimits, yLimits, xTicks, colorMap, yLabel)

% Create x axis label
xLabel = 'Iteration Number';

% Put the current subplot in focus
subplot(axHandle);

% Plot each iteration as a different color
dots = plot_tuning_curve(transpose(iterNumber), transpose(vecToPlot), ...
        'ReadoutIsLog', paramIsLog, ...
        'PLimits', xLimits, 'ReadOutLimits', yLimits, ...
        'PTicks', xTicks, 'ColorMap', colorMap, ...
        'PLabel', xLabel, 'ReadoutLabel', yLabel, ...
        'FigTitle', 'suppress', 'LegendLocation', 'suppress');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Convert error struct array to a table
errorTable = struct2table(errorStructs, 'AsArray', true);
% Make iterStr row names
errorTable.Properties.RowNames = iterStr;
% Add variables in the beginning
errorTable = addvars(errorTable, candLabel, cellName, ...
                        iterStr, 'Before', 1);
% Create candidate labels
candLabel = strcat(cellName, '_from_', iterStr);

subplot(3, 2, 1);
plot_tuning_curve(transpose(iterNumber), transpose(totalError), ...
        'PLimits', xLimits, 'ReadOutLimits', yLimits, ...
        'PTicks', xTicks, 'ColorMap', colorMap, ...
        'PLabel', 'Iteration Number', 'LegendLocation', 'suppress', ...
        'ReadoutLabel', 'Total Error', 'FigTitle', 'suppress');
% Plot the total error
subplot(3, 2, 1);
plot_tuning_curve(transpose(iterNumber), transpose(totalError), ...
        'PLimits', xLimits, 'ReadOutLimits', yLimits, ...
        'PTicks', xTicks, 'ColorMap', colorMap, ...
        'PLabel', 'Iteration Number', 'LegendLocation', 'suppress', ...
        'ReadoutLabel', 'Total Error', 'FigTitle', 'suppress');
% Plot the average sweep error
subplot(3, 2, 2);
plot_tuning_curve(transpose(iterNumber), transpose(avgSwpError), ...
        'PLimits', xLimits, 'ReadOutLimits', yLimits, ...
        'PTicks', xTicks, 'ColorMap', colorMap, ...
        'PLabel', 'Iteration Number', 'LegendLocation', 'suppress', ...
        'ReadoutLabel', 'Sweep Error', 'FigTitle', 'suppress');
% Plot the average LTS error
subplot(3, 2, 3);
plot_tuning_curve(transpose(iterNumber), transpose(ltsMatchError), ...
        'PLimits', xLimits, 'ReadOutLimits', yLimits, ...
        'PTicks', xTicks, 'ColorMap', colorMap, ...
        'PLabel', 'Iteration Number', 'LegendLocation', 'suppress', ...
        'ReadoutLabel', 'Match Error', 'FigTitle', 'suppress');
% Plot the average LTS amp error
subplot(3, 2, 4);
plot_tuning_curve(transpose(iterNumber), transpose(avgLtsAmpError), ...
        'PLimits', xLimits, 'ReadOutLimits', yLimits, ...
        'PTicks', xTicks, 'ColorMap', colorMap, ...
        'PLabel', 'Iteration Number', 'LegendLocation', 'suppress', ...
        'ReadoutLabel', 'Amp Error', 'FigTitle', 'suppress');
% Plot the average LTS time error
subplot(3, 2, 5);
plot_tuning_curve(transpose(iterNumber), transpose(avgLtsDelayError), ...
        'PLimits', xLimits, 'ReadOutLimits', yLimits, ...
        'PTicks', xTicks, 'ColorMap', colorMap, ...
        'PLabel', 'Iteration Number', 'LegendLocation', 'suppress', ...
        'ReadoutLabel', 'Time Error', 'FigTitle', 'suppress');
% Plot the average LTS slope error
subplot(3, 2, 6);
plot_tuning_curve(transpose(iterNumber), transpose(avgLtsSlopeError), ...
        'PLimits', xLimits, 'ReadOutLimits', yLimits, ...
        'PTicks', xTicks, 'ColorMap', colorMap, ...
        'PLabel', 'Iteration Number', 'LegendLocation', 'suppress', ...
        'ReadoutLabel', 'Slope Error', 'FigTitle', 'suppress');


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
