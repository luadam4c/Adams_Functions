function handles = plot_table_parallel (historyTable, varargin)
%% Plots the variables in a table in separate subplots using markers
% Usage: handles = plot_table_parallel (historyTable, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       handles     - a structure with fields:
%                       fig
%                       ax
%                       dots
%                   specified as a scalar structure
%
% Arguments:
%       historyTable    - a table of variables where each row is an iteration
%                       must be a 2D table
%       varargin    - 'VarsToPlot': variables to plot
%                   must be a numeric array,
%                       a string scalar or a character vector, 
%                       or a cell array of character vectors
%                   default == 'all' (no restrictions)
%                   - 'RowsToPlot': rows to plot
%                   must be a numeric array,
%                       a string scalar or a character vector, 
%                       or a cell array of character vectors
%                   default == 'all' (no restrictions)
%                   - 'XValues': x axis values corresponding to 
%                               each row of the table
%                   must be empty or a numeric vector
%                   default == rowsToPlot
%                   - 'VarIsLog': whether variable values are to be plotted 
%                               log-scaled
%                   must be a cell array or a numeric array of 
%                       logical/numeric binaries
%                   default == all false
%                   - 'XLimits': x axis limits
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == [min(xValues) - 1, max(xValues) + 1]
%                   - 'YLimits': y axis limits
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                       or a cell array of them
%                   default == set in plot_tuning_curve.m
%                   - 'XTicks': x tick values
%                   must be a numeric vector
%                   default == all x values
%                   - 'XLabel': label for the x axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                   default == 'Iteration Number'
%                   - 'YLabel': label(s) for the y axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == varsToPlot
%                   - 'ColorMap' - color map used
%                   must be a 2-D numeric array with 3 columns
%                   default == @lines
%                   - 'LegendLocation': location for legend
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'auto'      - use default
%                       'suppress'  - no legend
%                       anything else recognized by the legend() function
%                   default == 'suppress'
%                   - 'FigTitle': title for the figure
%                   must be a string scalar or a character vector
%                   default == none
%                   - 'FigNumber': figure number for creating figure
%                   must be a positive integer scalar
%                   default == []
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - Any other parameter-value pair for plot_tuning_curve()
%
% Requires:
%       cd/count_strings.m
%       cd/create_labels_from_numbers.m
%       cd/create_subplots.m
%       cd/create_error_for_nargin.m
%       cd/extract_vars.m
%       cd/find_first_match.m
%       cd/force_column_vector.m
%       cd/islegendlocation.m
%       cd/ispositiveintegervector.m
%       cd/match_format_vector_sets.m
%       cd/match_positions.m
%       cd/plot_tuning_curve.m
%       cd/save_all_figtypes.m
%
% Used by:
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_rank_neurons.m

% File History:
% 2019-12-29 Moved from m3ha_neuron_choose_best_params.m
% TODO: Rename as plot_comparison_table?
% TODO: Merge with plot_table.m

%% Hard-coded parameters
defaultXLabel = 'Iteration Number';

% TODO: Make optional argument
xTickLabels = {};

%% Default values for optional arguments
varsToPlotDefault = 'all';      % plot all variables by default
rowsToPlotDefault = 'all';      % plot all rows by default
xValuesDefault = [];            % set later
varIsLogDefault = [];           % set later
xLimitsDefault = [];            % set later
yLimitsDefault = [];            % set later
xTicksDefault = [];             % set later
xLabelDefault = '';             % set later
yLabelDefault = {};             % set later
colorMapDefault = [];           % set later
legendLocationDefault = 'suppress';
figTitleDefault = '';
figNumberDefault = [];
figNameDefault = '';
figTypesDefault = 'png';

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
addRequired(iP, 'historyTable', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'VarsToPlot', varsToPlotDefault, ...
    @(x) assert(ispositiveintegervector(x) || ischar(x) || ...
                    iscellstr(x) || isstring(x), ...
                ['VarsToPlot must be either a positive integer vector, ', ...
                    'a string array or a cell array of character arrays!']));
addParameter(iP, 'RowsToPlot', rowsToPlotDefault, ...
    @(x) assert(ispositiveintegervector(x) || ischar(x) || ...
                    iscellstr(x) || isstring(x), ...
                ['RowsToPlot must be either a positive integer vector, ', ...
                    'a string array or a cell array of character arrays!']));
addParameter(iP, 'XValues', xValuesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'VarIsLog', varIsLogDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric', 'cell'}, {'vector'}));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2 || iscell(x));
addParameter(iP, 'XTicks', xTicksDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'XLabel', xLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ColorMap', colorMapDefault);
addParameter(iP, 'LegendLocation', legendLocationDefault, ...
    @(x) all(islegendlocation(x, 'ValidateMode', true)));
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigNumber', figNumberDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'FigNumber must be a empty or a positive integer scalar!'));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, historyTable, varargin{:});
varsToPlot = iP.Results.VarsToPlot;
rowsToPlot = iP.Results.RowsToPlot;
xValues = iP.Results.XValues;
varIsLog = iP.Results.VarIsLog;
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;
xTicks = iP.Results.XTicks;
colorMap = iP.Results.ColorMap;
xLabel = iP.Results.XLabel;
yLabel = iP.Results.YLabel;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);
figTitle = iP.Results.FigTitle;
figNumber = iP.Results.FigNumber;
figName = iP.Results.FigName;
figTypes = iP.Results.FigTypes;

% Keep unmatched arguments for the plot_tuning_curve() function
otherArguments = iP.Unmatched;

%% Preparation
% Count the number of rows
nRowsOrig = height(historyTable);

% Get all row names
allRowNames = historyTable.Properties.RowNames;

% Restrict to rows to plot if requested
if ischar(rowsToPlot) && strcmp(rowsToPlot, 'all')
    rowsToPlot = transpose(1:nRowsOrig);
else
    % Restrict to those rows
    historyTable = historyTable(rowsToPlot, :);

    % Convert rowsToPlot to numeric values
    if ~isnumeric(rowsToPlot)
        if ~isempty(allRowNames)
            rowsToPlot = find_first_match(rowsToPlot, allRowNames, ...
                                'MatchMode', 'exact', 'IgnoreCase', false);
        else
            error('rowsToPlot can''t be text if row names are not present!');
        end
    end
end

% Count the new number of rows
nRows = numel(rowsToPlot);

% Create x values if not provided
if isempty(xValues)
    xValues = transpose(1:nRows);
end

% Decide on the variables to plot
if ischar(varsToPlot) && strcmp(varsToPlot, 'all')
    varsToPlot = historyTable.Properties.VariableNames;
end

% Count the number of variables
nVarsToPlot = count_strings(varsToPlot);

% Decide on the number of rows for subplots
nSubplotRows = ceil(sqrt(nVarsToPlot));

% Compute the number of columns
nSubplotColumns = ceil(nVarsToPlot/nSubplotRows);

% Decide on axis limits
if isempty(xLimits)
    xLimits = [min(xValues) - 1, max(xValues) + 1];
end

% Decide on the x-axis label
if isempty(xLabel)
    xLabel = defaultXLabel;
end

% Decide on the y-axis labels
if isempty(yLabel)
    yLabel = varsToPlot;
end

% Decide on color map
if isempty(colorMap)
    colorMap = {@lines};
end

% Decide on tick locations
if isempty(xTicks)
    xTicks = force_column_vector(xValues);
end

% Decide on tick labels
if isempty(xTickLabels)
    if ~isempty(allRowNames)
        rowLabels = allRowNames(rowsToPlot);
    else
        rowLabels = create_labels_from_numbers(rowsToPlot);
    end
    xTickLabels = {match_positions(rowLabels, xValues, xTicks)};
end

% Decide on whether to plot on a log scale
if isempty(varIsLog)
    varIsLog = repmat({false}, nVarsToPlot, 1);
elseif ~iscell(varIsLog)
    varIsLog = num2cell(varIsLog);
end

%% Do the job
% Extract variables from table
dataToPlot = extract_vars(historyTable, varsToPlot);

% Match the number of items with dataToPlot
[varIsLog, xLimits, yLimits, xTicks, xTickLabels, colorMap, yLabel] = ...
    argfun(@(x) match_format_vector_sets(x, dataToPlot), ...
            varIsLog, xLimits, yLimits, xTicks, xTickLabels, colorMap, yLabel);

% Decide whether to clear figure
if ~isempty(figName)
    clearFigure = true;
else
    clearFigure = false;
end

% Create subplots
[fig, ax] = create_subplots(nSubplotRows, nSubplotColumns, ...
                'FigNumber', figNumber, 'ClearFigure', clearFigure, ...
                'FigExpansion', [nSubplotColumns / 2, nSubplotRows / 3]);

% Plot each variable on a separate subplot
dots = cellfun(@(a, b, c, d, e, f, g, h, i) ...
                update_subplot(a, xValues, b, c, d, e, f, g, ...
                                xLabel, h, i, otherArguments), ...
                num2cell(ax), dataToPlot, varIsLog, xLimits, yLimits, ...
                xTicks, colorMap, yLabel, xTickLabels, 'UniformOutput', false);

% Create an overarching title
if ~isempty(figTitle)
    suptitle(figTitle);
end

% Generate a legend if requested
if ~strcmpi(legendLocation, 'suppress')
    % lgd = legend(dots, 'location', legendLocation);
    % set(lgd, 'AutoUpdate', 'off', 'Interpreter', 'none');
end

% Save figure
if ~isempty(figName)
    save_all_figtypes(fig, figName, figTypes);
end

%% Outputs
handles.fig = fig;
handles.ax = ax;
handles.dots = dots;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dots = update_subplot(axHandle, iterNumber, vecToPlot, ...
                                varIsLog, xLimits, yLimits, xTicks, ...
                                colorMap, xLabel, yLabel, xTickLabels, ...
                                otherArguments)

% Put the current subplot in focus
subplot(axHandle);

% Plot each iteration as a different color
dots = plot_tuning_curve(transpose(iterNumber), transpose(vecToPlot), ...
                        'ReadoutIsLog', varIsLog, ...
                        'PLimits', xLimits, 'ReadOutLimits', yLimits, ...
                        'PTicks', xTicks, 'PTickLabels', xTickLabels, ...
                        'PLabel', xLabel, 'ReadoutLabel', yLabel, ...
                        'ColorMap', colorMap, 'FigTitle', 'suppress', ...
                        'LegendLocation', 'suppress', otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
