function handles = plot_table_parallel (myTable, varargin)
%% Plots the variables in a table in separate subplots
% Usage: handles = plot_table_parallel (myTable, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load_examples
%       plot_table_parallel(myTableNumeric)
%
% Outputs:
%       handles     - a structure with fields:
%                       fig     - figure handle
%                       ax      - array of subplot axes handles
%                       supAx   - handle for the super axes for the title
%                       hTuning - array of handles structures from 
%                                   plot_tuning_curve()
%                   specified as a scalar structure
%
% Arguments:
%       myTable    - a table of variables where each row is an iteration
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
%                   - 'SubplotDimensions': subplot dimensions
%                   must be empty or a numeric vector
%                   default == []
%                   - 'RowValues': x axis values corresponding to 
%                               each row of the table
%                   must be empty or a numeric vector
%                   default == rowsToPlot
%                   - 'VarIsLog': whether variable values are to be plotted 
%                               log-scaled
%                   must be a cell array or a numeric array of 
%                       logical/numeric binaries
%                   default == all false
%                   - 'RowLimits': x axis limits
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == set in plot_tuning_curve.m
%                   - 'ReadoutLimits': y axis limits
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                       or a cell array of them
%                   default == set in plot_tuning_curve.m
%                   - 'RowTickLocs': x tick values
%                   must be a numeric vector
%                   default == all x values
%                   - 'RowTickLabels': x tick labels
%                   must be a cell array of character vectors/strings
%                   default == row names or times if provided
%                   - 'RowLabel': label for the x axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                   default == 'Row Number'
%                   - 'ReadoutLabel': label(s) for the y axis, 
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
%                   - 'AxTitles': title for each subplot
%                   must be a string scalar or a character vector
%                   default == none
%                   - 'FigNumber': figure number for creating figure
%                   must be a positive integer scalar
%                   default == []
%                   - 'ClearFigure': whether to clear figure
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true if 'FigNumber' provided 
%                               but false otherwise
%                   - 'AlwaysNew': whether to always create a new figure even if
%                                   figNumber is not passed in
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
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
%       cd/resize_subplots_for_labels.m
%       cd/save_all_figtypes.m
%
% Used by:
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_rank_neurons.m

% File History:
% 2019-12-29 Moved from m3ha_neuron_choose_best_params.m
% 2019-12-30 Changed the default x label to 'row number'
% 2020-02-06 Added 'AxTitles' as an optional argument
% 2025-10-08 Added 'ClearFigure' and 'AlwaysNew' as optional arguments
% 2025-10-09 rowLimits default is now set in plot_tuning_curve.m
% TODO: Merge with plot_table.m
% TODO: 

%% Hard-coded parameters
defaultRowLabel = 'Row Number';

%% Default values for optional arguments
varsToPlotDefault = 'all';      % plot all variables by default
rowsToPlotDefault = 'all';      % plot all rows by default
subplotDimensionsDefault = [];  % set later
rowValuesDefault = [];          % set later
varIsLogDefault = [];           % set later
rowLimitsDefault = [];            % set later
readoutLimitsDefault = [];            % set later
rowTickLocsDefault = [];        % set later
rowTickLabelsDefault = {};      % set later
rowLabelDefault = '';           % set later
readoutLabelDefault = {};       % set later
colorMapDefault = [];           % set later
legendLocationDefault = 'suppress';
figTitleDefault = '';
axTitlesDefault = {''};
figNumberDefault = [];
clearFigureDefault = [];        % set later
alwaysNewDefault = false;       % don't always create new figure by default
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
addRequired(iP, 'myTable', ...
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
addParameter(iP, 'SubplotDimensions', subplotDimensionsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'RowValues', rowValuesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'VarIsLog', varIsLogDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric', 'cell'}, {'vector'}));
addParameter(iP, 'RowLimits', rowLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'ReadoutLimits', readoutLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2 || iscell(x));
addParameter(iP, 'RowTickLocs', rowTickLocsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'RowTickLabels', rowTickLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'RowLabel', rowLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ReadoutLabel', readoutLabelDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ColorMap', colorMapDefault);
addParameter(iP, 'LegendLocation', legendLocationDefault, ...
    @(x) all(islegendlocation(x, 'ValidateMode', true)));
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'AxTitles', axTitlesDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'FigNumber', figNumberDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'FigNumber must be a empty or a positive integer scalar!'));
addParameter(iP, 'ClearFigure', clearFigureDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'AlwaysNew', alwaysNewDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, myTable, varargin{:});
varsToPlot = iP.Results.VarsToPlot;
rowsToPlot = iP.Results.RowsToPlot;
subplotDimensions = iP.Results.SubplotDimensions;
rowValues = iP.Results.RowValues;
varIsLog = iP.Results.VarIsLog;
rowLimits = iP.Results.RowLimits;
readoutLimits = iP.Results.ReadoutLimits;
rowTickLocs = iP.Results.RowTickLocs;
rowTickLabels = iP.Results.RowTickLabels;
colorMap = iP.Results.ColorMap;
rowLabel = iP.Results.RowLabel;
readoutLabel = iP.Results.ReadoutLabel;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);
figTitle = iP.Results.FigTitle;
axTitles = iP.Results.AxTitles;
figNumber = iP.Results.FigNumber;
clearFigure = iP.Results.ClearFigure;
alwaysNew = iP.Results.AlwaysNew;
figName = iP.Results.FigName;
figTypes = iP.Results.FigTypes;

% Keep unmatched arguments for the plot_tuning_curve() function
otherArguments = iP.Unmatched;

%% Preparation
% Count the number of rows
nRowsOrig = height(myTable);

% Get all row names
allRowNames = myTable.Properties.RowNames;

% Restrict to rows to plot if requested
if isempty(rowsToPlot) || ischar(rowsToPlot) && strcmp(rowsToPlot, 'all')
    rowsToPlot = transpose(1:nRowsOrig);
else
    % Restrict to those rows
    myTable = myTable(rowsToPlot, :);

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
if isempty(rowValues)
    rowValues = transpose(1:nRows);
end

% Decide on the variables to plot
if ischar(varsToPlot) && strcmp(varsToPlot, 'all')
    varsToPlot = force_column_vector(myTable.Properties.VariableNames);
end

% Count the number of variables
nVarsToPlot = count_strings(varsToPlot);

if isempty(subplotDimensions)
    % Decide on the number of rows for subplots
    nSubplotRows = ceil(sqrt(nVarsToPlot));

    % Compute the number of columns
    nSubplotColumns = ceil(nVarsToPlot/nSubplotRows);
else
    nSubplotRows = subplotDimensions(1);
    nSubplotColumns = subplotDimensions(2);
end

% Decide on the x-axis label
if isempty(rowLabel)
    rowLabel = defaultRowLabel;
end

% Decide on the y-axis labels
if isempty(readoutLabel)
    readoutLabel = varsToPlot;
end

% Decide on color map
if isempty(colorMap)
    colorMap = {@lines};
end

% Decide on tick locations
if isempty(rowTickLocs)
    rowTickLocs = force_column_vector(rowValues);
end

% Decide on tick labels
if isempty(rowTickLabels)
    if ~isempty(allRowNames)
        rowLabels = allRowNames(rowsToPlot);
    else
        rowLabels = create_labels_from_numbers(rowsToPlot);
    end
    rowTickLabels = match_positions(rowLabels, rowValues, rowTickLocs);
	rowTickLabels = {rowTickLabels};
elseif iscell(rowTickLabels) && ~iscell(rowTickLabels{1})
    % Ensure rowTickLabels is a cell array of cell arrays
    rowTickLabels = {rowTickLabels}; 
end

% Decide on whether to plot on a log scale
if isempty(varIsLog)
    varIsLog = repmat({false}, nVarsToPlot, 1);
elseif ~iscell(varIsLog)
    varIsLog = num2cell(varIsLog);
end

%% Do the job
% Extract variables from table
dataToPlot = extract_vars(myTable, varsToPlot);

% Match the number of items with dataToPlot
[varIsLog, rowLimits, readoutLimits, rowTickLocs, rowTickLabels, ...
        colorMap, readoutLabel, axTitles] = ...
    argfun(@(x) match_format_vector_sets(x, dataToPlot), ...
            varIsLog, rowLimits, readoutLimits, rowTickLocs, rowTickLabels, ...
            colorMap, readoutLabel, axTitles);

% Decide whether to clear figure
if isempty(clearFigure)
    if ~isempty(figName)
        clearFigure = true;
    else
        clearFigure = false;
    end
end

% Create subplots
[fig, ax] = create_subplots(nSubplotRows, nSubplotColumns, ...
                'FigNumber', figNumber, 'ClearFigure', clearFigure, ...
                'AlwaysNew', alwaysNew, ...
                'FigExpansion', [nSubplotColumns / 2, nSubplotRows / 3]);

% Only use as many subplots as needed
axToUse = ax(1:nVarsToPlot);

% Delete extra subplots
if numel(ax) > nVarsToPlot
    delete(ax(nVarsToPlot + 1:end));
end
            
% Plot each variable on a separate subplot
hTuning = cellfun(@(a, b, c, d, e, f, g, h, i, j) ...
                update_subplot(a, rowValues, b, c, d, e, f, g, ...
                                rowLabel, h, i, j, otherArguments), ...
                num2cell(axToUse), dataToPlot, varIsLog, rowLimits, readoutLimits, ...
                rowTickLocs, colorMap, readoutLabel, rowTickLabels, axTitles, ...
                'UniformOutput', true);

% Create an overarching title
if ~isempty(figTitle)
    [supAx, ax] = resize_subplots_for_labels('FigTitle', figTitle);
else
    supAx = gobjects;
end

% Generate a legend if requested
if ~strcmpi(legendLocation, 'suppress')
    % TODO: Create overarching legend
    % lgd = legend(hTuning, 'location', legendLocation);
    % set(lgd, 'AutoUpdate', 'off', 'Interpreter', 'none');
end

% Save figure
if ~isempty(figName)
    save_all_figtypes(fig, figName, figTypes);
end

%% Outputs
handles.fig = fig;
handles.ax = ax;
handles.supAx = supAx;
handles.hTuning = hTuning;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hTuning = update_subplot(axHandle, iterNumber, vecToPlot, ...
                                varIsLog, rowLimits, readoutLimits, rowTickLocs, ...
                                colorMap, rowLabel, readoutLabel, rowTickLabels, ...
                                axTitle, otherArguments)

% Place a title for the current subplot
if ~isempty(axTitle)
    title(axHandle, axTitle);
end

% Plot each iteration as a different color
hTuning = plot_tuning_curve(transpose(iterNumber), transpose(vecToPlot), ...
                        'ReadoutIsLog', varIsLog, ...
                        'PLimits', rowLimits, 'ReadoutLimits', readoutLimits, ...
                        'PTicks', rowTickLocs, 'PTickLabels', rowTickLabels, ...
                        'PLabel', rowLabel, 'ReadoutLabel', readoutLabel, ...
                        'ColorMap', colorMap, 'FigTitle', 'suppress', ...
                        'LegendLocation', 'suppress', ...
                        'AxesHandle', axHandle, otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Decide on axis limits
if isempty(rowLimits)
    rowLimits = [min(rowValues) - 1, max(rowValues) + 1];
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
