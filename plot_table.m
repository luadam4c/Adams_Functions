function handles = plot_table (myTable, varargin)
%% Plots variables (columns) in a table
% Usage: handles = plot_table (myTable, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load_examples
%       plot_table(myTableNumeric)
%       plot_table(myTableNumeric, 'PlotMode', 'parallel')
%       plot_table(myTableNumeric, 'PlotMode', 'separate')
%
% Outputs:
%       handles     - structure with fields dependent on 'PlotMode':
%                   For 'overlapped' mode (from plot_tuning_curve.m):
%                       fig         - figure handle
%                       ax          - axes handle
%                       curves      - handles to the plotted tuning curves
%                       confInts    - (optional) handles to confidence intervals
%                       boundaries  - (optional) handles to boundary lines
%                       selected    - (optional) handles to selected value markers
%                       averages    - (optional) handles to phase average lines
%                       avgWindows  - (optional) handles to average window bars
%                   For 'parallel' mode (from plot_table_parallel.m):
%                       fig     - figure handle
%                       ax      - array of subplot axes handles
%                       supAx   - handle for the super axes for the title
%                       hTuning - array of handles structures from 
%                                   plot_tuning_curve()
%                   For 'separate' mode (from plot_struct.m):
%                       A structure array where each element contains the
%                       handles for a single plot, with fields as described
%                       for the 'overlapped' mode.
%                   specified as a structure or structure array
%
% Arguments:
%       myTable     - a table with variables to plot
%                   must be a table
%       varargin    - 'PlotType': type of plot
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'tuning'    - circles
%                       'bar'       - horizontal bars
%                   default == 'tuning'
%                   - 'PlotMode': plotting mode for multiple columns
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'overlapped'    - overlapped in a single plot
%                       'parallel'      - in parallel in subPlots
%                       'separate'      - in separate figures
%                   default == 'overlapped'
%                   - 'VarsToPlot': variable (column) names of the table
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == plot all variables
%                   - 'VarLabels': variable labels
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == use distinct parts of variable names
%                   - 'DistinctParts': whether to extract distinct parts
%                                       or variable labels
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'RowsToPlot': rows to plot
%                   must be a numeric array,
%                       a string scalar or a character vector, 
%                       or a cell array of character vectors
%                   default == 'all' (no restrictions)
%                   - 'RowValues': x axis values corresponding to 
%                               each row of the table
%                   must be empty or a numeric vector
%                   default == rowsToPlot
%                   - 'PhaseVariables': variable (column) names for phases
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == none
%                   - 'Delimiter': delimiter used for extracting distinct parts
%                   must be a string scalar or a character vector
%                   default == '_'
%                   - 'ReadoutLabel': label for the readout
%                   must be a string scalar or a character vector
%                   default == set by plot_tuning_curve.m
%                   - 'TableLabel': label for the table
%                   must be a string scalar or a character vector
%                   default == either a common prefix from variable names
%                               or the input table variable name
%                   - 'RowLabel': label for the rows
%                   must be a string scalar or a character vector
%                   default == none ('suppress')
%                   - 'RowTickLocs': x tick locations
%                   must be a numeric vector
%                   default == 1:numel(rowTickLabels)
%                   - 'RowTickLabels': x tick labels
%                   must be a cell array of character vectors/strings
%                   default == row names or times if provided
%                   - 'FigTitle': title for the figure
%                   must be a string scalar or a character vector
%                   default == none
%                   - 'ClearFigure': whether to clear figure
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true if 'FigNumber' provided 
%                               but false otherwise
%                   - 'ShowFigure': whether to show figure
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'AlwaysNew': whether to always create a new figure even if
%                                   figNumber is not passed in
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OutFolder': output folder if FigNames not set
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - Any other parameter-value pair for the plot_struct() 
%                       or the plot_tuning_curve() 
%                       or the plot_table_parallel() function
%
% Requires:
%       cd/char2rgb.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/extract_common_directory.m
%       cd/extract_common_prefix.m
%       cd/extract_fileparts.m
%       cd/plot_struct.m
%       cd/plot_traces.m
%       cd/plot_tuning_curve.m
%
% Used by:
%       cd/plot_measures.m
%       cd/parse_multiunit.m
%       cd/plot_repetitive_protocols.m

% File History:
% 2018-12-18 Moved from plot_repetitive_protocols.m
% 2018-12-18 Now uses iP.KeepUnmatched
% 2018-12-18 Now uses extract_common_directory.m
% 2018-12-18 Now uses row names without processing if not file names
% 2019-03-17 Deal with timetables differently for RowNames
% 2019-03-17 Implemented plotting together (use plot_tuning_curve directly)
% 2019-03-17 Added 'PlotSeparately' as an optional argument
% 2019-03-25 Added 'PhaseVariables' as an optional argument
% 2019-05-08 Added 'PlotType' as an optional argument
% 2019-08-07 Added 'RowTickLabels' as an optional argument
% 2019-08-07 Added 'RowTickLocs' as an optional argument
% 2019-12-30 Changed 'PlotSeparately' to 'PlotMode'
% 2025-10-07 Implemented 'parallel' plot mode by calling plot_table_parallel.m
% 2025-10-07 Now uses create_labels_from_numbers.m
% 2025-10-08 Added 'ClearFigure' and 'AlwaysNew' as optional arguments
% 2025-10-17 Added 'ShowFigure' as an optional argument with default true
% TODO: Transfer 'VarIsLog', 'SubplotDimensions', 'FigTypes' from plot_table_parallel.m to plot_table.m
%           and update m3ha_neuron_choose_best_params.m and m3ha_rank_neurons.m to use plot_table.m
% TODO: Return handles to plots
% TODO: Pass in figNames or figNumbers when plotting separately
% 

%% Hard-coded parameters
validPlotTypes = {'tuning', 'bar'};
validPlotModes = {'overlapped', 'parallel', 'separate'};

lineSpecOverlapped = '-';
lineWidthOverlapped = 1;
markerEdgeColorOverlapped = char2rgb('DarkOrchid');
markerFaceColorOverlapped = char2rgb('LightSkyBlue');

lineSpecParallel = '-';
lineWidthParallel = 1;
markerEdgeColorParallel = char2rgb('DarkOrchid');
markerFaceColorParallel = char2rgb('LightSkyBlue');

lineSpecSeparate = 'o';
lineWidthSeparate = 1;
markerEdgeColorSeparate = char2rgb('DarkOrchid');
markerFaceColorSeparate = char2rgb('LightSkyBlue');

%% Default values for optional arguments
plotTypeDefault = 'tuning';
plotModeDefault = 'overlapped'; % plot columns overlapped by default
lineSpecDefault = '';
lineWidthDefault = [];
markerEdgeColorDefault = [];
markerFaceColorDefault = [];
varsToPlotDefault = {};         % plot all variables by default
varLabelsDefault = {};          % set later
distinctPartsDefault = true;    % extract distinct parts of variable names
                                %   by default
rowsToPlotDefault = {};         % plot all rows by default
rowValuesDefault = [];          % set later
phaseVariablesDefault = {};     % no phases by default
delimiterDefault = '_';         % use '_' as delimiter by default
readoutLabelDefault = '';       % set later
tableLabelDefault = '';         % set later
rowLabelDefault = 'suppress';   % No x label by default
rowTickLocsDefault = [];
rowTickLabelsDefault = {};
figTitleDefault = '';
clearFigureDefault = [];        % set later
showFigureDefault = true;       % show figure by default
alwaysNewDefault = [];          % set later
outFolderDefault = pwd;
figNameDefault = '';

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
    @(x) validateattributes(x, {'table', 'timetable'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PlotType', plotTypeDefault, ...
    @(x) any(validatestring(x, validPlotTypes)));
addParameter(iP, 'PlotMode', plotModeDefault, ...
    @(x) any(validatestring(x, validPlotModes)));
addParameter(iP, 'LineSpec', lineSpecDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'LineWidth', lineWidthDefault);
addParameter(iP, 'MarkerEdgeColor', markerEdgeColorDefault);
addParameter(iP, 'MarkerFaceColor', markerFaceColorDefault);
addParameter(iP, 'VarsToPlot', varsToPlotDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['VarsToPlot must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'VarLabels', varLabelsDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['VarLabels must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'DistinctParts', distinctPartsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RowsToPlot', rowsToPlotDefault, ...
    @(x) assert(ispositiveintegervector(x) || ischar(x) || ...
                    iscellstr(x) || isstring(x), ...
                ['RowsToPlot must be either a positive integer vector, ', ...
                    'a string array or a cell array of character arrays!']));
addParameter(iP, 'RowValues', rowValuesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'PhaseVariables', phaseVariablesDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['VarsToPlot must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ReadoutLabel', readoutLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TableLabel', tableLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'RowLabel', rowLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'RowTickLocs', rowTickLocsDefault, ...
    @(x) isempty(x) || isnumericvector(x));
addParameter(iP, 'RowTickLabels', rowTickLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ClearFigure', clearFigureDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ShowFigure', showFigureDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'AlwaysNew', alwaysNewDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, myTable, varargin{:});
plotType = validatestring(iP.Results.PlotType, validPlotTypes);
plotMode = validatestring(iP.Results.PlotMode, validPlotModes);
lineSpec = iP.Results.LineSpec;
lineWidth = iP.Results.LineWidth;
markerEdgeColor = iP.Results.MarkerEdgeColor;
markerFaceColor = iP.Results.MarkerFaceColor;
varsToPlot = iP.Results.VarsToPlot;
varLabels = iP.Results.VarLabels;
distinctParts = iP.Results.DistinctParts;
rowsToPlot = iP.Results.RowsToPlot;
rowValues = iP.Results.RowValues;
phaseVariables = iP.Results.PhaseVariables;
delimiter = iP.Results.Delimiter;
readoutLabel = iP.Results.ReadoutLabel;
tableLabel = iP.Results.TableLabel;
rowLabel = iP.Results.RowLabel;
rowTickLocs = iP.Results.RowTickLocs;
rowTickLabels = iP.Results.RowTickLabels;
figTitle = iP.Results.FigTitle;
clearFigure = iP.Results.ClearFigure;
showFigure = iP.Results.ShowFigure;
alwaysNew = iP.Results.AlwaysNew;
outFolder = iP.Results.OutFolder;
figName = iP.Results.FigName;

% Keep unmatched arguments for the plotting function
otherArguments = iP.Unmatched;

%% Preparation
% Check if output directory exists
check_dir(outFolder);

% Restrict to variables to plot or extract the variable names
if ~isempty(varsToPlot)
    tableToPlot = myTable(:, varsToPlot);
else
    tableToPlot = myTable;
    varsToPlot = myTable.Properties.VariableNames;
end

% Count the number of rows
nRowsOrig = height(myTable);

% Get all row names
if isprop(myTable.Properties, 'RowNames') && ...
        iscell(myTable.Properties.RowNames) && ...
        ~isempty(myTable.Properties.RowNames)
    % Get the row names
    allRowNames = myTable.Properties.RowNames;
else
    % Else set as empty
    allRowNames = {};
end

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

% If provided, make sure there are an equal number of phase variables
%   and extract the phase vectors for each variable
if ~isempty(phaseVariables)
    % Force as a column cell array
    phaseVariables = force_column_cell(phaseVariables);

    % Match the number of phase variables to the number of variables to plot
    phaseVariables = match_row_count(phaseVariables, numel(varsToPlot));

    % Extract the phase vectors
    phaseVectors = cellfun(@(x) myTable{:, x}, phaseVariables, ...
                            'UniformOutput', false);
else
    phaseVectors = {};
end

% Extract distinct parts if requested
if isempty(varLabels)
    if distinctParts
        varLabels = extract_fileparts(varsToPlot, 'distinct', 'Delimiter', delimiter);
    else
        varLabels = varsToPlot;
    end
end

% Decide on table label
if isempty(tableLabel)
    % First try to extract a common prefix from the variables to plot
    tableLabel = extract_common_prefix(varsToPlot, 'Delimiter', delimiter);

    % If no such prefix exists, use the table variable name
    if isempty(tableLabel)
        tableLabel = inputname(1);
    end
end

% Decide on x values
if isempty(rowValues)
    if istimetable(tableToPlot)
        % Extract time
        rowValues = tableToPlot.Properties.RowTimes;
    else
        % Count rows
        nRows = height(tableToPlot);

        % Use row numbers
        rowValues = transpose(1:nRows);
    end
end

% Decide on rowTickLabels
if isempty(rowTickLabels)
    if ~isempty(allRowNames)
        % If all row names are file names, process them
        %   Otherwise, just use the row names as the x tick labels
        if all(isfile(allRowNames))
            % Extract the distinct file bases
            rowTickLabels = extract_fileparts(allRowNames, 'distinct');

            % Replace all instances of '_' with '\_'
            rowTickLabels = replace(rowTickLabels, '_', '\_');
        else
            % Just use the row names
            rowTickLabels = allRowNames;
        end
    elseif istimetable(myTable)
        rowTickLabels = {};
    else
        % Use default x tick labels
        if isempty(rowLabel)
            rowTickLabels = create_labels_from_numbers(rowsToPlot, 'Prefix', 'Row #');
        else
            rowTickLabels = create_labels_from_numbers(rowsToPlot);
        end
    end
end

% Decide on lineSpec
if isempty(lineSpec)
    switch plotMode
    case 'overlapped'
        lineSpec = lineSpecOverlapped;
    case 'parallel'
        lineSpec = lineSpecParallel;
    case 'separate'
        lineSpec = lineSpecSeparate;
    end
end

% Decide on lineWidth
if isempty(lineWidth)
    switch plotMode
    case 'overlapped'
        lineWidth = lineWidthOverlapped;
    case 'parallel'
        lineWidth = lineWidthParallel;
    case 'separate'
        lineWidth = lineWidthSeparate;
    end
end

% Decide on markerEdgeColor
if isempty(markerEdgeColor)
    switch plotMode
    case 'overlapped'
        markerEdgeColor = markerEdgeColorOverlapped;
    case 'parallel'
        markerEdgeColor = markerEdgeColorParallel;
    case 'separate'
        markerEdgeColor = markerEdgeColorSeparate;
    end
end

% Decide on markerFaceColor
if isempty(markerFaceColor)
    switch plotMode
    case 'overlapped'
        markerFaceColor = markerFaceColorOverlapped;
    case 'parallel'
        markerFaceColor = markerFaceColorParallel;
    case 'separate'
        markerFaceColor = markerFaceColorSeparate;
    end
end

% Decide on whether to create a new figure
if isempty(alwaysNew)
    switch plotMode
    case {'overlapped', 'parallel'}
        alwaysNew = false;
    case 'separate'
        alwaysNew = true;
    end    
end

% Decide on figure title
if isempty(figTitle)
    switch plotMode
    case {'overlapped', 'parallel'}
        if istimetable(tableToPlot)
            figTitle = replace([tableLabel, ' over time'], '_', ' ');
        else
            figTitle = '';
        end
    case 'separate'
        figTitle = varLabels;
    end
end

% Decide on figure name(s)
if isempty(figName)
    switch plotMode
    case {'overlapped', 'parallel'}
        figName = fullfile(outFolder, [tableLabel, '.png']);
    case 'separate'
        figName = fullfile(outFolder, strcat(varLabels, '.png'));
    end
end

%% Do the job
switch plotMode
case 'overlapped'
    % Convert to an array
    if istimetable(tableToPlot)
        % Extract variables
        myArray = tableToPlot.Variables;
    else
        % Use table2array
        myArray = table2array(tableToPlot);
    end

    % Decide on readout label
    if isempty(readoutLabel)
        readoutLabel = replace(tableLabel, '_', ' ');
    end

    % Decide on the figure handle
    set_figure_properties('ClearFigure', clearFigure, 'AlwaysNew', alwaysNew, ...
                            'ShowFigure', showFigure);

    % Plot tuning curve(s)
    switch plotType
    case 'tuning'
        handles = plot_tuning_curve(rowValues, myArray, ...
                        'PhaseVectors', phaseVectors, ...
                        'PTicks', rowTickLocs, 'PTickLabels', rowTickLabels, ...
                        'PLabel', rowLabel, ...
                        'ReadoutLabel', readoutLabel, ...
                        'ColumnLabels', varLabels, ...
                        'FigTitle', figTitle, ...
                        'FigName', figName, ...
                        'LineSpec', lineSpec, 'LineWidth', lineWidth, ...
                        'MarkerEdgeColor', markerEdgeColor, ...
                        'MarkerFaceColor', markerFaceColor, ...
                        otherArguments);
    case 'bar'
        % TODO
    otherwise
        error('plotType unrecognized!')
    end
case 'parallel'
    % Decide on readout label
    if isempty(readoutLabel)
        readoutLabel = varLabels;
    end

    % Plot variables in parallel subplots
    handles = plot_table_parallel(tableToPlot, ...
                    'RowValues', rowValues, ...
                    'RowTickLocs', rowTickLocs, ...
                    'RowTickLabels', rowTickLabels, ...
                    'RowLabel', rowLabel, ...
                    'ReadoutLabel', readoutLabel, ...
                    'LineSpec', lineSpec, 'LineWidth', lineWidth, ...
                    'MarkerEdgeColor', markerEdgeColor, ...
                    'MarkerFaceColor', markerFaceColor, ...
                    'FigTitle', figTitle, ...
                    'AxTitles', varLabels, ...
                    'FigName', figName, ...
                    'ClearFigure', clearFigure, 'AlwaysNew', alwaysNew, ...
                    'ShowFigure', showFigure, ...
                    otherArguments);
case 'separate'
    % Convert to a structure array
    myStruct = table2struct(tableToPlot);

    % Plot fields
    handles = plot_struct(myStruct, ...
                        'PValues', rowValues, ...
                        'PhaseVectors', phaseVectors, ...
                        'PlotType', plotType, ...
                        'FieldLabels', varLabels, ...
                        'PTicks', rowTickLocs, 'PTickLabels', rowTickLabels, ...
                        'PLabel', rowLabel, ...
                        'LineSpec', lineSpec, 'LineWidth', lineWidth, ...
                        'MarkerEdgeColor', markerEdgeColor, ...
                        'MarkerFaceColor', markerFaceColor, ...
                        'OutFolder', outFolder, ...
                        'FigTitles', figTitle, ...
                        'FigNames', figName, ...
                        'ClearFigure', clearFigure, 'AlwaysNew', alwaysNew, ...
                        'ShowFigure', showFigure, ...
                         otherArguments);
otherwise
    error('plotMode unrecognized!');
end

%% Outputs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Get the file bases
[~, fileBases, ~] = ...
    cellfun(@(x) fileparts(x), newPaths, 'UniformOutput', false);

% Create x tick labels
rowTickLabels = cellfun(@(x) strrep(x, '_', '\_'), fileBases, ...
                        'UniformOutput', false);

% Will not work for durations data
if isempty(rowTickLocs) && ~isempty(rowTickLabels)
    rowTickLocs = (1:numel(rowTickLabels))';
end

% Does not work if rowTickLocs is not also set
elseif isfield(myTable.Properties, 'RowTimes')
    % Convert time to minutes
    timeVec = minutes(myTable.Properties.RowTimes);

    % Convert to a cell array of character vectors
    rowTickLabels = convert_to_char(timeVec);

%                   - 'PlotSeparately': whether to plot each column separately
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
plotSeparatelyDefault = false;  % plot variables together by default
addParameter(iP, 'PlotSeparately', plotSeparatelyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
plotSeparately = iP.Results.PlotSeparately;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
