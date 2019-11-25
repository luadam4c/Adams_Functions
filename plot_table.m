function figs = plot_table (table, varargin)
%% Plots all variables in a table as tuning curves
% Usage: figs = plot_table (table, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       figs        - figure handle(s) for the created figure(s)
%                   specified as a figure object handle column vector
% Arguments:
%       table       - a table with variables to plot
%                   must be a table
%       varargin    - 'PlotType': type of plot
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'tuning'    - circles
%                       'bar'       - horizontal bars
%                   default == 'tuning'
%                   - 'VariableNames': variable (column) names of the table
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == plot all variables
%                   - 'PlotSeparately': whether to plot each column separately
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'VarLabels': variable labels
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == use distinct parts of variable names
%                   - 'DistinctParts': whether to extract distinct parts
%                                       or variable labels
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
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
%                   - 'PLabel': label for the parameter
%                   must be a string scalar or a character vector
%                   default == none ('suppress')
%                   - 'PTicks': x tick values for the parameter values
%                   must be a numeric vector
%                   default == 1:numel(pTickLabels)
%                   - 'PTickLabels': x tick labels
%                   must be a cell array of character vectors/strings
%                   default == row names or times if provided
%                   - 'OutFolder': output folder if FigNames not set
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - Any other parameter-value pair for the plot_struct() 
%                       or the plot_tuning_curve() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_common_directory.m
%       cd/extract_common_prefix.m
%       cd/extract_fileparts.m
%       cd/plot_struct.m
%       cd/plot_tuning_curve.m
%       cd/struct2arglist.m
%       ~/Downloaded_Functions/rgb.m
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
% 2019-08-07 Added 'PTickLabels' as an optional argument
% 2019-08-07 Added 'PTicks' as an optional argument
% TODO: Return handles to plots
% TODO: Pass in figNames or figNumbers when plotting separately
% 

%% Hard-coded parameters
validPlotTypes = {'tuning', 'bar'};
lineSpecSeparate = 'o';
lineWidthSeparate = 1;
markerEdgeColorSeparate = rgb('DarkOrchid');
markerFaceColorSeparate = rgb('LightSkyBlue');

lineSpecTogether = '-';
lineWidthTogether = 1;
markerEdgeColorTogether = rgb('DarkOrchid');
markerFaceColorTogether = rgb('LightSkyBlue');

%% Default values for optional arguments
plotTypeDefault = 'tuning';
lineSpecDefault = '';
lineWidthDefault = [];
markerEdgeColorDefault = [];
markerFaceColorDefault = [];
variableNamesDefault = {};      % plot all variables by default
plotSeparatelyDefault = false;  % plot variables together by default
varLabelsDefault = {};          % set later
distinctPartsDefault = true;    % extract distinct parts of variable names
                                %   by default
phaseVariablesDefault = {};     % no phases by default
delimiterDefault = '_';         % use '_' as delimiter by default
readoutLabelDefault = '';       % set later
tableLabelDefault = '';         % set later
pLabelDefault = 'suppress';     % No x label by default
pTicksDefault = [];
pTickLabelsDefault = {};
outFolderDefault = pwd;
figNameDefault = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% If not compiled, add directories to search path for required functions
if ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for rgb.m
    addpath_custom(fullfile(functionsDirectory, 'Downloaded_Functions'));
end

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
addRequired(iP, 'table', ...
    @(x) validateattributes(x, {'table', 'timetable'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PlotType', plotTypeDefault, ...
    @(x) any(validatestring(x, validPlotTypes)));
addParameter(iP, 'LineSpec', lineSpecDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'LineWidth', lineWidthDefault);
addParameter(iP, 'MarkerEdgeColor', markerEdgeColorDefault);
addParameter(iP, 'MarkerFaceColor', markerFaceColorDefault);
addParameter(iP, 'VariableNames', variableNamesDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['VariableNames must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'PlotSeparately', plotSeparatelyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'VarLabels', varLabelsDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['VarLabels must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'DistinctParts', distinctPartsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PhaseVariables', phaseVariablesDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['VariableNames must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ReadoutLabel', readoutLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TableLabel', tableLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PLabel', pLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PTicks', pTicksDefault, ...
    @(x) isempty(x) || isnumericvector(x));
addParameter(iP, 'PTickLabels', pTickLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, table, varargin{:});
plotType = validatestring(iP.Results.PlotType, validPlotTypes);
lineSpec = iP.Results.LineSpec;
lineWidth = iP.Results.LineWidth;
markerEdgeColor = iP.Results.MarkerEdgeColor;
markerFaceColor = iP.Results.MarkerFaceColor;
varToPlot = iP.Results.VariableNames;
plotSeparately = iP.Results.PlotSeparately;
varLabels = iP.Results.VarLabels;
distinctParts = iP.Results.DistinctParts;
phaseVariables = iP.Results.PhaseVariables;
delimiter = iP.Results.Delimiter;
readoutLabel = iP.Results.ReadoutLabel;
tableLabel = iP.Results.TableLabel;
pLabel = iP.Results.PLabel;
pTicks = iP.Results.PTicks;
pTickLabels = iP.Results.PTickLabels;
outFolder = iP.Results.OutFolder;
figName = iP.Results.FigName;

% Keep unmatched arguments for the line() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Check if output directory exists
check_dir(outFolder);

% Restrict to variables to plot or extract the variable names
if ~isempty(varToPlot)
    tableToPlot = table(:, varToPlot);
else
    tableToPlot = table;
    varToPlot = table.Properties.VariableNames;
end

% If provided, make sure there are an equal number of phase variables
%   and extract the phase vectors for each variable
if ~isempty(phaseVariables)
    % Force as a column cell array
    phaseVariables = force_column_cell(phaseVariables);

    % Match the number of phase variables to the number of variables to plot
    phaseVariables = match_row_count(phaseVariables, numel(varToPlot));

    % Extract the phase vectors
    phaseVectors = cellfun(@(x) table{:, x}, phaseVariables, ...
                            'UniformOutput', false);
else
    phaseVectors = {};
end

% Extract distinct parts if requested
if distinctParts
    varLabels = extract_fileparts(varToPlot, 'distinct', 'Delimiter', delimiter);
else
    varLabels = varToPlot;
end

% Decide on table label
if isempty(tableLabel)
    % First try to extract a common prefix from the variables to plot
    tableLabel = extract_common_prefix(varToPlot, 'Delimiter', delimiter);

    % If no such prefix exists, use the table variable name
    if isempty(tableLabel)
        tableLabel = inputname(1);
    end
end

% Decide on pTickLabels
if isempty(pTickLabels)
    if isfield(table.Properties, 'RowNames') && ...
            iscell(table.Properties.RowNames)
        % Get the row names
        rowNames = table.Properties.RowNames;

        % If all row names are file names, process them
        %   Otherwise, just use the row names as the x tick labels
        if all(isfile(rowNames))
            % Extract the distinct file bases
            pTickLabels = extract_fileparts(rowNames, 'distinct');

            % Replace all instances of '_' with '\_'
            pTickLabels = replace(pTickLabels, '_', '\_');
        else
            % Just use the row names
            pTickLabels = rowNames;
        end
    else
        % Use default x tick labels
        pTickLabels = {};
    end
end

% Decide on lineSpec
if isempty(lineSpec)
    if plotSeparately
        lineSpec = lineSpecSeparate;
    else
        lineSpec = lineSpecTogether;
    end
end

% Decide on lineWidth
if isempty(lineWidth)
    if plotSeparately
        lineWidth = lineWidthSeparate;
    else
        lineWidth = lineWidthTogether;
    end
end

% Decide on markerEdgeColor
if isempty(markerEdgeColor)
    if plotSeparately
        markerEdgeColor = markerEdgeColorSeparate;
    else
        markerEdgeColor = markerEdgeColorTogether;
    end
end

% Decide on markerFaceColor
if isempty(markerFaceColor)
    if plotSeparately
        markerFaceColor = markerFaceColorSeparate;
    else
        markerFaceColor = markerFaceColorTogether;
    end
end

%% Do the job
if plotSeparately
    % Convert to a structure array
    myStruct = table2struct(tableToPlot);

    % Plot fields
    figs = plot_struct(myStruct, 'OutFolder', outFolder, ...
                        'PhaseVectors', phaseVectors, ...
                        'PlotType', plotType, ...
                        'FieldLabels', varLabels, ...
                        'PTicks', pTicks, 'PTickLabels', pTickLabels, ...
                        'PLabel', pLabel, ...
                        'LineSpec', lineSpec, 'LineWidth', lineWidth, ...
                        'MarkerEdgeColor', markerEdgeColor, ...
                        'MarkerFaceColor', markerFaceColor, ...
                         otherArguments{:});
else
    % Create a figure name if empty
    if isempty(figName)
        figName = fullfile(outFolder, [tableLabel, '.png']);
        figName = fullfile(outFolder, [tableLabel, '.epsc']);
    end

    % Convert to an array
    if istimetable(tableToPlot)
        % Extract variables
        myArray = tableToPlot.Variables;
    else
        % Use table2array
        myArray = table2array(tableToPlot);
    end

    % Decide on x values
    if istimetable(tableToPlot)
        % Extract time
        xValues = tableToPlot.Properties.RowTimes;
    else
        % Count rows
        nRows = height(tableToPlot);

        % Use row numbers
        xValues = transpose(1:nRows);
    end

    % Decide on readout label
    if isempty(readoutLabel)
        readoutLabel = replace(tableLabel, '_', ' ');
    end

    % Decide on figure title
    if istimetable(tableToPlot)
        figTitle = replace([tableLabel, ' over time'], '_', ' ');
    else
        figTitle = '';
    end

    % Clear the current figure
    clf;

    % Plot a tuning curve
    switch plotType
        case 'tuning'
            handles = plot_tuning_curve(xValues, myArray, 'FigName', figName, ...
                            'PhaseVectors', phaseVectors, ...
                            'PTicks', pTicks, 'PTickLabels', pTickLabels, ...
                            'PLabel', pLabel, ...
                            'ReadoutLabel', readoutLabel, ...
                            'ColumnLabels', varLabels, ...
                            'FigTitle', figTitle, ...
                            'LineSpec', lineSpec, 'LineWidth', lineWidth, ...
                            'MarkerEdgeColor', markerEdgeColor, ...
                            'MarkerFaceColor', markerFaceColor, ...
                            otherArguments{:});
            figs = handles.fig;
        case 'bar'
            % TODO
        otherwise
            error('plotType unrecognized!')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Get the file bases
[~, fileBases, ~] = ...
    cellfun(@(x) fileparts(x), newPaths, 'UniformOutput', false);

% Create x tick labels
pTickLabels = cellfun(@(x) strrep(x, '_', '\_'), fileBases, ...
                        'UniformOutput', false);

% Will not work for durations data
if isempty(pTicks) && ~isempty(pTickLabels)
    pTicks = (1:numel(pTickLabels))';
end

% Does not work if pTicks is not also set
elseif isfield(table.Properties, 'RowTimes')
    % Convert time to minutes
    timeVec = minutes(table.Properties.RowTimes);

    % Convert to a cell array of character vectors
    pTickLabels = convert_to_char(timeVec);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
