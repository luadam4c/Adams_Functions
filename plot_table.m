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
%       varargin    - 'VariableNames': variable (column) names of the table
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == plot all variables
%                   - 'XLabel': label for the parameter
%                   must be a string scalar or a character vector
%                   default == none ('suppress')
%                   - 'OutFolder': output folder if FigNames not set
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - Any other parameter-value pair for the plot_struct() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_common_directory.m
%       cd/extract_fileparts.m
%       cd/plot_struct.m
%       cd/struct2arglist.m
%       ~/Downloaded_Functions/rgb.m
%
% Used by:
%       cd/plot_protocols.m

% File History:
% 2018-12-18 Moved from plot_protocols.m
% 2018-12-18 Now uses iP.KeepUnmatched
% 2018-12-18 Now uses extract_common_directory.m
% 2018-12-18 Now uses row names without processing if not file names
% 2019-03-17 Deal with timetables differently for RowNames
% TODO: PlotTogether (use plot_tuning_curve directly)
% TODO: Return handles to plots
% 

%% If not compiled, add directories to search path for required functions
if ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for rgb.m
    addpath_custom(fullfile(functionsDirectory, 'Downloaded_Functions'));
end

%% Hard-coded parameters
plotTogether = false;

%% Default values for optional arguments
lineSpecDefault = 'o';
lineWidthDefault = 1;
markerEdgeColorDefault = rgb('DarkOrchid');
markerFaceColorDefault = rgb('LightSkyBlue');
variableNamesDefault = {};  % plot all variables by default
xLabelDefault = 'suppress'; % No x label by default
outFolderDefault = pwd;

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
addRequired(iP, 'table', ...
    @(x) validateattributes(x, {'table', 'timetable'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'LineSpec', lineSpecDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'LineWidth', lineWidthDefault);
addParameter(iP, 'MarkerEdgeColor', markerEdgeColorDefault);
addParameter(iP, 'MarkerFaceColor', markerFaceColorDefault);
addParameter(iP, 'VariableNames', variableNamesDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['VariableNames must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'XLabel', xLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, table, varargin{:});
lineSpec = iP.Results.LineSpec;
lineWidth = iP.Results.LineWidth;
markerEdgeColor = iP.Results.MarkerEdgeColor;
markerFaceColor = iP.Results.MarkerFaceColor;
varToPlot = iP.Results.VariableNames;
xLabel = iP.Results.XLabel;
outFolder = iP.Results.OutFolder;

% Keep unmatched arguments for the line() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Check if output directory exists
check_dir(outFolder);

% Restrict to variables to plot
if ~isempty(varToPlot)
    table = table(:, varToPlot);
end

% Decide on xTickLabels
if isfield(table.Properties, 'RowNames') && iscell(table.Properties.RowNames)
    % Get the row names
    rowNames = table.Properties.RowNames;

    % If all row names are file names, process them
    %   Otherwise, just use the row names as the x tick labels
    if all(isfile(rowNames))
        % Extract the distinct file bases
        xTickLabels = extract_fileparts(rowNames, 'distinct');

        % Replace all instances of '_' with '\_'
        xTickLabels = replace(xTickLabels, '_', '\_');
    else
        % Just use the row names
        xTickLabels = rowNames;
    end
elseif isfield(table.Properties, 'RowTimes')
    xTickLabels = convert_to_char(table.Properties.RowTimes);
else
    % Use default x tick labels
    xTickLabels = {};
end

% Extract the variable names
varNames = table.Properties.VariableNames;

% TODO: extractDistinct

%% Do the job
if plotTogether
    % TODO: Extract
    % fileArray = table2array(table)

    % 
    % plot_tuning_curve(fileArray, 'ColumnLabels', varNames, otherArguments{:});
else
    % Convert to a structure array
    fileStruct = table2struct(table);

    % Plot fields
    figs = plot_struct(fileStruct, 'OutFolder', outFolder, ...
                        'FieldLabels', varNames, ...
                        'XTickLabels', xTickLabels, 'XLabel', xLabel, ...
                        'LineSpec', lineSpec, 'LineWidth', lineWidth, ...
                        'MarkerEdgeColor', markerEdgeColor, ...
                        'MarkerFaceColor', markerFaceColor, ...
                        otherArguments{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Get the file bases
[~, fileBases, ~] = ...
    cellfun(@(x) fileparts(x), newPaths, 'UniformOutput', false);

% Create x tick labels
xTickLabels = cellfun(@(x) strrep(x, '_', '\_'), fileBases, ...
                        'UniformOutput', false);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
