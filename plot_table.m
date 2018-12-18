function h = plot_table (table, varargin)
%% Plots all variables in a table as tuning curves
% Usage: h = plot_table (table, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       h           - handle to created figure
%                   specified as a figure handle
% Arguments:
%       table       - a table with variables to plot
%                   must be a table
%       varargin    - 'VariableNames': TODO: Description of VariableNames
%                   must be a TODO
%                   default == TODO
%                   - 'XLabel': label for the parameter
%                   must be a string scalar or a character vector
%                   default == none ('suppress')
%                   - 'OutFolder': output folder if FigNames not set
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - Any other parameter-value pair for the plot() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/plot_struct.m
%
% Used by:
%       cd/plot_protocols.m

% File History:
% 2018-12-18 Moved from plot_protocols.m
% 2018-12-18 Now uses iP.KeepUnmatched
% 

%% Hard-coded parameters

%% Default values for optional arguments
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
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'VariableNames', variableNamesDefault, ...
    @(x) assert(iscellstr(x) || isstring(x), ...
                ['VariableNames must be a cell array of character arrays ', ...
                'or a string array!']));
addParameter(iP, 'XLabel', xLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, table, varargin{:});
varToPlot = iP.Results.VariableNames;
xLabel = iP.Results.XLabel;
outFolder = iP.Results.OutFolder;

% Keep unmatched arguments for the line() function
otherArguments = iP.Unmatched;

%% Preparation
% Check if output directory exists
check_dir(outFolder);

% Restrict to variables to plot
if ~isempty(varToPlot)
    table = table(:, varToPlot);
end

% Decide on xTickLabels
if iscell(table.Properties.RowNames)
    % Get the file names
    fileNames = table.Properties.RowNames;

    % Get the file bases
    [~, fileBases, ~] = ...
        cellfun(@(x) fileparts(x), fileNames, 'UniformOutput', false);

    % Create x tick labels
    xTickLabels = cellfun(@(x) strrep(x, '_', '\_'), fileBases, ...
                            'UniformOutput', false);
else
    xTickLabels = {};
end

%% Do the job
% Convert to a structure array
fileStruct = table2struct(table);

% Plot fields
h = plot_struct(fileStruct, 'OutFolder', outFolder, ...
            'XTickLabels', xTickLabels, 'XLabel', xLabel, otherArguments);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%