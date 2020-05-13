function combinedTable = combine_param_tables (paramPathsOrTables, varargin)
%% Combine parameter tables with a 'Value' column and row names as parameters
% Usage: combinedTable = combine_param_tables (paramPathsOrTables, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       combinedTable   - combined table of parameters
%                       specified as a table
%
% Arguments:
%       paramPathsOrTables  - parameter tables or file paths to them
%                               table must have a 'Value' column 
%                                   and row names
%                           must be a cell array of tables or a string array 
%                               or a cell array of character arrays
%       varargin    - 'NewRowNames': new row names
%                   must be a string array or a cell array of character arrays
%                   default == parameter path base names
%
% Requires:
%       cd/apply_over_cells.m
%       cd/array_fun.m
%       cd/create_error_for_nargin.m
%       cd/create_label_from_numbers.m
%       cd/extract_fileparts.m
%       cd/force_column_cell.m
%       cd/renamevars_custom.m
%       cd/transpose_table.m
%
% Used by:
%       cd/m3ha_network_analyze_spikes.m
%       cd/m3ha_simulate_population.m

% File History:
% 2020-04-09 Moved from in m3ha_network_analyze_spikes.m
% 

%% Hard-coded parameters

% TODO
saveFlag = false;

%% Default values for optional arguments
newRowNamesDefault = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'paramPathsOrTables');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'newRowNames', newRowNamesDefault);

% Read from the Input Parser
parse(iP, paramPathsOrTables, varargin{:});
newRowNames = iP.Results.newRowNames;

%% Preparation
% Force as a column cell array
paramPathsOrTables = force_column_cell(paramPathsOrTables);

% Parse the first argument
if iscellstr(paramPathsOrTables)
    % First argument is parameter paths
    paramPaths = paramPathsOrTables;

    % Load simulation parameters
    paramTables = array_fun(@(x) readtable(x, 'ReadRowNames', true), ...
                            paramPaths, 'UniformOutput', false);
else
    % First argument is parameter tables
    paramTables = paramPathsOrTables;

    % Count the number of tables
    nTables = numel(paramTables);

    % Create table paths
    paramPaths = create_label_from_numbers(1:nTables, 'Prefix', 'param_table_');
end

% Decide on new row names
if isempty(newRowNames)
    newRowNames = extract_fileparts(paramPaths, 'base');
end

%% Do the job
% Keep just the Value column in all tables
paramTables = array_fun(@(x) x(:, 'Value'), ...
                        paramTables, 'UniformOutput', false);

% Rename 'Value' by the condition string
paramTables = array_fun(@(x, y) renamevars_custom(x, 'Value', y), ...
                        paramTables, newRowNames, 'UniformOutput', false);

% Combine all tables
allParamsTable = apply_over_cells(@horzcat, paramTables);

% Initialize the oscillation table with simulation parameters
combinedTable = transpose_table(allParamsTable);

%% Save results
if saveFlag
    % TODO
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
