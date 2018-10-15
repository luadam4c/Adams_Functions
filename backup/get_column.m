function column = get_column (T, columnName)
%% Extract the column with a specific variable name from a table
% Usage: column = get_column (T, columnName)
% Outputs:
%       column      - the column as a numeric array or a cell array
%                   specified as a TODO
% Arguments:
%       T           - table to get column from
%                   must be a table
%       columnName  - name of table column
%                   must be a string scalar or a character vector
%
% Used by:    
%       cd/plot_all_abfs.m

% File History:
% 2018-10-03 Created by Adam Lu
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'table', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addRequired(iP, 'columnName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, table, columnName);

%% Do the job
% Slice the table to get just that column
columnTable = T(:, columnName);

% Convert to the appropriate data type
% TODO
column = table2cell(columnTable);
% column = table2array(T(:, columnName));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%