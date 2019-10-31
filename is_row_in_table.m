function isInTable = is_row_in_table (rowName, table, varargin)
%% Returns whether a row name is an existing row in a table
% Usage: isInTable = is_row_in_table (rowName, table, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       isInTable   - whether the variable is an existing row in the table
%                   specified as a logical scalar
% Arguments:
%       rowName     - row name to look for
%                   must be a string scalar or a character vector
%       table       - table to look in
%                   must be a table
%       varargin    - Any other parameter-value pair for is_var_in_table()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/is_var_in_table.m
%
% Used by:
%       cd/m3ha_network_update_dependent_params.m

% File History:
% 2019-10-31 Created by Adam Lu
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'rowName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'table', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Read from the Input Parser
parse(iP, rowName, table, varargin{:});

% Keep unmatched arguments for the is_var_in_table() function
otherArguments = iP.Unmatched;

%% Do the job
isInTable = is_var_in_table(rowName, table, 'RowInstead', true, otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
