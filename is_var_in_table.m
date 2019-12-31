function isInTable = is_var_in_table (varName, table, varargin)
%% Returns whether a variable name is an existing column in a table
% Usage: isInTable = is_var_in_table (varName, table, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       T = readtable('outages.csv');
%       is_var_in_table('Cause', T)
%
% Outputs:
%       isInTable   - whether the variable is an existing column in the table
%                   specified as a logical scalar
%
% Arguments:
%       varName     - variable (column) name to look for
%                   must be a string scalar or a character vector
%       table       - table to look in
%                   must be a table
%       varargin    - 'RowInstead': find in rows instead
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for the ismatch() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/ismatch.m
%
% Used by:
%       cd/is_field.m
%       cd/is_row_in_table.m
%       cd/m3ha_network_change_params.m
%       cd/m3ha_select_sweeps.m
%       cd/plot_measures.m
%       /media/adamX/m3ha/optimizer4compgabab/singleneuronfitting59.m

% File History:
% 2019-10-31 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
rowInsteadDefault = false;

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
addRequired(iP, 'varName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'table', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'RowInstead', rowInsteadDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varName, table, varargin{:});
rowInstead = iP.Results.RowInstead;

% Keep unmatched arguments for the ismatch() function
otherArguments = iP.Unmatched;

%% Do the job
% Get all variable names
if rowInstead
    allVarNames = table.Properties.RowNames;
else
    allVarNames = table.Properties.VariableNames;
end

% Test whether varName is a match
isInTable = any(ismatch(allVarNames, varName, otherArguments));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
