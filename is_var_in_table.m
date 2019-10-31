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
% Arguments:
%       varName     - variable (column) name to look for
%                   must be a string scalar or a character vector
%       table       - table to look in
%                   must be a table
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the ismatch() function
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/m3ha_network_change_params.m

% File History:
% 2019-10-31 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, varName, table, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the ismatch() function
otherArguments = iP.Unmatched;

%% Do the job
% Get all variable names of the table
allVarNames = table.Properties.VariableNames;

% Test whether varName is a match
isInTable = any(ismatch(allVarNames, varName, otherArguments));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%