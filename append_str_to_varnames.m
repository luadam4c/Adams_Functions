function myTable = append_str_to_varnames (myTable, str, varargin)
%% Appends a string to all variable names in a table
% Usage: myTable = append_str_to_varnames (myTable, str, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       myTable     - TODO: Description of myTable
%                   specified as a TODO
%
% Arguments:
%       myTable     - TODO: Description of myTable
%                   must be a TODO
%       str         - TODO: Description of myTable
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2020-08-19 Created by Adam Lu
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
addRequired(iP, 'myTable', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addRequired(iP, 'str', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, myTable, str, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Extract originial variable names
varNamesOrig = myTable.Properties.VariableNames;

% Append the string
varNamesNew = strcat(varNamesOrig, str);

% Update variable names
myTable.Properties.VariableNames = varNamesNew;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%