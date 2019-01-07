function helpText = print_help (functionName, varargin)
%% Prints and returns the documentation for a specific function
% Usage: helpText = print_help (functionName, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       helpText    - documentation text
%                   specified as a character array
% Arguments:
%       functionName    - the function name to look for documentation for
%                       must be a string scalar or a character vector
%       varargin    - 'ToPrint': whether to actually print the string
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/create_error_for_nargin.m

% File History:
% 2018-12-17 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
toPrintDefault = true;                  % default: print the string

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
addRequired(iP, 'functionName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ToPrint', toPrintDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, functionName, varargin{:});
toPrint = iP.Results.ToPrint;

%% Do the job
% Return the help text
helpText = evalc(sprintf('help %s', functionName));

% Print the help text
if toPrint
    disp(helpText)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%