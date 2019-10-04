function error_unrecognized (stringVar, string, funcName)
%% Throws an error for unrecognized string
% Usage: error_unrecognized (stringVar, string, funcName)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Arguments:    
%       stringVar   - variable name for the unrecognized string
%                   must be a string scalar or a character vector
%       string      - the unrecognized string
%                   must be a string scalar or a character vector
%       funcName    - the function name
%                   must be a string scalar or a character vector
%
% Used by:
%       cd/compute_combined_trace.m
%       cd/compute_weighted_average.m
%       cd/extract_subvectors.m

% File History:
% 2018-10-26 Created by Adam Lu
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'stringVar', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'string', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'funcName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, stringVar, string, funcName);

%% Do the job
error(['The value ''%s'' of the string variable ''%s'' is unrecognized!\n', ...
        'Please type ''help %s'' for usage'], stringVar, string, funcName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%