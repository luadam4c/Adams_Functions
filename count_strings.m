function nStrings = count_strings (strs, varargin)
%% Count the number of strings
% Usage: nStrings = count_strings (strs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       count_strings({'apple', 'banana'})
%       count_strings(["apple", "banana"])
%       count_strings('apple')
%       count_strings({{'apple', 'banana'}, ["Mark", "Katie", "Matt"]})
%
% Outputs:
%       nStrings    - number of strings
%                   specified as a numeric scalar
%
% Arguments:
%       strs        - strings to count
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/plot_table_parallel.m
%       cd/read_lines_from_files.m

% File History:
% 2019-12-26 Created by Shin-Shin Nien
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'strs', ...
    @(x) assert(iscell(x) || ischar(x) || isstring(x), ...
        ['strs5 must be a character array or a string array ', ...
            'or cell array of them!']));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, strs, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
if iscell(strs) && ~iscellstr(strs)
   nStrings = cellfun(@count_strings, strs);
else
    if ischar(strs)
        nStrings = 1;
    else 
        nStrings = numel(strs);
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
