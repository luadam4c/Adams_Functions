function subStrs = extract_substrings (strs, varargin)
%% Extracts substring(s) from strings
% Usage: subStrs = extract_substrings (strs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       extract_substrings('test')
%       extract_substrings('test01234', 'RegExp', '[\d]*')
%       extract_substrings('98765test01234', 'RegExp', '[\d]*')
%       strs = {'Many_A105034_later', 'Mary_B203491_now'};
%       extract_substrings(strs, 'RegExp', '[A-Z][0-9]{6}')
%       extract_substrings({'test', 'test23', '45test'}, 'RegExp', '[\d]*')
%
% Outputs:
%       subStrs     - substrings extracted
%                   specified as a character vector or a string vector
%                       or a cell array of character vectors
%
% Arguments:
%       strs        - strings to extract from
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%       varargin    - 'RegExp': regular expression to match
%                   must be a character vector or a string vector
%                       or a cell array of character vectors
%                   default == none
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/extract_common_prefix.m
%       cd/m3ha_extract_cell_name.m
%       cd/m3ha_extract_iteration_string.m
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2019-11-25 Created by Adam Lu
% 2019-11-28 Fixed the case when no substrings are found

%% Hard-coded parameters

%% Default values for optional arguments
regExpDefault = '';

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
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['strs5 must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'RegExp', regExpDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['regExp must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));

% Read from the Input Parser
parse(iP, strs, varargin{:});
regExp = iP.Results.RegExp;

% Check relationships between arguments
% TODO

%% Preparation
% TODO

%% Do the job
if ~isempty(regExp)
    % Match the regular expression
    matchedSubStrs = regexp(strs, regExp, 'match');

    % Extract the first match
    if iscellstr(matchedSubStrs)
        subStrs = extract_first_substr(matchedSubStrs);
    else
        subStrs = cellfun(@extract_first_substr, matchedSubStrs, ...
                            'UniformOutput', false);
    end
else
    subStrs = strs;
end

%% Output results
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subStr = extract_first_substr(matchedSubStrs)

if isempty(matchedSubStrs)
    subStr = '';
else
    subStr = matchedSubStrs{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%