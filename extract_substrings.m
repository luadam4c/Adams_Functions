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
%       extract_substrings({'test45', 'test23'}, 'ForceSingleOutput', true)
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
%                   - 'FromBaseName': whether to restrict to the base name
%                                       (remove everything before last filesep
%                                        and after last .)
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ForceSingleOutput': whether to force as a single output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_common_prefix.m
%
% Used by:
%       cd/compute_activation_profile.m
%       cd/extract_common_prefix.m
%       cd/m3ha_extract_cell_name.m
%       cd/m3ha_extract_iteration_string.m
%       cd/m3ha_extract_sweep_name.m
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2019-11-25 Created by Adam Lu
% 2019-11-28 Fixed the case when no substrings are found

%% Hard-coded parameters

%% Default values for optional arguments
regExpDefault = '';
fromBaseNameDefault = false;
forceSingleOutputDefault = false;

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
addParameter(iP, 'FromBaseName', fromBaseNameDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ForceSingleOutput', forceSingleOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, strs, varargin{:});
regExp = iP.Results.RegExp;
fromBaseName = iP.Results.FromBaseName;
forceSingleOutput = iP.Results.ForceSingleOutput;

%% Do the job
% Restrict to base name if requested
if fromBaseName
    strs = extract_fileparts(strs, 'base');
end

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

% Reduce to one cell name if requested
if forceSingleOutput
    % Extract the common prefix
    subStrs = extract_common_prefix(subStrs);
    
    % If doesn't exist, change it to 'mixed'
    if isempty(subStrs)
        fprintf('Warning: Common substring doesn''t exist!\n');
        subStrs = 'mixed';
    end
end

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