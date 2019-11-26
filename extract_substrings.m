function subStrs = extract_substrings (strs, varargin)
%% Extracts substrings from strings
% Usage: subStrs = extract_substrings (strs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       strs = {'Many_A105034_later', 'Mary_B203491_now'};
%       extract_substrings(strs, 'RegExp', '[A-Z][0-9]{6}')
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
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_plot_figure02.m

% File History:
% 2019-11-25 Created by Adam Lu
% 

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
    matchedCell = regexp(strs, regExp, 'match');

    % Extract the first match
    subStrs = cellfun(@(x) x{1}, matchedCell, ...
                                'UniformOutput', false);
end

%% Output results
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%