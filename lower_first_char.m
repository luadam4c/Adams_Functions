function strs = lower_first_char (strs, varargin)
%% Converts the first character of each string to lower case
% Usage: strs = lower_first_char (strs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       strs        - returned strings
%                   specified as a TODO
%
% Arguments:
%       strs        - strings
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%
% Used by:
%       cd/m3ha_plot_figure07.m

% File History:
% 2020-03-06 Created by Adam Lu
% TODO: Make 'UpperInstead' an optional argument

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
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['strs must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, strs, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% TODO: apply_to_each_str.m
if isempty(strs)
    return
elseif ischar(strs)
    strs = lower_first_char_helper(strs);
elseif iscell(strs)
    strs = cellfun(@lower_first_char_helper, strs, 'UniformOutput', false);
elseif isstring(strs)
    strs = arrayfun(@lower_first_char_helper, strs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = lower_first_char_helper (str)

% Convert to character array
if isstring(str)
    str = str{1};
    wasString = true;
else
    wasString = false;
end

% Lower the first character
if numel(str) >= 2
    str = [lower(str(1)), str(2:end)];
else
    str = lower(str);
end

% Convert back to string if necessary
if wasString
    str = string(str);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%