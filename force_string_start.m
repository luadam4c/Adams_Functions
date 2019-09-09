function newStr = force_string_start (oldStr, subStr, varargin)
%% Force the string to start with a certain substring
% Usage: newStr = force_string_start (oldStr, subStr, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       force_string_start('dog', '/')
%       force_string_start("dog", "_")
%       force_string_start("_dog", "_")
%       force_string_start("dog", '_')
%       force_string_start("_dog", '_')
%       force_string_start(pwd, '.')
%       force_string_start("", '!', 'OnlyIfNonempty', true)
%       force_string_start("", "_", 'OnlyIfNonempty', true)
%       suffix = force_string_start(suffix, "_", 'OnlyIfNonempty', true)
%
% Outputs:
%       newStr      - resulting string
%                   specified as a string scalar or a character vector
%
% Arguments:    
%       oldStr      - original string
%                   must be a string scalar or a character vector
%       subStr      - substring to end with
%                   must be a string scalar or a character vector
%       varargin    - 'OnlyIfNonempty': whether to append substring
%                                       only if original string is non-empty
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for force_string_end()
%
% Requires:
%       cd/force_string_end.m
%
% Used by:
%       cd/all_files.m
%       cd/create_labels_from_numbers.m

% File History:
% 2019-09-06 Modified from force_string_end.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
onlyIfNonemptyDefault = false;      % append even if empty by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'oldStr', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['oldStr must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addRequired(iP, 'subStr', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['oldStr must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OnlyIfNonempty', onlyIfNonemptyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, oldStr, subStr, varargin{:});
onlyIfNonempty = iP.Results.OnlyIfNonempty;

% Keep unmatched arguments for the force_string_end() function
otherArguments = iP.Unmatched;

%% Do the job
newStr = force_string_end(oldStr, subStr, 'StringStartInstead', true, ...
                        'OnlyIfNonempty', onlyIfNonempty, otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
