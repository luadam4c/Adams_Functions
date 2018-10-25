function newStr = force_string_end (oldStr, subStr, varargin)
%% Force the string to end with a certain substring
% Usage: newStr = force_string_end (oldStr, subStr, varargin)
% Explanation:
%       TODO
% Example(s):
%       force_string_end('dog', '/')
%       force_string_end("dog", "_")
%       force_string_end("dog", '_')
%       force_string_end("", '!', 'OnlyIfNonempty', true)
%       force_string_end("", "_", 'OnlyIfNonempty', true)
%       prefix = force_string_end(prefix, "_", 'OnlyIfNonempty', true)
% Outputs:
%       newStr      - resulting string
%                   specified as a string scalar or a character vector
% Arguments:    
%       oldStr      - original string
%                   must be a string scalar or a character vector
%       subStr      - substring to end with
%                   must be a string scalar or a character vector
%       varargin    - 'OnlyIfNonempty': whether to append substring
%                                       only if original string is non-empty
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Used by:    
%       cd/create_simulation_output_filenames.m
%       cd/run_neuron.m
%       cd/m3ha_create_simulation_params.m
%       cd/m3ha_create_single_neuron_commands.m

% File History:
% 2018-10-21 Created by Adam Lu
% TODO: Deal with the case when substr is more than one character
% 

%% Hard-coded parameters

%% Default values for optional arguments
onlyIfNonemptyDefault = false;  % append even if empty by default

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

% Add required inputs to the Input Parser
addRequired(iP, 'oldStr', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'subStr', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OnlyIfNonempty', onlyIfNonemptyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, oldStr, subStr, varargin{:});
onlyIfNonempty = iP.Results.OnlyIfNonempty;

%% Do the job
% Return original string if empty and requested so
if onlyIfNonempty && strcmp(oldStr, '')
    newStr = oldStr;
    return
end

% Look for the substring at the end of the old string
startIndex = regexp(oldStr, strcat(subStr, '$'));

% If not found, append the substring to the old string
if isempty(startIndex)
    newStr = strcat(oldStr, subStr);
else
    newStr = oldStr;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%   Note: One must use == '' if oldStr is a string type (in double quotes)

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%