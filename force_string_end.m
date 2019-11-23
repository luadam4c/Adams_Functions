function newStr = force_string_end (oldStr, subStr, varargin)
%% Force the string to end with a certain substring
% Usage: newStr = force_string_end (oldStr, subStr, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       force_string_end('dog', '/')
%       force_string_end("dog", "_")
%       force_string_end("dog_", "_")
%       force_string_end("dog", '_')
%       force_string_end("dog_", '_')
%       force_string_end(pwd, '.')
%       force_string_end("", '!', 'OnlyIfNonempty', true)
%       force_string_end("", "_", 'OnlyIfNonempty', true)
%       prefix = force_string_end(prefix, "_", 'OnlyIfNonempty', true)
%
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
%                   - 'StringStartInstead': whether to prepend substring
%                                               instead
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Used by:
%       cd/archive_dependent_scripts.m
%       cd/combine_swd_sheets.m
%       cd/create_labels_from_numbers.m
%       cd/create_simulation_output_filenames.m
%       cd/force_string_start.m
%       cd/run_neuron.m
%       cd/match_format_vector_sets.m
%       cd/m3ha_neuron_create_simulation_params.m
%       cd/m3ha_neuron_create_TC_commands.m
%       cd/m3ha_import_raw_traces.m
%       cd/parse_spike2_mat.m
%       cd/plot_traces_spike2_mat.m
%       cd/save_all_zooms.m
%       cd/test_var_difference.m

% File History:
% 2018-10-21 Created by Adam Lu
% 2019-01-01 Now allows oldStr and subStr to be cell arrays
% 2019-06-03 Now escapes the metacharacter .
% 2019-09-06 Added 'StringStartInstead' as an optional argument
% TODO: Escape all metacharacters for regexp
% TODO: Deal with the case when substr is more than one character
% TODO: force_string_start.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
onlyIfNonemptyDefault = false;      % append even if empty by default
stringStartInsteadDefault = false;  % force_string_end by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

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
addParameter(iP, 'StringStartInstead', stringStartInsteadDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, oldStr, subStr, varargin{:});
onlyIfNonempty = iP.Results.OnlyIfNonempty;
stringStartInstead = iP.Results.StringStartInstead;

%% Preparation
[oldStr, subStr] = match_format_vector_sets(oldStr, subStr, ...
                    'TreatCellStrAsArray', false, ...
                    'TreatCellAsArray', false);

%% Do the job
if iscell(oldStr)
    newStr = cellfun(@(x, y) force_string_end_helper(x, y, ...
                                    onlyIfNonempty, stringStartInstead), ...
                    oldStr, subStr, 'UniformOutput', false);
else
    newStr = force_string_end_helper(oldStr, subStr, ...
                                    onlyIfNonempty, stringStartInstead);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newStr = force_string_end_helper (oldStr, subStr, ...
                                            onlyIfNonempty, stringStartInstead)

% Return original string if empty and requested so
%   TODO: isempty_custom.m to consider "" as empty
if onlyIfNonempty && any(strcmp(oldStr, {'', ""}))
    newStr = oldStr;
    return
end

% Excape special characters for regexp
subStrRegExp = replace(subStr, '.', '\.');

% Form the regular expression to match
if stringStartInstead
    % Substr must occur in the beginning
    regExp = strcat('^', subStrRegExp);
else
    % Substr must occur at the end
    regExp = strcat(subStrRegExp, '$');
end

% Look for the substring at the end of the old string
startIndex = regexp(oldStr, regExp, 'ONCE');

% If not found, append the substring to the old string
if isempty(startIndex)
    if stringStartInstead
        % Prepend instead
        newStr = [subStr, oldStr];
    else
        newStr = [oldStr, subStr];
    end
else
    newStr = oldStr;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%   Note: One must use == '' if oldStr is a string type (in double quotes)

strcmp(oldStr, '')

% Note: the following will not work if a space is in between
newStr = strcat(subStr, oldStr);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
