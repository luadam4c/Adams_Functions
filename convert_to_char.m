function strs = convert_to_char(data, varargin)
%% Converts other data types to character arrays or a cell array of character arrays
% Usage: strs = convert_to_char(data, varargin)
% Explanation:
%       TODO
% Example(s):
%       convert_to_char(1:10)
%       convert_to_char({1:10, 3:5})
%       convert_to_char({'dog', 'cat'})
%       convert_to_char(["dog", "cat"])
%       convert_to_char({"dog", "cat"})
%       convert_to_char({{'dog', 'cat'}, "fly"})
%       convert_to_char({{'dog', 'cat'}, "fly"}, 'SingleOutput', true)
% Outputs:
%       strs        - strings
%                   specified as a character array 
%                       or a cell array of character arrays
% Arguments:
%       data        - data
%       varargin    - 'SingleOutput': whether to output a single character array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Delimiter': used to delimit separate entries
%                   must be a character array
%                   default == '_'
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/create_labels_from_numbers.m
%       cd/test_difference.m

% File History:
% 2018-12-27 Created by Adam Lu
% 2019-01-11 Now accepts any data type
% 2019-01-11 Added 'SingleOutput' and 'Delimiter' as optional arguments
% TODO: Make a convert_to_string.m for string array outputs
%           that can take non-scalar arguments
% 

%% Hard-coded parameters

%% Default values for optional arguments
singleOutputDefault = false;    % accept cell array outputs by default
delimiterDefault = '_';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'data');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SingleOutput', singleOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, data, varargin{:});
singleOutput = iP.Results.SingleOutput;
delimiter = iP.Results.Delimiter;               % Examples: ',' '/'

%% Preparation
% Do nothing if already a cell array of character arrays
%   or a character array
if ischar(data) || ~singleOutput && iscellstr(data)
    strs = data;
    return
end

%% Do the job
if numel(data) > 1
    if singleOutput && (iscellstr(data) || isstring(data))
        % TODO: Use print_cellstr.m?
        strs = char(strjoin(data, delimiter));
    elseif iscell(data)
        strs = cellfun(@(x) convert_to_char_helper(x, singleOutput, ...
                                                    delimiter), ...
                        data, 'UniformOutput', false);
    else
        strs = arrayfun(@(x) convert_to_char_helper(x, singleOutput, ...
                                                    delimiter), ...
                        data, 'UniformOutput', false);
    end
else
    strs = convert_to_char_helper(data);
end

if singleOutput && ~ischar(strs)
    strs = convert_to_char(strs, 'SingleOutput', singleOutput);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = convert_to_char_helper(x, singleOutput, delimiter)

% If there is more than one element, apply the function to this element
if numel(x) > 1
    str = convert_to_char(x, 'SingleOutput', singleOutput, ...
                            'Delimiter', delimiter);
    return
end

% Otherwise, convert based on data type
if isnumeric(x)
    str = num2str(x);
elseif islogical(x)
    str = char(string(x));
elseif isdatetime(x)
    str = datestr(x);
elseif isduration(x) || iscellstr(x) || isstring(x)
    str = char(x);
else
    str = char(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
addRequired(iP, 'data', ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'3d'}));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%