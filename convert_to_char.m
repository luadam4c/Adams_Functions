function strs = convert_to_char(numbers, varargin)
%% Converts other data types to character arrays
% Usage: strs = convert_to_char(numbers, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       strs        - strings
%                   specified as a character array 
%                       or a cell array of character arrays
% Arguments:
%       numbers     - numbers
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/create_labels_from_numbers.m

% File History:
% 2018-12-27 Created by Adam Lu
% TODO: Make a convert_to_string.m for string array outputs
%           that can take non-scalar arguments
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
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'numbers', ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'3d'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, numbers, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
if numel(numbers) > 1
    strs = arrayfun(@convert_to_char_helper, numbers, ...
                    'UniformOutput', false);
else
    strs = convert_to_char_helper(numbers);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = convert_to_char_helper(x)

if isnumeric(x)
    str = num2str(x);
elseif islogical(x)
    str = char(string(x));
elseif isdatetime(x)
    str = datestr(x);
elseif isduration(x)
    str = char(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%