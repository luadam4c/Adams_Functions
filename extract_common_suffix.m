function suffix = extract_common_suffix (strs, varargin)
%% Extracts the common suffix of a cell array of strings
% Usage: suffix = extract_common_suffix (strs, varargin)
% Explanation:
%       TODO
% Example(s):
%       extract_common_suffix({'a_b_e_c_d.m', 'a_c_d.m', 'a_b_c_d.m'})
%
% Outputs:
%       suffix      - the common suffix
%                   specified as a character vector
% Arguments:
%       strs        - strings
%                   must be a character vector or a string vector
%                       or a cell array of character vectors
%       varargin    - 'Delimiter': delimiter used
%                   must be a string scalar or a character vector
%                   default == '_'
%                   - 'KeepDelimiter': whether to keep the preceding delimiter
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_common_prefix.m
%
% Used by:
%       cd/convert_units.m
%       cd/extract_distinct_fileparts.m
%       cd/plot_swd_histogram.m

% File History:
% 2018-12-26 Adapted from extract_common_prefix.m
% 2018-12-27 Added 'KeepDelimiter' as an optional argument


%% Hard-coded parameters

%% Default values for optional arguments
delimiterDefault = '_';         % use '_' as the delimiter by default
keepDelimiterDefault = false;   % don't keep the preceding delimiter by default

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
addParameter(iP, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'KeepDelimiter', keepDelimiterDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, strs, varargin{:});
delimiter = iP.Results.Delimiter;
keepDelimiter = iP.Results.KeepDelimiter;

%% Do the job
% Use extract_common_prefix.m
suffix = extract_common_prefix(strs, 'Delimiter', delimiter, ...
                        'KeepDelimiter', keepDelimiter, 'SuffixInstead', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
