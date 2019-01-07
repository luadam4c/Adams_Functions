function distinctParts = extract_distinct_fileparts (paths, varargin)
%% Extracts distinct file parts (removes common parent directory and common suffix)
% Usage: distinctParts = extract_distinct_fileparts (paths, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       distinctParts   - distinct parts of the file paths
%                       specified as a TODO
% Arguments:
%       paths       - file paths
%                   must be a character vector or a string vector
%                       or a cell array of character vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_common_directory.m
%       cd/extract_fileparts.m
%       cd/extract_common_suffix.m
%
% Used by:
%       cd/extract_fileparts.m
%       cd/plot_swd_histogram.m

% File History:
% 2018-12-27 Moved from extract_fileparts.m
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

% Add required inputs to the Input Parser
addRequired(iP, 'paths', ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['paths must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, paths, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Extract the common parent directory
commonParent = extract_common_directory(paths, 'KeepFileSep', true);

% Extract everything after the common parent directory
relativePaths = extractAfter(paths, commonParent);

% Extract the file extensions
fileExt = extract_fileparts(relativePaths, 'extension');

% Remove the file extensions
distinctParts = extractBefore(relativePaths, fileExt);

% Extract the common suffix
commonSuffix = extract_common_suffix(distinctParts, 'KeepDelimiter', true);

% Remove the common suffix
distinctParts = extractBefore(distinctParts, commonSuffix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%