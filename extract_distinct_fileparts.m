function distinctParts = extract_distinct_fileparts (paths, varargin)
%% Extracts distinct file parts (removes common parent directory, common prefix and common suffix)
% Usage: distinctParts = extract_distinct_fileparts (paths, varargin)
% Explanation:
%       TODO
% Example(s):
%       [~, paths] = all_files;
%       extract_distinct_fileparts(paths)
%       extract_distinct_fileparts({'a+b+c', 'a+d+c'}, 'Delimiter', '+')
%
% Outputs:
%       distinctParts   - distinct parts of the file paths
%                       specified as a TODO
% Arguments:
%       paths       - file paths
%                   must be a character vector or a string vector
%                       or a cell array of character vectors
%       varargin    - 'Delimiter': delimiter used for prefix and/or suffix
%                   must be a string scalar or a character vector
%                   default == '_'
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_common_directory.m
%       cd/extract_fileparts.m
%       cd/extract_common_suffix.m
%       cd/isemptycell.m
%
% Used by:
%       cd/extract_fileparts.m
%       cd/plot_swd_histogram.m

% File History:
% 2018-12-27 Moved from extract_fileparts.m
% 2019-01-23 Fixed the case where commonSuffix is empty
% 2019-03-17 Now removes common prefix as well
% 2019-03-17 Added 'Delimiter' as an optional argument
% 

%% Hard-coded parameters

%% Default values for optional arguments
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

% Add required inputs to the Input Parser
addRequired(iP, 'paths', ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['paths must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, paths, varargin{:});
delimiter = iP.Results.Delimiter;

%% Do the job
% Extract the common parent directory
commonParent = extract_common_directory(paths, 'KeepFileSep', true);

% Extract everything after the common parent directory
if ~isempty(commonParent)
    relativePaths = extractAfter(paths, commonParent);
else
    relativePaths = paths;
end

% Extract the file extensions
fileExt = extract_fileparts(relativePaths, 'extension');

% Remove the file extensions if not empty
if ~isempty(fileExt) && ~all(isemptycell(fileExt))
    distinctParts = extractBefore(relativePaths, fileExt);
else
    distinctParts = relativePaths;    
end

% Extract the common suffix
commonSuffix = extract_common_suffix(distinctParts, 'Delimiter', delimiter, ...
                                     'KeepDelimiter', true);

% Remove the common suffix if not empty
if ~isempty(commonSuffix)
    distinctParts = extractBefore(distinctParts, commonSuffix);
end

% Extract the common prefix
commonPrefix = extract_common_prefix(distinctParts, 'Delimiter', delimiter, ...
                                     'KeepDelimiter', true);

% Remove the common prefix if not empty
if ~isempty(commonPrefix)
    distinctParts = extractAfter(distinctParts, commonPrefix);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%