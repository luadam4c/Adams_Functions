function parts = extract_fileparts (paths, partType, varargin)
%% Extracts directories, bases, extensions, distinct parts or the common directory from file paths, treating any path without an extension as a directory
% Usage: parts = extract_fileparts (paths, partType, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       parts       - parts extracted
%                   specified as a character array 
%                       or a cell array of character arrays
% Arguments:
%       paths       - file paths
%                   must be a character vector or a string vector
%                       or a cell array of character vectors
%       partType    - type of the file part to extract
%                   must be an unambiguous, case-insensitive match to one of:
%                       'commondirectory' - common directory across file(s)
%                       'commonsuffix'    - common suffix across file(s)
%                       'distinct'  - distinct parts across file(s)
%                       'directory' - directory containing the file(s)
%                       'base'      - file base name without the extension
%                       'extension' - file extension including the leading '.'
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_common_directory.m
%
% Used by:
%       cd/extract_common_directory.m
%       cd/plot_table.m

% File History:
% 2018-12-18 Created by Adam Lu
% TODO: Make the first argument accept a files structure array too
% 

%% Hard-coded parameters
validPartTypes = {'commondirectory', 'commonsuffix', 'distinct', ...
                    'directory', 'base', 'extension'};

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'paths', ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['paths must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addRequired(iP, 'partType', ...
    @(x) any(validatestring(x, validPartTypes)));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, paths, partType, varargin{:});
% param1 = iP.Results.param1;

% Validate partType
partType = validatestring(partType, validPartTypes);

%% Do the job
switch partType
case {'directory', 'base', 'extension'}
    parts = extract_simple_fileparts(paths, partType);
case 'commondirectory'
    parts = extract_common_directory(paths, varargin{:});
case 'distinct'
    parts = extract_distinct_parts(paths);
otherwise
    error('partType unrecognized!!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fileDir, fileBase, fileExtension] = fileparts_custom (filePath)
%% Same as fileparts but treats anything without an extension as a directory

% Use the original fileparts to get the extension
[fileDirTentative, fileBaseTentative, fileExtension] = fileparts(filePath);

% Check if there is an extension
if isempty(fileExtension)
    % If there is no extension, it is a directory
    fileDir = fullfile(fileDirTentative, fileBaseTentative);
    fileBase = '';
else
    % If there is an extension, it is a file
    fileDir = fileDirTentative;
    fileBase = fileBaseTentative;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parts = extract_simple_fileparts(paths, partType)
% Separate the paths with fileparts_custom()

if iscell(paths)
    [fileDirs, fileBases, fileExtensions] = ...
        cellfun(@(x) fileparts_custom(x), paths, 'UniformOutput', false);
else
    [fileDirs, fileBases, fileExtensions] = fileparts_custom(paths);
end

% Return results
switch partType
    case 'directory'
        parts = fileDirs;
    case 'base'
        parts = fileBases;
    case 'extension'
        parts = fileExtensions;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function distinctParts = extract_distinct_parts(paths)
% Extract the distinct parts from a set of file paths

% Extract the common parent directory
commonParent = extract_common_directory(paths);

% Extract everything after the common parent directory
relativePaths = extractAfter(paths, commonParent);

% Extract the file extensions
fileExt = extract_fileparts(relativePaths, 'extension');

% Remove the file extensions
distinctParts = extractBefore(relativePaths, fileExt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
