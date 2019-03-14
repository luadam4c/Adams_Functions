function parts = extract_fileparts (paths, partType, varargin)
%% Extracts directories, bases, extensions, distinct parts or the common directory from file paths, treating any path without an extension as a directory
% Usage: parts = extract_fileparts (paths, partType, varargin)
% Explanation:
%       TODO
% Example(s):
%       [~, paths] = all_files('Directory', pwd);
%       extract_fileparts(paths, 'commondirectory')
%       extract_fileparts(paths, 'commonprefix')
%       extract_fileparts(paths, 'commonsuffix')
%       extract_fileparts(paths, 'distinct')
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
%                       'commonprefix'    - common prefix across file(s)
%                       'commonsuffix'    - common suffix across file(s)
%                       'distinct'  - distinct parts across file(s)
%                       'directory' - directory containing the file(s)
%                       'base'      - file base name without the extension
%                       'extension' - file extension including the leading '.'
%       varargin    - 'Delimiter': delimiter used for file suffices
%                   must be a string scalar or a character vector
%                   default == '_'
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_common_directory.m
%       cd/extract_common_suffix.m
%       cd/extract_distinct_fileparts.m
%
% Used by:
%       cd/all_filebases.m
%       cd/extract_common_directory.m
%       cd/parse_all_multiunit.m
%       cd/plot_swd_histogram.m
%       cd/plot_table.m

% File History:
% 2018-12-18 Created by Adam Lu
% 2018-12-26 Added 'commonsuffix' as a part type
% 2018-12-27 Moved code to extract_distinct_fileparts.m
% 2019-03-14 Added 'commonprefix' as a part type
% TODO: Make the first argument accept a files structure array too
% 

%% Hard-coded parameters
validPartTypes = {'commondirectory', 'commonprefix', 'commonsuffix', ...
                    'distinct', 'directory', 'base', 'extension'};

%% Default values for optional arguments
delimiterDefault = '_';

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
addParameter(iP, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, paths, partType, varargin{:});
delimiter = iP.Results.Delimiter;

% Validate partType
partType = validatestring(partType, validPartTypes);

%% Preparation


%% Do the job
switch partType
case {'directory', 'base', 'extension'}
    parts = extract_simple_fileparts(paths, partType);
case 'commondirectory'
    parts = extract_common_directory(paths, varargin{:});
case {'commonprefix', 'commonsuffix'}
    % First, extract file bases
    fileBases = extract_simple_fileparts(paths, 'base');

    % Next, extract file prefixes or suffices
    switch partType
        case 'commonprefix'
            parts = extract_common_prefix(fileBases, 'Delimiter', delimiter);
        case 'commonsuffix'
            parts = extract_common_suffix(fileBases, 'Delimiter', delimiter);
    end
case 'distinct'
    parts = extract_distinct_fileparts(paths);
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

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
