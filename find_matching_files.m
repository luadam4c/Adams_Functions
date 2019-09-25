function varargout = find_matching_files (fileParts, varargin)
%% Finds matching files from file parts
% Usage: [files, fullPaths, distinctParts] = find_matching_files (fileParts, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [~, matPaths] = all_files('Ext', 'mat');
%       [csvFiles, csvPaths] = find_matching_files(matPaths, 'Extension', 'csv');
%       [wmvFiles, wmvPaths] = find_matching_files(matPaths, 'Extension', 'wmv');
%
% Outputs:
%       files       - file structure(s) for the files
%                   specified as a structure array with fields:
%                       name
%                       folder
%                       date
%                       bytes
%                       isdir
%                       datenum
%       fullPaths   - full path(s) to the files
%                   specified as a column cell array of character vectors
%
% Arguments:
%       fileParts   - file parts to match (can be )
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%       varargin    - 'PartType': part type to match
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'Prefix'    - solid line
%                       'Keyword'   - dashed line
%                       'Suffix'    - dotted line
%                       'Extension' - extension
%                   default == 'Keyword'
%                   - Any other parameter-value pair for all_files()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_distinct_fileparts.m
%
% Used by:
%       cd/load_matching_sheets.m

% File History:
% 2019-09-25 Created by Adam Lu
% 

%% Hard-coded parameters
validPartTypes = {'Prefix', 'Keyword', 'Suffix', 'Extension'};

%% Default values for optional arguments
partTypeDefault = 'Keyword';

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
addRequired(iP, 'fileParts', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['strs5 must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PartType', partTypeDefault, ...
    @(x) any(validatestring(x, validPartTypes)));

% Read from the Input Parser
parse(iP, fileParts, varargin{:});
partType = validatestring(iP.Results.PartType, validPartTypes);

% Keep unmatched arguments for the all_files() function
otherArguments = iP.Unmatched;

%% Extract distinct parts
% Extract distinct file parts
distinctParts = extract_distinct_fileparts(fileParts);

%% Do the job
% Find one matching file for each file part
[filesCell, fullPaths] = ...
    cellfun(@(x) all_files(partType, x, 'MaxNum', 1, ...
                        'ForceCellOutput', false, otherArguments), ...
            distinctParts, 'UniformOutput', false);

% Try to convert to an array
%   Note: this fails if a cell is empty
try
    files = cellfun(@(x) x, filesCell);
catch
    disp([mfilename, ': Some files were not found!']);
    files = filesCell;
end

% Get first output
varargout{1} = files;
varargout{2} = fullPaths;
varargout{3} = distinctParts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%