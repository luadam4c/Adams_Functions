function sliceBases = all_slice_bases (varargin)
%% Retrieves all unique slice bases from the .abf files in the directory
% Usage: sliceBases = all_slice_bases (varargin)
% Explanation:
%       TODO
% Example(s):
%       sliceBases = all_slice_bases;
%       sliceBases = all_slice_bases();
% Outputs:
%       sliceBases  - unique slice bases
%                   specified as a cell array of character vectors
% Arguments:
%       varargin    - 'Directory': the directory to search in
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'Extension': file extension to limit to
%                   must be a string scalar or a character vector
%                   default == abf
%                   - 'ForceCellOutput': whether to force output as a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SortBy': how to sort the files
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'name'  - by file name
%                       'date'  - by modification date
%                       'bytes' - by file size in bytes
%                   default == 'date'
%                   - 'RegExp': regular expression to extract
%                   must be a string scalar or a character vector
%                   default == '.*slice[0-9]*'
%                   - Any other parameter-value pair for the all_files() function
%
% Requires:
%       cd/struct2arglist.m
%
% Used by:
%       cd/combine_data_from_same_slice.m

% File History:
% 2019-07-24 Moved from combine_data_from_same_slice.m
% 

%% Hard-coded parameters
validSortBys = {'name', 'date', 'bytes'};

%% Default values for optional arguments
directoryDefault = pwd;
extensionDefault = 'abf';           % look in .abf file names by default
forceCellOutputDefault = true;      % force output as a cell array by default
sortByDefault = 'date';             % sort by date by default
regExpDefault = '.*slice[0-9]*';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Extension', extensionDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ForceCellOutput', forceCellOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SortBy', sortByDefault, ...
    @(x) any(validatestring(x, validSortBys)));
addParameter(iP, 'RegExp', regExpDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, varargin{:});
directory = iP.Results.Directory;
extension = iP.Results.Extension;
forceCellOutput = iP.Results.ForceCellOutput;
sortBy = validatestring(iP.Results.SortBy, validSortBys);
regExp = iP.Results.RegExp;

% Keep unmatched arguments for the all_files() function
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
% Get all the abf file names
[~, allAbfPaths] = ...
    all_files('Directory', directory, 'Extension', extension, ...
              'SortBy', sortBy, 'ForceCellOutput', forceCellOutput, ...
              otherArguments{:});

% Extract all slice names
allSliceNames = extract_fileparts(allAbfPaths, 'base', 'RegExp', regExp);

% Get unique slice bases
sliceBases = unique(allSliceNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%