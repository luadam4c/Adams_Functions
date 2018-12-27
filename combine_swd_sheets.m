function [combinedTables, combinedPaths] = combine_swd_sheets (varargin)
%% Combines all files ending with '_SWDs.csv' under a directory
% Usage: [combinedTables, combinedPaths] = combine_swd_sheets (varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       swdSheetFiles   - combined SWD table(s)
%                       specified as a 2-D table or a cell array of 2-D tables
%       combinedPaths   - full path(s) to the combined SWD spreadheet file
%                       specified as a character vector 
%                           or a cell array of character vectors
% Arguments:
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Directory': directory to look for SWD table files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'SheetType': sheet type;
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'csv'
%                   
% Requires:
%       cd/all_swd_sheets.m
%       cd/vertcat_spreadsheets.m
%
% Used by:
%       cd/parse_all_swds.m

% File History:
% 2018-12-26 Modified from all_swd_sheets.m
% TODO: Use 'startTimeOrig', 'endTimeOrig' & 'durationOrig'
% TODO: Add sweepStartTime to 'startTime', 'endTime' & 'duration'
% TODO: Make 'SweepStartTime' an optional argument 
%       and detect from corresponding abfPath
% TODO: Make 'OutFolder' an optional argument 
% 

%% Hard-coded parameters
combinedStr = '_combined';

%% Default values for optional arguments
verboseDefault = true;
directoryDefault = '';          % set later
sheetTypeDefault = 'csv';       % default spreadsheet type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, varargin{:});
verbose = iP.Results.Verbose;
directory = iP.Results.Directory;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);

%% Preparation
% Find all SWD spreadsheet files in the directory
[~, swdSheetPaths] = ...
    all_swd_sheets('Verbose', verbose, 'Directory', directory, ...
                    'SheetType', sheetType);

% Get the directory name
dirName = fileparts(directory);

% Extract a common suffix across all files
commonSuffix = extract_fileparts(paths, 'commonsuffix');

% Create a suffix for the combined SWD file
combinedSuffix = [commonSuffix, combinedStr];

% Create a combined SWD file name
combinedPath = fullfile(directory, [dirName, combinedSuffix, '.', sheetType]);

% Find corresponding .abf file(s)
% TODO: Use extract_swd_filebases.m
% TODO: Remove common suffix and 

%% Do the job
% Concatenate the SWD sheets
combinedTable = vertcat_spreadsheets(swdSheetPaths);

% TODO: Add corresponding sweepStartTime to 'startTime', 'endTime' & 'duration'
%       based on abfPath

%% Save the table in a file
writetable(combinedTable, combinedPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%