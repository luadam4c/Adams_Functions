function [outTable, outPath] = combine_swd_sheets (varargin)
%% Combines all files ending with '_SWDs.csv' under a directory
% Usage: [outTable, outPath] = combine_swd_sheets (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       outTable = combine_swd_sheets('Keyword', 'test3AtNight_200Hz')
%
% Outputs:
%       outTable   - combined SWD table
%                       specified as a 2-D table
%       outPath    - full path to the combined SWD spreadheet file
%                       specified as a character vector
%
% Arguments:
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Directory': directory to look for SWD table files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'Keyword': keyword for file bases
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'OutPrefix': output file prefix
%                   must be a string scalar or a character vector
%                   default == keyword or directory name
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
%       /home/Matlab/plethRO1/plethRO1_analyze.m

% File History:
% 2018-12-26 Modified from all_swd_sheets.m
% 2019-09-08 Added 'Keyword' as an optional argument
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
keywordDefault = '';
outPrefixDefault = '';
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
addParameter(iP, 'Keyword', keywordDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutPrefix', outPrefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, varargin{:});
verbose = iP.Results.Verbose;
directory = iP.Results.Directory;
keyword = iP.Results.Keyword;
outPrefix = iP.Results.OutPrefix;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);

%% Preparation
% Find all SWD spreadsheet files in the directory
[~, swdSheetPaths] = ...
    all_swd_sheets('Verbose', verbose, 'Directory', directory, ...
                    'Keyword', keyword, 'SheetType', sheetType);

% Decide on the combined SWD file prefix
if isempty(outPrefix)
    if ~isempty(keyword)
        outPrefix = keyword;
    else
        outPrefix = fileparts(directory);
    end
end

% Extract a common suffix across all files
commonSuffix = extract_fileparts(swdSheetPaths, 'commonsuffix');

% Make sure it start with '_'
commonSuffix = force_string_start(commonSuffix, '_', 'OnlyIfNonempty', true);

% Create a suffix for the combined SWD file
%   Note: this is necessary for this output file to be distinguished 
%           from input files
combinedSuffix = [commonSuffix, combinedStr];

% Create a combined SWD file name
outPath = fullfile(directory, [outPrefix, combinedSuffix, '.', sheetType]);

% Find corresponding trace file(s)
% TODO: Use extract_swd_filebases.m
% TODO: Remove common suffix and 

%% Do the job
% Concatenate the SWD sheets
outTable = vertcat_spreadsheets(swdSheetPaths);

% TODO: Add corresponding sweepStartTime to 'startTime', 'endTime' & 'duration'
%       based on trace file

%% Save the table in a file
writetable(outTable, outPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%