function [outTables, outPaths] = combine_swd_sheets (varargin)
%% Combines all files ending with '_SWDs.csv' under a directory
% Usage: [outTables, outPaths] = combine_swd_sheets (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       outTables = combine_swd_sheets('Keyword', 'test3AtNight_200Hz')
%
% Outputs:
%       outTables   - combined SWD table(s)
%                       specified as a 2-D table
%       outPaths    - full path(s) to the combined SWD spreadheet file
%                       specified as a character vector
%
% Arguments:
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Directory': directory to look for SWD table files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'OutPrefixes': output file prefixes
%                       Note: These are the strings before the piece string
%                   must be a string scalar or a character vector
%                   default == all those detected
%                   - 'Keyword': keyword for file prefixes
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'SheetType': sheet type;
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'csv'
%
% Requires:
%       cd/all_swd_sheets.m
%       cd/construct_fullpath.m
%       cd/extract_fileparts.m
%       cd/force_column_cell.m
%       cd/force_string_start.m
%       cd/force_string_end.m
%       cd/vertcat_spreadsheets.m
%
% Used by:
%       cd/parse_all_swds.m
%       /home/Matlab/plethR01/plethR01_analyze.m

% File History:
% 2018-12-26 Modified from all_swd_sheets.m
% 2019-09-08 Added 'Keyword' as an optional argument
% 2019-09-10 Now combines only groups of sheets with pieceStr in the name
% 2019-09-10 Now omits the combinedStr
% TODO: Use 'startTimeOrig', 'endTimeOrig' & 'durationOrig'
% TODO: Add sweepStartTime to 'startTime', 'endTime' & 'duration'
% TODO: Make 'SweepStartTime' an optional argument 
%       and detect from corresponding abfPath
% TODO: Make 'OutFolder' an optional argument 
% 

%% Hard-coded parameters
% The following must be consistent with write_data_atf.m
pieceStr = '_piece';            % string in file names that separate pieces

%% Default values for optional arguments
verboseDefault = true;
directoryDefault = '';          % set later
outPrefixesDefault = '';
keywordDefault = '';
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
addParameter(iP, 'OutPrefixes', outPrefixesDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Keyword', keywordDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, varargin{:});
verbose = iP.Results.Verbose;
directory = iP.Results.Directory;
outPrefixes = iP.Results.OutPrefixes;
keyword = iP.Results.Keyword;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);

%% Preparation
% Decide on the combined SWD file prefixes
%   Note: These are the strings before the piece string
if isempty(outPrefixes)
    % Make sure the keyword ends with the piece string
    %   Note: this assumes that the keyword is always before the piece string
    keywordWithPiece = force_string_end(keyword, pieceStr);

    % Find all SWD spreadsheet files with the piece string in the directory
    [~, swdSheetPaths] = ...
        all_swd_sheets('Verbose', verbose, 'Directory', directory, ...
                        'Keyword', keywordWithPiece, 'SheetType', sheetType);

    % Extract everything before the piece string
    outPrefixes = extractBefore(swdSheetPaths, pieceStr);
end

% Force as a cell array
outPrefixes = force_column_cell(outPrefixes);

% Get unique prefixes
uniquePrefixes = unique(outPrefixes);

% If there is none found, return
if numel(uniquePrefixes) == 0
    outTables = {};
    outPaths = {};
    return;
end

%% Do the job
% Combine the spreadsheets
[outTables, outPaths] = ...
    cellfun(@(x) combine_swd_sheets_helper(x, directory, sheetType, ...
                                            pieceStr, verbose), ...
                    uniquePrefixes, 'UniformOutput', false);

%% Outputs
if iscell(outTables) && numel(outTables) == 1
    outTables = outTables{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outTable, outPath] = ...
                combine_swd_sheets_helper (outPrefix, directory, ...
                                            sheetType, pieceStr, verbose)

% Display message
if verbose
    fprintf('Combining SWD spreadsheets with prefix %s ...\n', outPrefix);
end

% Extract just the prefix
keywordThis = extract_fileparts(outPrefix, 'dirbase');

% Make sure the keyword ends with the piece string
%   Note: this assumes that the keyword is always before the piece string
keywordWithPieceThis = force_string_end(keywordThis, pieceStr);

% Find all SWD spreadsheet files with the output prefix and piece string 
%   in the directory, sorting by name
[~, swdSheetPathsThis] = ...
    all_swd_sheets('Verbose', verbose, 'Directory', directory, ...
                    'Keyword', keywordWithPieceThis, 'SheetType', sheetType, ...
                    'SortBy', 'name');

if isempty(swdSheetPathsThis)
    outTable = table.empty;
    outPath = '';
    return
end

% Extract a common suffix across all files
outSuffix = extract_fileparts(swdSheetPathsThis, 'commonsuffix');

% Make sure it start with '_'
outSuffixWithUnderScore = ...
    force_string_start(outSuffix, '_', 'OnlyIfNonempty', true);

% Create a combined SWD file name
outName = [outPrefix, outSuffixWithUnderScore, '.', sheetType];

% Make sure it's a full path
outPath = construct_fullpath(outName, 'Directory', directory);

% Find corresponding trace file(s)
% TODO: Use extract_swd_tracefiles.m
% TODO:     Remove common suffix

%% Do the job
% Concatenate the SWD sheets
outTable = vertcat_spreadsheets(swdSheetPathsThis);

% TODO: Add corresponding sweepStartTime to 'startTime', 'endTime' & 'duration'
%       based on trace file

%% Save the table in a file
writetable(outTable, outPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

combinedStr = '_combined';
% Create a suffix for the combined SWD file
%   Note: this is necessary for this output file to be distinguished 
%           from input files
combinedSuffix = [commonSuffix, combinedStr];
% Create a combined SWD file name
outPath = fullfile(directory, [outPrefixes, combinedSuffix, '.', sheetType]);

if isempty(outPrefixes)
    if ~isempty(keyword)
        outPrefixes = keyword;
    else
        outPrefixes = fileparts(directory);
    end
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
