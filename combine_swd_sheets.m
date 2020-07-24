function [outTables, outPaths] = combine_swd_sheets (varargin)
%% Combines all files ending with '_SWDs.csv' and with '_piece' in the name under a directory
% Usage: [outTables, outPaths] = combine_swd_sheets (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       outTables = combine_swd_sheets;
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
%                   - 'PieceStr': string in file names that separate pieces
%                   must be a string scalar or a character vector
%                   default == '_piece'
%                   - 'Keyword': keyword for file prefixes
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'FileStartTime': start time for the original file
%                   must be a numeric scalar
%                   default == [] (nothing to add)
%                   - 'SheetType': sheet type;
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'csv'
%
% Requires:
%       cd/all_files.m
%       cd/all_swd_sheets.m
%       cd/construct_fullpath.m
%       cd/extract_fileparts.m
%       cd/force_column_cell.m
%       cd/force_string_start.m
%       cd/force_string_end.m
%       cd/is_var_in_table.m
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
% 2020-07-23 Made 'PieceStr' an optional argument 
% 2020-07-23 Made 'FileStartTime' an optional argument 
% TODO: Use 'startTimeOrig', 'endTimeOrig' & 'durationOrig'
% TODO: Detect 'FileStartTime' from corresponding abfPath
% TODO: Make 'OutFolder' an optional argument 
% 

%% Hard-coded parameters
% The following must be consistent with write_data_atf.m

%% Default values for optional arguments
verboseDefault = true;
directoryDefault = '';          % set later
outPrefixesDefault = '';
pieceStrDefault = '_piece';     % string in file names that separate pieces
keywordDefault = '';
fileStartTimeDefault = [];
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
addParameter(iP, 'PieceStr', pieceStrDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Keyword', keywordDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileStartTime', fileStartTimeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, varargin{:});
verbose = iP.Results.Verbose;
directory = iP.Results.Directory;
outPrefixes = iP.Results.OutPrefixes;
pieceStr = iP.Results.PieceStr;
keyword = iP.Results.Keyword;
fileStartTime = iP.Results.FileStartTime;
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


    % Try again if empty
    if isempty(swdSheetPaths)
        % Construct regular expression
        regExp = sprintf('^.*%s.*%s.*%s$', keyword, pieceStr, sheetType);

        % Find all SWD spreadsheet files with the piece string in the directory
        [~, swdSheetPaths] = ...
            all_files('Verbose', verbose, 'Directory', directory, ...
                        'Regexp', regExp);
    end

    % Return if nothing is found
    if isempty(swdSheetPaths)
        fprintf('Warning: No SWD sheets to combine with keyword %s in %s!\n', ...
                    keywordWithPiece, directory);
        outTables = {};
        outPaths = {};
        return
    end

    % Extract everything before the piece string
    if ~isempty(pieceStr)
        outPrefixes = extractBefore(swdSheetPaths, pieceStr);
    else
        outPrefixes = keyword;
    end
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
                                            pieceStr, fileStartTime, verbose), ...
                    uniquePrefixes, 'UniformOutput', false);

%% Outputs
if iscell(outTables) && numel(outTables) == 1
    outTables = outTables{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outTable, outPath] = ...
                combine_swd_sheets_helper (outPrefix, directory, sheetType, ...
                                            pieceStr, fileStartTime, verbose)

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
%   TODO: Sorting by name is probably not ok for more than 9 pieces!
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
if iscell(swdSheetPathsThis) && numel(swdSheetPathsThis) > 1
    outSuffix = extract_fileparts(swdSheetPathsThis, 'commonsuffix');
else
    afterPrefix = extractAfter(swdSheetPathsThis, keywordThis);
    outSuffix = extract_fileparts(afterPrefix, 'base');
end

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

% Add corresponding fileStartTime to 'startTime', 'endTime' & 'traceStartTime'
if ~isempty(fileStartTime)
    outTable = add_to_var_in_table(outTable, 'startTime', fileStartTime);
    outTable = add_to_var_in_table(outTable, 'endTime', fileStartTime);
    outTable = add_to_var_in_table(outTable, 'traceStartTime', fileStartTime);
end

%% Save the table in a file
writetable(outTable, outPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function myTable = add_to_var_in_table(myTable, varName, amount)

% Do only if varName is in the table
if is_var_in_table(varName, myTable)
    % Extract old variable value
    oldVar = myTable.(varName);

    % Add amount to old variable
    myTable.(varName) = oldVar + amount;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
