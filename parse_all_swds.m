function [swdTables, swdSheetPaths, ...
            swdCombinedTables, swdCombinedCsvFiles] = parse_all_swds (varargin)
%% Parses all Assyst, Sayli and manual SWD files in the current directory
% Usage: [swdTables, swdSheetPaths, ...
%           swdCombinedTables, swdCombinedCsvFiles] = parse_all_swds (varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Side Effects:
%       TODO
% Outputs:
%       TODO
%       swdTables       - SWD tables
%                       specified as a cell array of 2D tables
%       swdSheetPaths   - SWD table file names
%                       specified as a cell array of character vectors
%                       default == detected in swdFolder
% Arguments:
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'ToCombine': whether to combine spreadsheets
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Recursive': whether to search recursively
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OutFolder': the name of the directory in which 
%                                       plots will be placed
%                   must be a string scalar or a character vector
%                   default == same as the folder for the original SWD file
%                   - 'ManualFolder': directory to look for manual output files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'SayliFolder': directory to look for Sayli output files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'AssystFolder': directory to look for Assyst output files
%                   must be a string scalar or a character vector
%                   default == pwd
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
%                   - Any other parameter-value pair for parse_atf_swd() 
%                       or parse_assyst_swd()
%
% Requires:
%       cd/all_files.m
%       cd/combine_swd_sheets.m
%       cd/issheettype.m
%       cd/parse_atf_swd.m
%       cd/parse_assyst_swd.m
%
% Used by:
%       cd/combine_swd_pleth_data.m
%       cd/plot_swd_raster.m
%       /home/Matlab/plethR01/plethR01_analyze.m

% File History:
% 2018-12-26 Modified from plot_swd_raster.m
% 2019-09-08 Added 'ToCombine' as an optional parameter
% 2019-11-19 Added 'Recursive' as an optional parameter
% 2020-07-15 Now does not search recursively by default
% 2020-07-23 Made 'PieceStr' and 'Keyword' optional parameters

%% Hard-coded constants

%% Hard-coded parameters that must be consistent with abf2mat.m
assystStr = '_Assyst';          % string in file names for Assyst files

%{
sayliStr = '_Sayli';            % string in file names for Sayli files
animalStr = '_rat'; %'_animal'  % string in file names for animals
channelStr = '_channel';        % string in file names that separate channels
sweepStr = '_sweep';            % string in file names that separate sweeps
%}

%% Default values for optional arguments
verboseDefault = true;
toCombineDefault = true;
recursiveDefault = false;       % search recursively by default
outFolderDefault = '';          % set later
manualFolderDefault = '';       % set later
sayliFolderDefault = '';        % set later
assystFolderDefault = '';       % set later
pieceStrDefault = '_piece';     % string in file names that separate pieces
keywordDefault = '';
fileStartTimeDefault = [];
sheetTypeDefault = 'csv';       % default spreadsheet type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ToCombine', toCombineDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Recursive', recursiveDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ManualFolder', manualFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SayliFolder', sayliFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'AssystFolder', assystFolderDefault, ...
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
toCombine = iP.Results.ToCombine;
recursive = iP.Results.Recursive;
outFolder = iP.Results.OutFolder;
manualFolder = iP.Results.ManualFolder;
sayliFolder = iP.Results.SayliFolder;
assystFolder = iP.Results.AssystFolder;
pieceStr = iP.Results.PieceStr;
keyword = iP.Results.Keyword;
fileStartTime = iP.Results.FileStartTime;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);

% Keep unmatched arguments for parse_atf_swd() or parse_assyst_swd()
otherArguments = iP.Unmatched;

%% Preparation
% Set dependent argument defaults
if isempty(manualFolder)
    manualFolder = pwd;
end
if isempty(sayliFolder)
    sayliFolder = pwd;
end
if isempty(assystFolder)
    assystFolder = pwd;
end

% Find all .atf files in the manualFolder
[manualFiles, manualPaths] = ...
    all_files('Verbose', verbose, 'Recursive', recursive, 'WarnFlag', false, ...
                'ForceCellOutput', true, 'Directory', manualFolder, ...
                'Keyword', keyword, 'Extension', '.atf');

% Find all Assyst output files in the assystFolder
[assystFiles, assystPaths] = ...
    all_files('Verbose', verbose, 'Recursive', recursive, 'WarnFlag', false, ...
                'ForceCellOutput', true, 'Directory', assystFolder, ...
                'Keyword', keyword, 'Suffix', assystStr, 'Extension', '.txt');

% Count the number of files
nManualPaths = numel(manualPaths);
nAssystPaths = numel(assystPaths);

%% Parse all scored ATF file
% Apply parse_atf_swd.m to each .atf file
swdManualTables = cell(nManualPaths, 1);
swdManualCsvFiles = cell(nManualPaths, 1);
for iFile = 1:nManualPaths
    [swdManualTables{iFile}, swdManualCsvFiles{iFile}] = ...
        parse_atf_swd(manualPaths{iFile}, 'OutFolder', outFolder, ...
                        'SheetType', sheetType, otherArguments);
end

%% Parse all scored Assyst output file
% Apply parse_assyst_swd.m to each Assyst output file
swdAssystTables = cell(nAssystPaths, 1);
swdAssystCsvFiles = cell(nAssystPaths, 1);
parfor iFile = 1:nAssystPaths
    [swdAssystTables{iFile}, swdAssystCsvFiles{iFile}] = ...
        parse_assyst_swd(assystPaths{iFile}, 'OutFolder', outFolder, ...
                        'SheetType', sheetType, otherArguments);
end

%% Combine SWD sheets
if toCombine
    % Get all unique SWD folders
    allManualFolders = transpose(unique({manualFiles.folder}));
    allAssystFolders = transpose(unique({assystFiles.folder}));
    allSwdFolders = vertcat(allManualFolders, allAssystFolders);

    % Count the number of unique SWD folders
    nSwdFolders = numel(allSwdFolders);

    % Vertically concatenate all SWD files in the same SWD folders
    swdCombinedTables = cell(nSwdFolders, 1);
    swdCombinedCsvFiles = cell(nSwdFolders, 1);
    for iFolder = 1:nSwdFolders
    % parfor iFolder = 1:nSwdFolders
        [swdCombinedTables{iFolder}, swdCombinedCsvFiles{iFolder}] = ...
            combine_swd_sheets('Directory', allSwdFolders{iFolder}, ...
                                'Verbose', verbose, 'SheetType', sheetType, ...
                                'Keyword', keyword, 'PieceStr', pieceStr, ...
                                'FileStartTime', fileStartTime);
    end
else
    swdCombinedTables = table.empty;
    swdCombinedCsvFiles = {};
end

%% Output resultsswdSheetPaths
% Put all SWD tables together 
swdTables = vertcat(swdManualTables, swdAssystTables);
swdSheetPaths = vertcat(swdManualCsvFiles, swdAssystCsvFiles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
