function [swdSheetPaths, respSheetPaths] = combine_swd_resp_data (varargin)
%% Parses and concatenates SWD and resp files from different source data .atf and .mat files
% Usage: [swdSheetPaths, respSheetPaths] = combine_swd_resp_data (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       swdSheetPaths    - concatenated SWD file paths
%                       specified as a cell array of character vectors
%
% Arguments:
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/all_files.m
%       cd/extract_fileparts.m
%       cd/find_matching_files.m
%       cd/modify_table.m
%       cd/parse_all_swds.m
%       cd/read_data_atf.m
%       cd/unique_custom.m
%       cd/vertcat_spreadsheets.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2020-08-16 Created by Adam Lu
% TODO: input parser

%% Hard coded parameters
swdSuffix = 'Manual_SWDs';
respSuffix = 'resp_table';
timeInstantsStr = 'timeInstants';
directory = pwd;
pieceStr = '_piece';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Decide on files to concatenate
% Find all .atf files in the manualFolder
[~, atfPaths] = ...
    all_files('Extension', '.atf', 'WarnFlag', false, ...
                'ForceCellOutput', true, 'SortBy', 'name');

% Find all _scored.atf files in the manualFolder
[~, scoredAtfPaths] = ...
    all_files('Extension', '.atf', 'Suffix', '_scored', ...
                'WarnFlag', false, 'ForceCellOutput', true, 'SortBy', 'name');

% Find all data .atf files
dataAftPaths = setdiff(atfPaths, scoredAtfPaths);

% Extract the file base
dataAftPathBases = extract_fileparts(dataAftPaths, 'base');

% Remove '_piece' and find the unique labels
uniqueAftBases = ...
    unique_custom(extractBefore_if_exists(dataAftPathBases, pieceStr));

%% Find the file start times
% Find the last file for each unique label
[~, lastAftPaths] = ...
    cellfun(@(x) find_last_match(x, dataAftPaths), ...
            uniqueAftBases, 'UniformOutput', false);

% Decide on file start times
if numel(uniqueAftBases) > 1
    % Extract the atfParams
    [~, atfParams] = read_data_atf('FilePaths', lastAftPaths);

    % Extract the file end times
    atfTable = struct2table(atfParams);
    timeEndSec = atfTable.timeEndSec;

    % Construct the file start times
    fileStartTimeSec = [0; timeEndSec(1:(end-1))]; 
else
    fileStartTimeSec = 0; 
end

% Construct piece strings
pieceStrs = cellfun(@(x) find_pattern(x, pieceStr), ...
                    lastAftPaths, 'UniformOutput', false);

%% Directories and files
% All subdirectories
subDirs = {directory};
% [~, subDirs] = all_subdirs('Prefix', 'cage');

% Define final file prefixes
finalPrefix = extract_fileparts(directory, 'base');

%% Concatenate resp files
% Find matching resp sheets from each subdirectory
[~, respSheets] = ...
    cellfun(@(x) find_matching_files(uniqueAftBases, 'Directory', x, ...
                        'Extension', 'csv', 'Suffix', respSuffix, ...
                        'ReturnEmpty', true), ...
            subDirs, 'UniformOutput', false);

% Concatenate respiratory measure spreadsheets
if ~isemptycell(respSheets)
    % Extract subdirectory names
    % subDirNames = extract_fileparts(subDirs, 'dirbase');

    % Construct final file names
    % respSheetNames = strcat(finalPrefix, '_', subDirNames, '_', respSuffix, '.csv');
    respSheetNames = strcat(finalPrefix, '_', respSuffix, '.csv');

    % Construct final file paths
    respSheetPaths = fullfile(subDirs, respSheetNames);

    % Concatenate respiratory measure sheets
    cellfun(@(x, y) add_time_and_vertcat(x, y, ...
                        fileStartTimeSec, timeInstantsStr), ...
            respSheets, respSheetPaths, 'UniformOutput', false);
end

%% Parse all scored .atf files
cellfun(@(a, b, c) parse_all_swds('Recursive', false, 'Keyword', a, ...
                                'FileStartTime', b, 'PieceStr', c), ...
        uniqueAftBases, num2cell(fileStartTimeSec), pieceStrs, ...
        'UniformOutput', false);

%% Concatenate SWD files
% Find matching SWD sheets from each subdirectory
[~, swdSheets] = ...
    cellfun(@(x) find_matching_files(uniqueAftBases, 'Directory', x, ...
                        'Extension', 'csv', 'Suffix', swdSuffix, ...
                        'ReturnEmpty', true), ...
            subDirs, 'UniformOutput', false);

% Concatenate SWD sheets
if ~isemptycell(swdSheets)
    % Extract subdirectory names
    % subDirNames = extract_fileparts(subDirs, 'dirbase');

    % Construct final file names
    % swdSheetNames = strcat(finalPrefix, '_', subDirNames, '_', swdSuffix, '.csv');
    swdSheetNames = strcat(finalPrefix, '_', swdSuffix, '.csv');

    % Construct final file paths
    swdSheetPaths = fullfile(subDirs, swdSheetNames);

    % Concatenate SWD sheets
    cellfun(@(x, y) vertcat_spreadsheets(x, 'OutputFileName', y), ...
            swdSheets, swdSheetPaths, 'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idxMatch, match] = find_last_match(str, list)
%% TODO: Integrate into find_first_match.m

idxMatch = find(contains(list, str), 1, 'last');

match = list{idxMatch};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = extractBefore_if_exists (str, endStr)
%% TODO: Pull out as its own function

if iscell(str)
    str = cellfun(@(x) extractBefore_if_exists_helper(x, endStr), ...
                    str, 'UniformOutput', false);
else
    str = arrayfun(@(x) extractBefore_if_exists_helper(x, endStr), str);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = extractBefore_if_exists_helper (str, endStr)

if contains(str, endStr)
    str = extractBefore(str, endStr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = find_pattern (str, pattern)

if contains(str, pattern)
    output = pattern;
else
    output = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function add_time_and_vertcat (respSheets, respSheetPath, ...
                                fileStartTimeSec, timeInstantsStr)

% Extract tables
respTables = cellfun(@readtable, respSheets, 'UniformOutput', false);

% Add file start times to timeInstantsStr
respTables = cellfun(@(x, y) modify_table(x, @(a) a + y, ...
                                    'VariableNames', timeInstantsStr), ...
                    respTables, num2cell(fileStartTimeSec), ...
                    'UniformOutput', false);

% Vertically concatenate tables
vertcat_spreadsheets(respTables, 'OutputFileName', respSheetPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
