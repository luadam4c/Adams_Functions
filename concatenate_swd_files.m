function [outFilePaths] = concatenate_swd_files (varargin)
%% Parses and concatenates SWD files from different source data atf files
% Usage: [outFilePaths] = concatenate_swd_files (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       outFilePaths    - concatenated SWD file paths
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

%% Parse all scored .atf files
cellfun(@(a, b, c) parse_all_swds('Recursive', false, 'Keyword', a, ...
                                'FileStartTime', b, 'PieceStr', c), ...
        uniqueAftBases, num2cell(fileStartTimeSec), pieceStrs, ...
        'UniformOutput', false);

%% Concatenate files
% All subdirectories
subDirs = {directory};
% [~, subDirs] = all_subdirs('Prefix', 'cage');

% Find matching SWD sheets from each subdirectory
[~, swdSheets] = ...
    cellfun(@(x) find_matching_files(uniqueAftBases, 'Directory', x, ...
                        'Extension', 'csv', 'Suffix', swdSuffix), ...
        subDirs, 'UniformOutput', false);

% Define final file prefixes
finalPrefix = extract_fileparts(directory, 'base');

% Extract subdirectory names
% subDirNames = extract_fileparts(subDirs, 'dirbase');

% Construct final file names
% finalNames = strcat(finalPrefix, '_', subDirNames, '_', swdSuffix, '.csv');
finalNames = strcat(finalPrefix, '_', swdSuffix, '.csv');

% Construct final file paths
outFilePaths = fullfile(subDirs, finalNames);

% Concatenate SWD sheets
cellfun(@(x, y) vertcat_spreadsheets(x, 'OutputFileName', y), ...
        swdSheets, outFilePaths, 'UniformOutput', false);

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

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%