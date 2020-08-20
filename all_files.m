function varargout = all_files (varargin)
%% Returns all the files in a given directory (optionally recursive) that matches a prefix, keyword, suffix or extension
% Usage: [files, fullPaths] = all_files (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [files, fullPaths] = all_files;
%       [files, fullPaths] = all_files('SortBy', 'date');
%       [files, fullPaths] = all_files('Recursive', true);
%       files = all_files('Dir', {'parse_lts', 'plot_relative_events'}, 'Recursive', true)
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
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'WarnFlag': whether to warn if no files found
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Recursive': whether to search recursively
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'SubDirInstead': whether to look for 
%                                       subdirectories instead
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ForceCellOutput': whether to force output as a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Directory': the directory(ies) to search in
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == pwd
%                   - 'Prefix': prefix the file name must have
%                   must be a string scalar or a character vector
%                   default == no limits
%                   - 'Keyword': keyword the file name must contain
%                   must be a string scalar or a character vector
%                   default == no limits
%                   - 'Suffix': suffix the file name must have
%                   must be a string scalar or a character vector
%                   default == no limits
%                   - 'Extension': file extension to limit to
%                   must be a string scalar or a character vector
%                   default == no limits
%                   - 'RegExp': regular expression to limit to
%                   must be a string scalar or a character vector
%                   default == no limits
%                   - 'SortBy': how to sort the files
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'name'  - by file name
%                       'date'  - by modification date
%                       'bytes' - by file size in bytes
%                   default == 'name'
%                   - 'MaxNum': maximum number of files to find
%                   must be empty or a positive integer scalar
%                   default == [] (no restriction)
%
% Requires:
%       cd/construct_and_check_fullpath.m
%       cd/extract_fullpaths.m
%       cd/print_cellstr.m
%       cd/print_or_show_message.m
%
% Used by:
%       cd/all_subdirs.m
%       cd/all_swd_sheets.m
%       cd/atf2sheet.m
%       cd/combine_data_from_same_slice.m
%       cd/combine_swd_sheets.m
%       cd/compile_mod_files.m
%       cd/combine_swd_resp_data.m
%       cd/create_pleth_EEG_movies.m
%       cd/find_matching_files.m
%       cd/force_string_start.m
%       cd/read_matching_sheets.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_network_analyze_spikes.m
%       cd/m3ha_network_single_neuron.m
%       cd/m3ha_network_plot_essential.m
%       cd/m3ha_simulate_population.m
%       cd/m3ha_pfiles2csv.m
%       cd/m3ha_plot_figure02.m
%       cd/m3ha_plot_figure08.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/metabolismR01_plot_chevrons.m
%       cd/parse_all_abfs.m
%       cd/parse_all_multiunit.m
%       cd/parse_all_swds.m
%       cd/parse_current_family.m
%       cd/plot_all_abfs.m
%       cd/plot_calcium_imaging_traces.m
%       cd/plot_repetitive_protocols.m
%       cd/plot_traces_EEG.m
%       cd/read_data_atf.m
%       cd/update_file_base_in_matfiles.m
%       ~/plethR01/plethR01_analyze.m
%       ~/FluoroSNNAP/FluroSNNAP.m
%       /media/adamX/m3ha/optimizer4gabab/singleneuronfitting63.m
%       

% File History:
% 2018-10-04 Modified from all_subdirs.m
% 2018-11-21 Added 'Prefix', 'Keyword', 'Suffix', 'RegExp' as optional arguments
% 2018-11-26 Added 'Recursive' as an optional flag
% 2018-12-26 Added 'ForceCellOutput' as an optional argument
% 2019-03-15 Fixed the case when extension is not provided
% 2019-05-16 Added 'WarnFlag' as an optional flag
% 2019-05-21 Added 'SortBy' as an optional argument
% 2019-09-05 Add 'MaxNum' as an optional argument
% 2019-11-24 Now allows 'Directory' to be multiple directories
% 2019-12-13 Added 'SubDirInstead' as an optional argument 
% 2019-12-31 Fixed usage of the 'Prefix' option
% 2019-01-22 Added usage of add_escape_char()
% 2019-01-30 Fixed 'SubDirInstead'
% 2019-01-30 Now sorts by 'datenum' if the user wants to sort by 'date'
% 2019-01-30 Now allows a '.' to be in the prefix, keyword, suffix or extension
% TODO: Fix bug when a dot is in the folder name

%% Hard-coded parameters
validSortBys = {'name', 'date', 'bytes', 'datenum'};

%% Default values for optional arguments
verboseDefault = false;         % don't print to standard output by default
warnFlagDefault = true;         % warn if no files found by default
recursiveDefault = false;       % don't search recursively by default
subDirInsteadDefault = false;   % look for files by default
forceCellOutputDefault = false; % don't force output as a cell array by default
directoryDefault = '';          % construct_and_check_fullpath('') == pwd
prefixDefault = '';             % set later
keywordDefault = '';            % set later
suffixDefault = '';             % set later
extensionDefault = '';          % set later
regExpDefault = '';             % set later
sortByDefault = 'name';         % sort by name by default
maxNumDefault = [];             % no restriction by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'WarnFlag', warnFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Recursive', recursiveDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SubDirInstead', subDirInsteadDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ForceCellOutput', forceCellOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['Directory must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Keyword', keywordDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Suffix', suffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Extension', extensionDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'RegExp', regExpDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SortBy', sortByDefault, ...
    @(x) any(validatestring(x, validSortBys)));
addParameter(iP, 'MaxNum', maxNumDefault, ...       % maximum number of indices
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'MaxNum must be either empty or a positive integer scalar!'));

% Read from the Input Parser
parse(iP, varargin{:});
verbose = iP.Results.Verbose;
warnFlag = iP.Results.WarnFlag;
recursive = iP.Results.Recursive;
subDirInstead = iP.Results.SubDirInstead;
forceCellOutput = iP.Results.ForceCellOutput;
directory = iP.Results.Directory;
prefix = iP.Results.Prefix;
keyword = iP.Results.Keyword;
suffix = iP.Results.Suffix;
extension = iP.Results.Extension;
regExp = iP.Results.RegExp;
sortBy = validatestring(iP.Results.SortBy, validSortBys);
maxNum = iP.Results.MaxNum;

% Make sure the directory is an existing full path
[directory, dirExists] = construct_and_check_fullpath(directory);
if ~all(dirExists)
    varargout{1} = [];
    varargout{2} = {};
    return
end

%% Preparation
% Make sure the extension begins with a dot
extension = force_string_start(extension, '.', 'OnlyIfNonempty', true);

if subDirInstead
    detectType = 'subdirectories';
else
    detectType = 'files';
end

%% Find files
% Get or check the regular expression to match
if isempty(regExp)
    % Add an escape character for special characters
    [prefix, keyword, suffix, extension] = ...
        argfun(@(x) add_escape_char(x), prefix, keyword, suffix, extension);

    if ~isempty(extension)
        % Match the prefix, keyword, suffix and extension
        regExp = sprintf('^%s.*%s.*%s%s$', prefix, keyword, suffix, extension);
    else
        if subDirInstead
            % Match the prefix, keyword, suffix
            regExp = sprintf('^%s.*%s.*%s$', prefix, keyword, suffix);
        else
            % Match the prefix, keyword, suffix
            regExp = sprintf('^%s.*%s.*%s[.].*$', prefix, keyword, suffix);
        end
    end
else
    % Display warning if an extension is provided
    if ~isempty(prefix)
        fprintf('Warning: A regular expression will override the prefix!\n');
    end
    if ~isempty(keyword)
        fprintf('Warning: A regular expression will override the keyword!\n');
    end
    if ~isempty(suffix)
        fprintf('Warning: A regular expression will override the suffix!\n');
    end
    if ~isempty(extension)
        fprintf('Warning: A regular expression will override the extension!\n');
    end
end

% Get a list of all files and subdirectories in this directory 
filesOrDirs = all_files_or_dirs(directory, recursive);

% Get a logical vector that tells which entries are directories
isDir = transpose([filesOrDirs.isdir]);

% Get all file or directory names
names = transpose({filesOrDirs.name});

% Get a logical vector that tells which entries matches the regular expression
if ~isempty(regExp)
    % Test whether each matches the regular expression
    isMatch = cellfun(@any, regexpi(names, regExp));
else
    % All files will be considered matched
    isMatch = true(size(filesOrDirs));
end

% Get a logical vector that tells which entries are irrelevant ('.' or '..')
isIrrelevant = cellfun(@(x) any(strcmp(x, {'.', '..'})), names);

% Keep only those that are not directories and 
%   are matches to the regular expression
if subDirInstead
    files = filesOrDirs(isDir & ~isIrrelevant & isMatch);
else
    files = filesOrDirs(~isDir & isMatch);
end

% Sort by date or bytes if requested
if ~strcmpi(sortBy, 'name')
    % If the user wants to sort by date, actually sort by datenum
    if strcmpi(sortBy, 'date')
        sortBy = 'datenum';
    end

    % Convert the struct array to a table
    filesTable = struct2table(files);

    % Sort the table by the requested field
    filesTableSorted = sortrows(filesTable, sortBy); 

    % Change it back to a struct array
    files = table2struct(filesTableSorted);
end

% Restrict to maximum number of files
if ~isempty(maxNum)
    % Find the index of the last file tp return
    idxEnd = min(numel(files), maxNum);

    % Restrict to those files
    files = files(1:idxEnd);
end

% Get first output
varargout{1} = files;

% Extract the full paths
if nargout >= 2
    varargout{2} = extract_fullpaths(files, 'ForceCellOutput', forceCellOutput);
end

%% Print to standard output
% Count the number of files
nFiles = numel(files);

% Print appropriate message
directoryStr = print_cellstr(directory, 'OmitNewline', true, ...
                            'OmitBraces', true, 'OmitQuotes', true, ...
                            'ToPrint', false);
if nFiles == 0 && warnFlag
    message = sprintf('No %s with pattern %s found in %s!!\n', ...
                        detectType, regExp, directoryStr);
    mTitle = 'No files found warning';
    icon = 'warn';
    print_or_show_message(message, 'MTitle', mTitle, 'Icon', icon, ...
                            'MessageMode', 'show', 'Verbose', verbose, ...
                            'CreateMode', 'replace');
elseif verbose
    fprintf('%d %s with pattern %s found in %s!\n', ...
            nFiles, detectType, regExp, directoryStr);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filesOrDirs = all_files_or_dirs (directory, recursive)

% Look under subdirectories if requested
if recursive
    directoryForDir = fullfile(directory, '**');
else
    directoryForDir = directory;
end

% Use dir repeatedly
if iscell(directoryForDir)
    % Use dir repeatedly
    filesOrDirsAll = cellfun(@dir, directoryForDir, 'UniformOutput', false);
elseif isstring(directoryForDir)
    % Use dir repeatedly
    filesOrDirsAll = arrayfun(@dir, directoryForDir, 'UniformOutput', false);
else
    filesOrDirsAll = {dir(directoryForDir)};
end

% Vertically concatenate all File objects
filesOrDirs = vertcat(filesOrDirsAll{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = add_escape_char (str)
%% Adds an escape character for all special characters in the string
% TODO: Pull out as a function
% TODO for SHINSHIN: Examine usage of regexp for other special characters

specialChars = {'[', ']', '.'};
specialCharsEscaped = {'\[', '\]', '[.]'};

str = replace(str, specialChars, specialCharsEscaped);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
