function [his, fig] = plot_swd_histogram (varargin)
%% Plots SWD start times in a histogram
% Usage: [his, fig] = plot_swd_histogram (varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       his        - histogram object created
%                   specified as a histogram object handle
%       fig         - figure containing the histogram object
%                   specified as a figure object handle
% Arguments:
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'AbsoluteTime': whether to use absolute time
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Individually': whether tables are plotted individually
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'SwdTables': SWD tables
%                   must be a 2D table or a cell array of 2D tables
%                   default == read from swdSheetPaths
%                   - 'SwdSheetPaths': SWD table file names
%                   must be a string scalar or a character vector
%                       or a cell array of character vectors
%                   default == detected in swdFolder
%                   - 'SwdFolder': directory to look for SWD table files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'OutFolder': the name of the directory in which 
%                                       plots will be placed
%                   must be a string scalar or a character vector
%                   default == a subdirectory named by {fileName}_traces in pwd
%                   - 'Suffix': suffix the file name must have
%                   must be a string scalar or a character vector
%                   default == '_SWDs_combined'
%                   - 'SheetType': sheet type;
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'csv'
%                   - Any other parameter-value pair for the histogram() function
%
% Requires:
%       cd/apply_iteratively.m
%       cd/load_swd_sheets.m
%       cd/extract_common_directory.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2018-12-27 Created by Adam Lu
% TODO: 'UseAbsoluteTime' as an optional argument
% 

%% Hard-coded parameters
combinedSwdStr = '_SWDs_combined';

%% Default values for optional arguments
verboseDefault = true;
absoluteTimeDefault = true;
individuallyDefault = false;    % plot all tables together by default
swdTablesDefault = '';          % set later
swdSheetPathsDefault = '';      % set later
swdFolderDefault = '';          % set later
outFolderDefault = '';          % set later
outFolderDefault = '';          % set later
suffixDefault = '';             % set later
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
addParameter(iP, 'AbsoluteTime', absoluteTimeDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Individually', individuallyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SwdTables', swdTablesDefault, ...
    @(x) assert(isempty(x) || istable(x) || ischar(x) || iscell(x), ...
        'SwdTables must be a table or a cell array of tables!'));
addParameter(iP, 'SwdSheetPaths', swdSheetPathsDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['SwdSheetPaths must be a a character array, a string array ', ...
        'or cell array of character arrays!']));
addParameter(iP, 'SwdFolder', swdFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Suffix', suffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, varargin{:});
verbose = iP.Results.Verbose;
absoluteTime = iP.Results.AbsoluteTime;
individually = iP.Results.Individually;
swdTables = iP.Results.SwdTables;
swdSheetPaths = iP.Results.SwdSheetPaths;
swdFolder = iP.Results.SwdFolder;
outFolder = iP.Results.OutFolder;
suffix = iP.Results.Suffix;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);

% Keep unmatched arguments for the histogram() function
otherArguments = iP.Unmatched;

%% Preparation
% Set default string to recognize SWD spreadsheets
if isempty(suffix)
    suffix = combinedSwdStr;
end

% Find all combined SWD spreadsheet files in the directory
if isempty(swdTables)
    % Load all tables from either paths or folder
    [swdTables, swdSheetPaths] = ...
        load_swd_sheets('Verbose', verbose, 'FilePaths', swdSheetPaths, ...
                        'Directory', swdFolder, 'SheetType', sheetType, ...
                        'Suffix', suffix);
end

% Decide on the output folder based on swdSheetPaths
if isempty(outFolder)
    if ~isempty(swdSheetPaths)
        outFolder = extract_common_directory(swdSheetPaths);
    else
        outFolder = pwd;
    end
end
if verbose
    fprintf('Outfolder is %s ...\n', outFolder);
end

%% Plot all SWD tables
if iscell(swdTables) && individually
    [his, fig] = cellfun(@(x) plot_swd_histogram_helper(x, absoluteTime, ...
                                outFolder, swdSheetPaths, otherArguments), ...
                        swdTables);
else
    [his, fig] = plot_swd_histogram_helper(swdTables, absoluteTime, ...
                                outFolder, swdSheetPaths, otherArguments);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [his, fig] = plot_swd_histogram_helper(swdTables, absoluteTime, ...
                                    outFolder, swdSheetPaths, otherArguments)
%% Plots a histogram from an SWD table

%% Preparation
% Extract the start times in hours (or absolute datetime)
startTimes = extract_start_times(swdTables, absoluteTime);

% Find the minimum and maximum times
minTime = apply_iteratively(@min, startTimes);
maxTime = apply_iteratively(@max, startTimes);

% Find the nearest hours
if isdatetime(minTime)
    minHour = dateshift(minTime, 'start', 'hour');
    maxHour = dateshift(maxTime, 'start', 'hour');
else
    minHour = floor(minTime);
    maxHour = ceil(maxTime);
end

% Create bin edges
if isdatetime(minHour)
    binEdges = minHour:hours(1):maxHour;
else
    binEdges = minHour:1:maxHour;
end

%% Plot the histogram
% Create and clear figure
fig = figure('WindowState','maximized');
clf;

% Plot the histogram
% his = histogram(startTimes, binEdges, otherArguments);
his = histogram(startTimes, binEdges);

% Set other plot properties
xlim([minHour, maxHour])
% title('SWD start times for %s', swdSheetPaths)

% Save figure



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function startTimes = extract_start_times (swdTables, absoluteTime)
%% Extracts start times

SEC_PER_HR = 3600;

% Decide on the variable name
if absoluteTime
    % In absolute time
    startTimeStr = 'startTimeOrig';
else
    % In seconds
    startTimeStr = 'startTime';
end

% Extract the start times
% TODO: Make this a function extract_vars.m
if iscell(swdTables)
    startTimes = cellfun(@(x) x.(startTimeStr), swdTables, ...
                            'UniformOutput', false);
else
    startTimes = swdTable.(startTimeStr);
end

% Convert cell arrays to arrays filled with NaN
% TODO: Make this a function force_matrix.m
if iscell(startTimes)
    % Extract vectors padded on the right
    startTimes = extract_subvectors(startTimes, 'AlignMethod', 'leftAdjustPad');

    % Put together as an array
    startTimes = horzcat(startTimes{:});
end

% Convert to hours if not a datetime array
if ~isdatetime(startTimes)
    startTimes = startTimes / SEC_PER_HR;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

minHour = floor(minTime / SEC_PER_HR) * SEC_PER_HR;
maxHour = ceil(minTime / SEC_PER_HR) * SEC_PER_HR;

% Extract the start times
if iscell(swdTables)
    startTimes = cellfun(@(x) extract_start_times_helper(x, startTimeStr), ...
                            swdTables, 'UniformOutput', false);
else
    startTimes = extract_start_times_helper(swdTables, startTimeStr);
end
function startTimes = extract_start_times_helper (swdTable, startTimeStr)
% Extract the start times
startTimes = swdTable.(startTimeStr);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%