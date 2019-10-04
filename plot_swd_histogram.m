function [his, lines, fig] = plot_swd_histogram (varargin)
%% Plots SWD start times in a histogram
% Usage: [his, lines, fig] = plot_swd_histogram (varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       his         - histogram object created
%                   specified as a histogram object handle
%       lines       - vertical lines for infusion start times
%                   specified as a primitive line object handle array
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
%                   default == set in all_swd_sheets.m
%                   - 'SheetType': sheet type;
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'csv'
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == {'png', 'epsc2'}
%                   - Any other parameter-value pair for the histogram() function
%
% Requires:
%       cd/apply_iteratively.m
%       cd/create_default_grouping.m
%       cd/create_grouping_by_vectors.m
%       cd/extract_common_directory.m
%       cd/extract_common_suffix.m
%       cd/extract_distinct_fileparts.m
%       cd/force_matrix.m
%       cd/load_swd_sheets.m
%       cd/plot_grouped_histogram.m
%       cd/plot_vertical_line.m
%       cd/set_figure_properties.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2018-12-27 Created by Adam Lu
% 2019-09-10 Updated suffix default
% 2019-10-03 Added 'FigTypes' as an optional argument
% 

%% Hard-coded parameters
swdStr = '_SWDs';

% TODO: Make these optional arguments
recordingStartHrs = 16;         % time that recording started each day (hours)
infusionStartHrs = 20;          % time that infusion started each day (hours)

%% Default values for optional arguments
verboseDefault = true;
absoluteTimeDefault = true;     % use absolute time by default
individuallyDefault = false;    % plot all tables together by default
swdTablesDefault = '';          % set later
swdSheetPathsDefault = '';      % set later
swdFolderDefault = '';          % set later
outFolderDefault = '';          % set later
suffixDefault = '';             % set in all_swd_sheets.m
sheetTypeDefault = 'csv';       % default spreadsheet type
figTypesDefault = {'png', 'epsc2'};

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
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

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
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the histogram() function
otherArguments = iP.Unmatched;

%% Preparation
% Set default string to recognize SWD spreadsheets
if isempty(suffix)
    suffix = swdStr;
end

% Find all combined SWD spreadsheet files in the directory
if isempty(swdTables)
    % Load all tables from either paths or folder
    [swdTables, swdSheetPaths] = ...
        load_swd_sheets('Verbose', verbose, 'FilePaths', swdSheetPaths, ...
                        'Directory', swdFolder, 'SheetType', sheetType, ...
                        'Suffix', suffix);

    % Exit if nothing is loaded
    if isempty(swdTables)
        his = gobjects(1);
        lines = gobjects(1);
        fig = gobjects(1);
        return
    end
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
    [his, lines, fig] = ...
        cellfun(@(x, y) plot_swd_histogram_helper(x, y, ...
                                outFolder, absoluteTime, ...
                                recordingStartHrs, infusionStartHrs, ...
                                figTypes, otherArguments), ...
                swdTables, swdSheetPaths);
else
    [his, lines, fig] = ...
        plot_swd_histogram_helper(swdTables, swdSheetPaths, ...
                                outFolder, absoluteTime, ...
                                recordingStartHrs, infusionStartHrs, ...
                                figTypes, otherArguments);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [his, lines, fig] = ...
                plot_swd_histogram_helper(swdTables, swdSheetPaths, ...
                                        outFolder, absoluteTime, ...
                                        recordingStartHrs, infusionStartHrs, ...
                                        figTypes, otherArguments)
%% Plots a histogram from an SWD table

%% Preparation
% Extract the start times in hours (or absolute datetime)
startTimes = extract_start_times(swdTables, absoluteTime);

% Create grouping for the start times
if isdatetime(startTimes)
    % Shift all date times back by recordingStartHrs hours
    startTimesForGrouping = startTimes - hours(recordingStartHrs);

    % Find the date that the recording was started
    dates = dateshift(startTimesForGrouping, 'start', 'day');

    % Use the date as the group value
    grouping = arrayfun(@convert_to_datestr, dates, 'UniformOutput', false);

    % Simplify the grouping and generate grouping labels
    [grouping, ~, groupingLabels] = create_default_grouping('Grouping', grouping); 
else
    % Use the column number
    grouping = create_grouping_by_vectors(startTimes);

    % Extract distinct file parts
    distinctParts = extract_distinct_fileparts(swdSheetPaths);

    % Set grouping labels
    groupingLabels = distinctParts;
end

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
    binEdges = transpose(minHour:hours(1):maxHour);
else
    binEdges = transpose(minHour:1:maxHour);
end

% Choose the histogram style
if isdatetime(minHour)
    style = 'stacked';
else
    style = 'side-by-side';
end

% Set x axis limits
xLimits = [minHour, maxHour];

% Set x axis label
if isdatetime(minHour)
    xLabel = 'Date-time binned by hour';
else
    xLabel = 'Hours since start of recording';
end

% Set y axis label
yLabel = 'SWD Count';

% Set figure base name
if iscell(swdSheetPaths) && numel(swdSheetPaths) > 1
    % Extract common suffix
    commonSuffix = extract_common_suffix(swdSheetPaths);

    % Extract just the figure base
    figBase = extract_fileparts(commonSuffix, 'base');

    % Add '_all'
    figBase = [figBase, '_all'];
else
    % Extract file base
    figBase = extract_fileparts(swdSheetPaths, 'base');
end

% Set figure title
figTitle = sprintf('SWD start times for %s', strrep(figBase, '_', '\_'));

% Set figure name
figName = fullfile(outFolder, [figBase, '_histogram']);

% Compute the infusion start times
if isdatetime(minHour)
    infusionStartTimes = unique(dates) + hours(infusionStartHrs);
else
    % Start time is relative to recording start in hours
    infusionStartTimes = infusionStartHrs - recordingStartHrs;
end

%% Plot the histogram
% Create and clear figure
fig = set_figure_properties('AlwaysNew', true, 'WindowState', 'maximized');

% Plot the histogram
[his, fig] = plot_grouped_histogram(startTimes, grouping, 'Edges', binEdges, ...
                            'Style', style, 'XLimits', xLimits, ...
                            'XLabel', xLabel, 'YLabel', yLabel, ...
                            'GroupingLabels', groupingLabels, ...
                            'FigTitle', figTitle, 'FigHandle', fig, ...
                            otherArguments);

% Draw vertical lines for infusion start
lines = plot_vertical_line(infusionStartTimes, 'LineStyle', '--', ...
                            'Color', 'r', 'LineWidth', 1);

% Save the new figure
save_all_figtypes(fig, figName, figTypes);

% Close figure
% close(fig);

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
% TODO: Make this accept a partial match to startTimeStr
if iscell(swdTables)
    startTimes = cellfun(@(x) x.(startTimeStr), swdTables, ...
                            'UniformOutput', false);
else
    startTimes = swdTables.(startTimeStr);
end

% Convert cell arrays to arrays padded with NaNs
if iscell(startTimes)
    % Extract vectors padded on the right and put together as a matrix
    startTimes = force_matrix(startTimes, 'AlignMethod', 'leftAdjustPad');
end

% Convert to hours if not a datetime array
% TODO: convert_units.m
if isnum(startTimes) && ~isdatetime(startTimes) || ...
        iscellnumeric(startTimes) && ~any(cellfun(@isdatetime, startTimes))
    if iscell(startTimes)
        startTimes = cellfun(@(x) x / SEC_PER_HR, startTimes, ...
                            'UniformOutput', false);
    else
        startTimes = startTimes / SEC_PER_HR;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dateStr = convert_to_datestr(dateTime)

if isnat(dateTime)
    dateStr = '';
else
    dateStr = datestr(dateTime, 'yyyymmdd');
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

% his = histogram(startTimes, binEdges, otherArguments);
his = histogram(startTimes, binEdges);

xlim([minHour, maxHour])

startTimesForGrouping = ...
    arrayfun(@(x) addtodate(x, -recordingStartHrs, 'hour'), startTimes);

% Extract vectors padded on the right
startTimes = extract_subvectors(startTimes, 'AlignMethod', 'leftAdjustPad');

% Put together as an array
startTimes = horzcat(startTimes{:});

combinedSwdStr = '_SWDs_combined';

fig = figure('WindowState', 'maximized');
clf;

groupingLabels = replace(distinctParts, '_', '\_');

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%