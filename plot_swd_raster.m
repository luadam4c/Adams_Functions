function plot_swd_raster (varargin)
%% Plot SWDs for each channel in the current directory
% Usage: plot_swd_raster (varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Side Effects:
%       TODO
% Arguments:
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
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
%                   - 'ManualFolder': directory to look for manual output files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'SayliFolder': directory to look for Sayli output files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'AssystFolder': directory to look for Assyst output files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'SheetType': sheet type;
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'csv'
%                   - 'BarWidth': bar width relative to y value increments (0~1)
%                   must be a positive scalar between 0 and 1
%                   default == 0.6
%                   - 'LineStyle': line style of bars
%                   must be an unambiguous, case-insensitive match to one of: 
%                       '-'     - solid line
%                       '--'    - dashed line
%                       ':'     - dotted line
%                       '-.'    - dash-dotted line
%                       'none'  - no line
%                   default == '-'
%                   - 'LineWidth': line width of bars
%                   must be a positive scalar
%                   default == 2
%
% Requires:
%       cd/create_labels_from_numbers.m
%       cd/extract_common_parent.m
%       cd/find_ind_str_in_cell.m
%       cd/islinestyle.m
%       cd/issheettype.m
%       cd/plot_raster.m
%       cd/print_or_show_message.m
%
% Used by:
%       cd/plot_EEG.m
%

% File History:
% 2018-05-16 Created by Adam Lu
% 2018-05-17 Now makes ../atffiles the default atfFolder
% 2018-05-17 Changed atfFolder -> manualFolder
% 2018-11-22 Modified from plot_EEG_event_raster.m
% 2018-11-22 Changed default sheettype to .csv
% 2018-11-22 Changed startTimeStr to 'startTime'
% 2018-11-22 Changed default ManualFolder to pwd
% 2018-11-27 Now looks for SWDs.csv files recursively in SwdFolder
% 2018-12-17 Now uses create_labels_from_numbers.m
% TODO: If there is only one group, make each trace a different color
% TODO: Reanalyze data from ManualFolder, SayliFolder, AssystFolder
%           if provided
% TODO: Make the function plot_event_raster by replacing 'SWD' with 'event'
% 

%% Hard-coded constants

%% Hard-coded parameters that must be consistent with abf2mat.m
% TODO: Make these default values for optional arguments
swdStr = '_SWDs';               % string in file names for SWD spreadsheets
paramsStr = '_params';          % string in file names for params spreadsheets
startTimeStr1 = 'startTime';
startTimeStr2 = 'startTimes';
manualStr = '_Manual';          % string in file names for manual files
assystStr = '_Assyst';          % string in file names for Assyst files
sayliStr = '_Sayli';            % string in file names for Sayli files
animalStr = '_rat'; %'_animal'  % string in file names for animals
channelStr = '_channel';        % string in file names that separate channels
pieceStr = '_piece';            % string in file names that separate pieces
sweepStr = '_sweep';            % string in file names that separate sweeps
possibleSuffices = {manualStr, assystStr, sayliStr, ...
                    animalStr, channelStr, pieceStr, sweepStr};

%% Default values for optional arguments
verboseDefault = true;
swdTablesDefault = '';          % set later
swdSheetPathsDefault = '';      % set later
swdFolderDefault = '';          % set later
outFolderDefault = '';          % set later
manualFolderDefault = '';       % set later
sayliFolderDefault = '';        % set later
assystFolderDefault = '';       % set later
sheetTypeDefault = 'csv';       % default spreadsheet type
barWidthDefault = 0.6;          % default bar width rel to y value increments
lineStyleDefault = '-';         % default line style of bars
lineWidthDefault = 2;           % default line width of bars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
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
addParameter(iP, 'ManualFolder', manualFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SayliFolder', sayliFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'AssystFolder', assystFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));
addParameter(iP, 'BarWidth', barWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));
addParameter(iP, 'LineStyle', lineStyleDefault, ...
    @(x) all(islinestyle(x, 'ValidateMode', true)));
addParameter(iP, 'LineWidth', lineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));

% Read from the Input Parser
parse(iP, varargin{:});
verbose = iP.Results.Verbose;
swdTables = iP.Results.SwdTables;
swdSheetPaths = iP.Results.SwdSheetPaths;
swdFolder = iP.Results.SwdFolder;
outFolder = iP.Results.OutFolder;
manualFolder = iP.Results.ManualFolder;
sayliFolder = iP.Results.SayliFolder;
assystFolder = iP.Results.AssystFolder;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);
barWidth = iP.Results.BarWidth;
[~, lineStyle] = islinestyle(iP.Results.LineStyle, 'ValidateMode', true);
lineWidth = iP.Results.LineWidth;

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

%% Preparation
% Read in SWD tables if not provided
if isempty(swdTables)
    if ~isempty(swdSheetPaths)
        % Use full paths if provided
        [swdTables, swdSheetPaths] = ...
            load_swd_sheets('Verbose', verbose, 'FilePaths', swdSheetPaths, ...
                            'SheetType', sheetType);
    else
        % Otherwise look in the swdFolder
        [swdTables, swdSheetPaths] = ...
            load_swd_sheets('Verbose', verbose, 'Directory', swdFolder, ...
                            'SheetType', sheetType);
    end
else
    % Display warning message?
    % TODO
end

% Decide on the output folder based on swdSheetPaths
if isempty(outFolder)
    if ~isempty(swdSheetPaths)
        if iscell(swdSheetPaths)
            outFolder = extract_common_parent(swdSheetPaths);
        else
            outFolder = fileparts(swdSheetPaths);
        end
    else
        outFolder = pwd;
    end
end
if verbose
    fprintf('Outfolder is %s ...\n', outFolder);
end

% Extract the original data file bases and algorithm labels for each SWD table
if ~isempty(swdSheetPaths)
    % Get all file bases
    [~, swdSheetBases, ~] = cellfun(@fileparts, swdSheetPaths, ...
                                    'UniformOutput', false);

    % Remove everything beyond possible suffices
    %   This will be the original data file bases
    dataFileBases = swdSheetBases;
    for iSuffix = 1:numel(possibleSuffices)
        % Get the current suffix
        thisSuffix = possibleSuffices{iSuffix};

        % Split the file bases by the current suffix
        tempCells = cellfun(@(x) strsplit(x, thisSuffix, ...
                            'DelimiterType', 'RegularExpression'), ...
                        dataFileBases, 'UniformOutput', false);

        % Read all strings up to before the suffix and make them new file bases
        dataFileBases = cellfun(@(x) x{1}, tempCells, 'UniformOutput', false);
    end
else
    % Count the number of tables
    nTables = numel(swdTables);

    % Construct swdSheetBases
    swdSheetBases = create_labels_from_numbers(1:nTables, ...
                            'Prefix', 'unnamed_sweep', 'Suffix', '_SWDs');

    % Construct dataFileBases
    dataFileBases = repmat({'unnamed'}, nTables, 1);
end

% Get unique file bases
uniqueDataFileBases = unique(dataFileBases);

% Get the number of unique data file bases
nDataFileBases = numel(uniqueDataFileBases);

% Collect nSheets, labels, yLabels, eventTimes & timeLimits for each file base
nSheets = zeros(nDataFileBases, 1);
labels = cell(nDataFileBases, 1);
yLabels = cell(nDataFileBases, 1);
eventTimes = cell(nDataFileBases, 1);
timeLimits = zeros(nDataFileBases, 2);
parfor iBase = 1:nDataFileBases
    % Get the current file base
    dataFileBase = uniqueDataFileBases{iBase};

    % Use data file base as label
    %   Add escape character for underscores in the labels
    labels{iBase} = replace(dataFileBase, '_', '\_');

    % Look for all SWD spreadsheets for this file base
    indSheetsThisBase = find_ind_str_in_cell(dataFileBase, swdSheetBases, ...
                                            'SearchMode', 'substrings');

    % Get the number of SWD spreadsheets for this file base
    nSheetsThisBase = length(indSheetsThisBase);

    % Store in a vector
    nSheets(iBase) = nSheetsThisBase;

    % Extract the yLabels, eventTimes & timeLimits for each SWD spreadsheet
    yLabelsThisBase = cell(nSheetsThisBase, 1);
    eventTimesThisBase = cell(nSheetsThisBase, 1);
    timeLimitsThisBase = zeros(nSheetsThisBase, 2);
    for iSheet = 1:length(indSheetsThisBase)
        % Get the original index of this sheet
        idxSheetThis = indSheetsThisBase(iSheet);

        % Get the current SWD table
        swdsTable = swdTables{idxSheetThis};

        % Get the current SWD file base
        swdSheetBaseThis = swdSheetBases{idxSheetThis};

        % Get everything in between the data file base and swdStr
        %   and make it the Y Tick Label
        yLabelsThisBase{iSheet} = replace(swdSheetBaseThis, ...
                                    {[dataFileBase, '_'], swdStr, '_'}, ...
                                    {'', '', '\_'});

        % Extract the event start times and place in cell array
        if any(strcmp(startTimeStr1, fieldnames(swdsTable)))
            eventTimesThisSheet = swdsTable.(startTimeStr1);
        elseif any(strcmp(startTimeStr2, fieldnames(swdsTable)))
            eventTimesThisSheet = swdsTable.(startTimeStr2);
        else 
            message = sprintf('No %s field or %s field found for %s!', ...
                        startTimeStr1, startTimeStr2, swdSheetBaseThis);
            mTitle = 'Missing field warning';
            icon = 'warn';
            print_or_show_message(message, 'MTitle', mTitle, 'Icon', icon, ...
                                    'MessageMode', 'show', 'Verbose', true);
            continue;
        end

        eventTimesThisBase{iSheet} = eventTimesThisSheet;

        % Get the time limits from the original data file
        %   or from user input or just use minimum and maximum start time
        % TODO
        timeLimitsThisBase(iSheet, :) = [min(eventTimesThisSheet), ...
                                    max(eventTimesThisSheet)];
    end

    % Summarize the Y tick labels for this file base
    yLabels{iBase} = yLabelsThisBase;

    % Find the number of events for each SWD spreadsheet
    nEventsThisBase = cellfun(@length, eventTimesThisBase);

    % Find the maximum number of events
    maxNEvents = max(nEventsThisBase);

    % Put the event start times for this file base in the same array
    %   (those from each spreadsheet in a separate column)
    %   and pad the shorter vectors with NaN
    eventTimesArray = ones(maxNEvents, nSheetsThisBase) * NaN;
    for iSheet = 1:nSheetsThisBase
        % Get the number of events in this spreadsheet
        nEventsThisSheet = nEventsThisBase(iSheet);

        % Copy over the event times for this spreadsheet
        eventTimesArray(1:nEventsThisSheet, iSheet) = ...
            eventTimesThisBase{iSheet};
    end

    % Store in the eventTimes cell array
    eventTimes{iBase} = eventTimesArray;

    % Summarize the time limits for this file base
    if ~all(isnan(timeLimitsThisBase))
        timeLimits(iBase, :) = [min(timeLimitsThisBase(:, 1)), ...
                                max(timeLimitsThisBase(:, 2))];
    else
        timeLimits(iBase, :) = [NaN, NaN];
    end
end

% Find the x-axis limits that spans all time limits
if ~all(isnan(timeLimits))
    xLimits = [min(timeLimits(:, 1)), max(timeLimits(:, 2))];
else
    xLimits = [];
end

% Compute the total number of sheets
nSheetsAll = sum(nSheets);

% Put the Y Tick Labels in a single cell array
yTickLabels = vertcat(yLabels{:});

% Make sure the number of sheets match the number of labels
if numel(yTickLabels) ~= nSheetsAll
    message = 'The number of sheets don''t match the number of labels!';
    mTitle = 'Numbers mismatch error';
    icon = 'error';
    print_or_show_message(message, 'MTitle', mTitle, 'Icon', icon, ...
                            'MessageMode', 'wait', 'Verbose', true);
    return;
end

% Get outfolder name
[~, outFolderName] = fileparts(outFolder);

% Decide on the figure name
figName = fullfile(outFolder, [outFolderName, '_SWDs_raster.png']);

%% Plot the raster plot
% Create and clear figure
h = figure(15342);
clf;

% Plot the raster plot
[hLines, eventTimes, yEnds, yTicksTable] = ...
    plot_raster(eventTimes, 'BarWidth', barWidth, ...
                'LineStyle', lineStyle, 'LineWidth', lineWidth, ...
                'Labels', labels, 'YTickLabels', yTickLabels);

% Compute y axis limits
yLimits = [min(yTicksTable.locs) - 1, max(yTicksTable.locs) + 1];

% Get the handles to the first line objects for each group
hFirstLines = [];
labelsToShow = cell(0, 1);
for iBase = 1:nDataFileBases
    if ~isempty(hLines{iBase})
        hFirstLines = [hFirstLines, hLines{iBase}(1)];
        labelsToShow = [labelsToShow, labels{iBase}];
    end
end

% Set other plot properties and save figure
if ~isempty(xLimits)
    xlim(xLimits);
end
ylim(yLimits);
xlabel('Time (s)');
ylabel('Sweep #');
title(['Event start times for ', replace(outFolderName, '_', '\_')]);
legend(hFirstLines, labelsToShow, 'Location', 'SouthOutside');
saveas(h, figName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Get the folder name one level above
temp = strsplit(swdFolder, filesep);
folderName = temp{end-1};

tempCell = textscan(swdSheetBaseThis, [dataFileBase, '_%s', swdStr]);
yLabelsThisBase{iSheet} = tempCell{1};

swdFolder = fileparts(swdSheetPaths{1});

% Read in the table for this SWD file
swdsTable = readtable(swdsPath);

% Read in the tables from the files
swdTables = cellfun(@readtable, swdSheetPaths, 'UniformOutput', false);

% Make sure swdSheetPaths and swdFolder are defined
if isempty(swdSheetPaths)
    %% Find all files ending with '_SWDs.csv' under the SWD folder recursively
    [~, swdSheetPaths] = all_swd_sheets('Verbose', verbose, ...
                                        'Directory', swdFolder, ...
                                        'SheetType', sheetType);

    % Exit function if no spreadsheet files are found
    if isempty(swdSheetPaths)
        return;
    end
end

% Decide on the SWD folder to look in if not provided
if isempty(swdFolder)
    if ~isempty(swdSheetPaths)
        if iscell(swdSheetPaths)
            swdFolder = fileparts(swdSheetPaths{1});
            % TODO: Consider the case where the immediate folders are different
        else
            swdFolder = fileparts(swdSheetPaths);
        end
    else
        swdFolder = pwd;
    end
end

% Get the current SWD spreadsheet name
swdsPath = swdSheetPaths{idxSheetThis};

message = sprintf('No %s field or %s field found for %s!', ...
            startTimeStr1, startTimeStr2, swdsPath);

    swdSheetBases = arrayfun(@(x) ['unnamed_sweep', num2str(x), '_SWDs'], ...
                            transpose(1:nTables), 'UniformOutput', false);

%}
