function valueTable = compute_and_plot_values_online (valueFunc, varargin)
%% Computes and plots a value whenever a new .abf file is completed
% Usage: valueTable = compute_and_plot_values_online (valueFunc, varargin)
% Explanation:
%       Computes and plots values on the fly
% Example(s):
%       compute_and_plot_values_online(@compute_oscillation_duration)
% Outputs:
%       valueTable  - a table with two columns:
%                       inputFiles (inPathVarStr)
%                       valueUnits (valueStr)
%                   specified as a table
% Arguments:
%       valueFunc   - function used to compute values from an abf file
%                   must be a function handle that takes a file name
%                       as the first argument and return a numeric scalar
%       varargin    - 'ValueStr': string for computed values
%                   must be a string scalar or a character vector
%                   default == 'valueUnits'
%                   - 'InFolder': directory to read files from
%                   must be a string scalar or a character vector
%                   default == uigetdir()
%                   - 'InFileExt': input file extension
%                   must be a string scalar or a character vector
%                   default == 'abf'
%                   - 'PauseTime': pause time when waiting in seconds
%                   must be a nonnegative scalar
%                   default == 1 second
%                   - 'MinFileSizeBytes': minimum file size to be analyzed 
%                                           in bytes
%                   must be a nonnegative scalar
%                   default == 10000 bytes
%                   - 'AvgFunc': function used to compute average values
%                   must be a function handle that takes a numeric vector
%                       as the first argument and return a numeric scalar
%                   default == @(x) compute_phase_average(x, ...
%                                   'NToAverage', nToAverage, ...
%                                   'SelectionMethod', selectionMethod, ...
%                                   'MaxRange2Mean', maxRange2Mean)
%                   - 'NLast': number of values in the last of a phase
%                   must be a nonnegative integer scalar
%                   default == 10
%                   - 'NToAverage': number of values to average
%                   must be a positive integer scalar
%                   default == 5
%                   - 'SelectionMethod': the selection method
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'notNaN'        - select any non-NaN value
%                       'maxRange2Mean' - select vales so that the maximum 
%                                           range is within a percentage 
%                                           of the mean
%                   default == 'notNaN'
%                   - 'MaxRange2Mean': maximum percentage of range versus mean
%                   must be a nonnegative scalar
%                   default == 40%
%                   - 'ValueColor': color of marker edges and faces
%                   default == 'b'
%                   - 'Marker': type of markers
%                   default == 'o'
%                   - 'MarkerSize': size of markers
%                   default == 8
%                   - 'LineWidth': line width of markers
%                   default == 2
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == uses compute_axis_limits.m
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == uses compute_axis_limits.m
%                   - 'XLabel': label for the time axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector
%                   default == 'File Name'
%                   - 'YLabel': label(s) for the y axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector
%                   default == valueStr
%                   - 'XTickAngle': angle of X ticks in degrees
%                   default == 75
%                   - 'MaxXTicks': maximum number of X ticks
%                   default == 20
%                   - 'FigTitle': title for the figure
%                   must be a string scalar or a character vector
%                   default == ['Online values for ', yLabel]
%                   - 'InPathVarStr': input file path variable string
%                   must be a string scalar or a character vector
%                   default == 'inputFiles'
%                   - 'SaveSheetFlag': whether to save values in a spreadsheet
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveMatFlag': whether to save values in a .mat file
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'SheetType': sheet type; 
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'xlsx'
%                   - 'OutFolder': directory to place output files
%                   must be a string scalar or a character vector
%                   default == same as inFolder
%                   - Any other parameter-value pair for the plot() function
%
% Requires:
%       cd/all_files.m
%       cd/compute_phase_average.m
%       cd/compute_axis_limits.m
%       cd/create_error_for_nargin.m
%       cd/create_indices.m
%       cd/create_time_stamp.m
%       cd/extract_fileparts.m
%       cd/find_in_strings.m
%       cd/issheettype.m
%       cd/plot_horizontal_line.m
%       cd/struct2arglist.m
%
% Used by:
%       ~/Online_Stats/onlineOmight_interface.m

% File History:
% 2019-05-14 Adapted from onlineOmight_interface.m
% 2019-05-15 Finish input parser
% 2019-05-15 Now plots crosses for points used for averaging
% TODO: Save figure at the end saveFigFlag

%% Hard-coded parameters
validSelectionMethods = {'notNaN', 'maxRange2Mean'};

%% TODO: Make optional parameters
acquisFileExt = 'rsv';

%% Default values for optional arguments
valueStrDefault = 'valueUnits'; % default string for computed values
inFolderDefault = '';           % set later
inFileExtDefault = 'abf';       % read in .abf files by default
pauseTimeDefault = 1;           % pause for 1 second when waiting by default
minFileSizeBytesDefault = 10000;% analyze files > 10000 bytes by default
avgFuncDefault = [];            % set later
nLastDefault = 10;              % 10 values in last of phase by default
nToAverageDefault = 5;          % average over 5 values by default
selectionMethodDefault = 'maxRange2Mean';
                                % select values to average using 
                                %   maxRange2Mean by default
maxRange2MeanDefault = 40;      % range is not more than 40% of mean by default
valueColorDefault = 'b';        % plot blue markers by default
markerDefault = 'o';            % plot circles by default
markerSizeDefault = 8;          % marker size of 8 by default
xLimitsDefault = [];            % set later
yLimitsDefault = [];            % set later
xLabelDefault = 'File Name';    % default x-axis label
yLabelDefault = '';             % set later
lineWidthDefault = 2;           % marker line width of 2 by default
xTickAngleDefault = 75;         % x ticks are at 75 degrees angle by default
maxXTicksDefault = 20;          % maximum of 20 X ticks by default
figTitleDefault = '';           % set later
inPathVarStrDefault = 'inputFiles'; % default string for input files
saveSheetFlagDefault = true;    % save values in a spreadsheet by default
saveMatFlagDefault = false;     % don't save values in a .mat file by default
sheetTypeDefault = 'csv';       % default spreadsheet type is csv
outFolderDefault = '';          % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'valueFunc', ...
    @(x) validateattributes(x, {'function_handle'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ValueStr', valueStrDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'InFolder', inFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'InFileExt', inFileExtDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PauseTime', pauseTimeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'MinFileSizeBytes', minFileSizeBytesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'AvgFunc', avgFuncDefault, ...
    @(x) validateattributes(x, {'function_handle'}, {'2d'}));
addParameter(iP, 'NLast', nLastDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'integer', 'scalar'}));
addParameter(iP, 'NToAverage', nToAverageDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
addParameter(iP, 'SelectionMethod', selectionMethodDefault, ...
    @(x) any(validatestring(x, validSelectionMethods)));
addParameter(iP, 'MaxRange2Mean', maxRange2MeanDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'ValueColor', valueColorDefault);
addParameter(iP, 'Marker', markerDefault);
addParameter(iP, 'MarkerSize', markerSizeDefault);
addParameter(iP, 'LineWidth', lineWidthDefault);
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || iscell(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'XTickAngle', xTickAngleDefault);
addParameter(iP, 'MaxXTicks', maxXTicksDefault);
addParameter(iP, 'XLabel', xLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'InPathVarStr', inPathVarStrDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SaveSheetFlag', saveSheetFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveMatFlag', saveMatFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, valueFunc, varargin{:});
valueStr = iP.Results.ValueStr;
inFolder = iP.Results.InFolder;
inFileExt = iP.Results.InFileExt;
pauseTime =  iP.Results.PauseTime;
minFileSizeBytes =  iP.Results.MinFileSizeBytes;
avgFunc = iP.Results.AvgFunc;
nLast =  iP.Results.NLast;
nToAverage =  iP.Results.NToAverage;
selectionMethod = validatestring(iP.Results.SelectionMethod, ...
                                    validSelectionMethods);
maxRange2Mean = iP.Results.MaxRange2Mean;
valueColor = iP.Results.ValueColor;
marker = iP.Results.Marker;
markerSize = iP.Results.MarkerSize;
lineWidth = iP.Results.LineWidth;
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;
xTickAngle = iP.Results.XTickAngle;
maxXTicks = iP.Results.MaxXTicks;
xLabel = iP.Results.XLabel;
yLabel = iP.Results.YLabel;
figTitle = iP.Results.FigTitle;
inPathVarStr = iP.Results.InPathVarStr;
saveSheetFlag = iP.Results.SaveSheetFlag;
saveMatFlag = iP.Results.SaveMatFlag;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);
outFolder = iP.Results.OutFolder;

% Keep unmatched arguments for the plot() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Set default averaging function
if isempty(avgFunc)
    avgFunc = @(x) compute_phase_average(x, 'ReturnLastTrial', true, ...
                                    'NToAverage', nToAverage, ...
                                    'SelectionMethod', selectionMethod, ...
                                    'MaxRange2Mean', maxRange2Mean);
end

% Select a directory that contains .abf files
if isempty(inFolder)
    inFolder = uigetdir();
end

% Decide on the output directory
if isempty(outFolder)
    outFolder = inFolder;
end

% Set the default x-axis labels
if isempty(yLabel)
    yLabel = valueStr;
end

if isempty(figTitle)
    figTitle = ['Online values for ', yLabel];
end

% Create value labels and average value labels
valueLabel = yLabel;
avgValueLabel = ['Average ', valueLabel, ' over last ', ...
                num2str(nToAverage), ' sweeps'];

%% Initialize variables
% The .abf files already analyzed
abfPathsAnalyzed = {};

% The computed values
values = [];

% The .abf file path previously under acquisition
pathPrevUnderAcquis = '';

% The last file number analyzed
iFile = 0;

%% Create a figure
% Create the figure
fig = figure(1);

% Clear the figure
clf(fig);

% Set x axis limits
xlim('manual');
if ~isempty(xLimits) && ...
    ~(ischar(xLimits) && ~strcmpi(xLimits, 'suppress'))
    xlim(xLimits);
end

% Set y axis limits
if ~isempty(yLimits) && ...
    ~(ischar(yLimits) && ~strcmpi(yLimits, 'suppress'))
    ylim(yLimits); ylim('manual');
end

% Set x axis label
if ~strcmpi(xLabel, 'suppress')
    xlabel(xLabel);
end

% Set y axis label
if ~strcmpi(yLabel, 'suppress')
    ylabel(yLabel);
end

% Set x tick angle
if ~isempty(xTickAngle)
    xtickangle(xTickAngle);
end

% Set figure title
if ~strcmpi(figTitle, 'suppress')
    title(figTitle);
end

% Hold on
hold on;

% Initialize lines
lines = gobjects(3, 1);

% Initialize points used for the average
pointsUsed = gobjects(1, 1);

%% Analyze as long as the figure is open
while ishandle(fig)
    % Get all input files from this directory
    [abfFiles, abfPaths] = ...
        all_files('Directory', inFolder, 'Extension', inFileExt, ...
                    'ForceCellOutput', true, 'WarnFlag', false);

    % Get all acquisition files from this directory
    [~, acquisPaths] = ...
        all_files('Directory', inFolder, 'Extension', acquisFileExt, ...
                    'ForceCellOutput', true, 'WarnFlag', false);

    % Create the full path to the corresponding input file
    inPathsUnderAcquis = replace(acquisPaths, acquisFileExt, inFileExt);

    % Determine which files have not been analyzed and not under acquisition
    [abfPathsNotAnalyzed, iPathsNotAnalyzed] = ...
        setdiff(abfPaths, union(abfPathsAnalyzed, inPathsUnderAcquis));

    %   Retrieve the corresponding file structures
    abfFilesNotAnalyzed = abfFiles(iPathsNotAnalyzed);

    % Count the number of files not analyzed
    nFilesNotAnalyzed = numel(abfPathsNotAnalyzed);

    % If there are no files not analyzed and not under acquisition, 
    %   pause and continue
    if nFilesNotAnalyzed == 0
        % If there is a file under acquisition, 
        %   update the path previously under acquisition
        if ~isempty(inPathsUnderAcquis)
            pathPrevUnderAcquis = inPathsUnderAcquis{1};
        end

        % Pause and continue
        pause(pauseTime); continue
    end

    % Determine whether there is a file previously under acquisition
    idxToAnalyze = find_in_strings(pathPrevUnderAcquis, abfPathsNotAnalyzed, ...
                                    'ReturnNan', false, 'MaxNum', 1);

    % If the file is not found
    if isempty(idxToAnalyze)
        % If there are no paths under acquisition, 
        %   assume that it has been deleted
        %   and choose the first unanalyzed file to be 
        %   "previously under acquisition"
        if isempty(inPathsUnderAcquis)
            pathPrevUnderAcquis = abfPathsNotAnalyzed{1};
        end

        % Pause and continue
        pause(pauseTime); continue
    end

    % Otherwise, analyze the file that was previously under acquisition
    abfFileToAnalyze = abfFilesNotAnalyzed(idxToAnalyze);
    abfPathToAnalyze = abfPathsNotAnalyzed{idxToAnalyze};
    iFile = iFile + 1;

    % Extract the file size
    abfFileSize = abfFileToAnalyze.bytes;
    
    % Only apply the analysis function if a file size threshold is met
    if abfFileSize >= minFileSizeBytes
        % Apply the analysis function
        newValue = valueFunc(abfPathToAnalyze);
    else
        % TODO: Fix bug here?
        newValue = NaN;
    end

    % Add to the list of already analyzed .abf files
    abfPathsAnalyzed = [abfPathsAnalyzed; {abfPathToAnalyze}];

    % Add to the list of computed values
    values = [values; newValue];

    % Compute the starting and ending indices of last values
    idxLastStart = max(iFile - nLast + 1, 1);
    idxLastEnd = iFile;

    % Extract the last values
    lastValues = values(idxLastStart:idxLastEnd);

    % Compute a new average value of last values
    [newAvg, indSelectedTemp] = avgFunc(lastValues);

    % Get the indices of the selected values
    indSelected = (idxLastStart - 1) + indSelectedTemp;

    % Compute a new limits
    upperLimit = newAvg * (1 + maxRange2Mean / 200);
    lowerLimit = newAvg * (1 - maxRange2Mean / 200);

    % Print new value and average value
    fprintf("%s: %g\n", valueLabel, newValue);
    fprintf("%s: %g\n", avgValueLabel, newAvg);

    % Plot new value
    points(iFile) = ...
        plot(iFile, newValue, 'Color', valueColor, 'Marker', marker, ...
                'MarkerSize', markerSize, 'LineWidth', lineWidth, ...
                otherArguments{:});

    % Decide on x ticks
    indTicks = create_indices([1, iFile], 'MaxNum', maxXTicks);

    % Update x ticks
    xticks(indTicks);

    % Extract file bases
    fileBasesAnalyzed = extract_fileparts(abfPathsAnalyzed, 'distinct');

    % Update x ticks
    xticklabels(fileBasesAnalyzed(indTicks));

    % Expand x limits
    xlim([0, iFile + 1]);

    % Expand y limits
    if ~isempty(yLimits)
        ylim(compute_axis_limits(values, 'y'));
    end

    % Delete previous horizontal line(s) if any
    if any(ishandle(lines))
        delete(lines)
    end

    % Plot horizontal line(s)
    %   Note: Do this after x limits is updated
    lines(1, 1) = ...
        plot_horizontal_line(newAvg, 'Color', 'g', 'LineStyle', '--', ...
                                'LineWidth', 2);
    lines(2:3, 1) = ...
        plot_horizontal_line([upperLimit, lowerLimit], ...
                                'Color', 'r', 'LineStyle', '--', ...
                                'LineWidth', 2);

    % Delete previous crosses if any
    if any(ishandle(pointsUsed))
        delete(pointsUsed)
    end

    % Plot red crosses
    if ~all(isnan(indSelected))
        pointsUsed = ...
            plot(indSelected, values(indSelected), 'rx', ...
                    'MarkerSize', markerSize, 'LineWidth', lineWidth);
    end

    % Update plot
    drawnow;
end

%% Output
% Create a table
valueTable = table(abfPathsAnalyzed, values, ...
                    'VariableNames', {inPathVarStr, valueStr});

%% Save active data
if saveSheetFlag || saveMatFlag
    % Create a time stamp
    timeStamp = create_time_stamp('FormatOut', 'yyyymmddTHHMM');

    % Create the file base for output files
    valuesOutFileBase = ...
        fullfile(outFolder, [timeStamp, '_online_values_', valueStr]);
end

if saveSheetFlag
    % Create a spreadsheet path
    valuesSheetPath = [valuesOutFileBase, '.', sheetType];

    % Save values
    writetable(valueTable, valuesSheetPath);
end

if saveMatFlag
    % Create a matfile path
    valuesMatFilePath = [valuesOutFileBase, '.mat'];

    % Save values
    save(valuesMatFilePath, 'values', '-v7.3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% If there is only one file not analyzed, 
%   check if it is under acquisition 
%   (if so, a corresponding .rsv file would exist)
if nFilesNotAnalyzed == 1
    % Get the full path to the file not analyzed
    pathNotAnalyzed = abfPathToAnalyze{1};

    % Update the path under acquisition
    pathPrevUnderAcquis = pathNotAnalyzed;

    % Pause and continue
    pause(pauseTime); continue
end

% Choose a new file to be under acquisition
pathPrevUnderAcquis = abfPathsNotAnalyzed{1};

% Update the files not analyzed
abfPathsNotAnalyzed = setdiff(abfPathsNotAnalyzed, {pathPrevUnderAcquis});

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%