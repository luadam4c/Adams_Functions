function compute_and_plot_values_online (valueFunc, varargin)
%% Computes and plots a value whenever a new .abf file is completed
% Usage: compute_and_plot_values_online (valueFunc, avgFunc (opt), varargin)
% Explanation:
%       TODO
% Example(s):
%       compute_and_plot_values_online(@compute_oscillation_duration)
% Arguments:
%       valueFunc   - function used to compute values from an abf file
%                   must be a function handle that takes a file name
%                       as the first argument and return a numeric scalar
%       avgFunc     - (opt) function used to compute average values
%                   must be a function handle that takes a numeric vector
%                       as the first argument and return a numeric scalar
%                   default == @(x) compute_phase_average(x, ...
%                                   'NToAverage', nToAverage, ...
%                                   'SelectionMethod', selectionMethod, ...
%                                   'MaxRange2Mean', maxRange2Mean)
%       varargin    - 'SelectionMethod': the selection method
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'notNaN'        - select any non-NaN value
%                       'maxRange2Mean' - select vales so that the maximum 
%                                           range is within a percentage 
%                                           of the mean
%                   default == 'maxRange2Mean'
%                   - 'MaxRange2Mean': maximum percentage of range versus mean
%                   must be a nonnegative scalar
%                   default == 40%
%                   - 'SaveSheetFlag': whether to save values in a spreadsheet
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveMatFlag': whether to save values in a .mat file
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
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
%       cd/plot_horizontal_line.m
%       cd/struct2arglist.m
%
% Used by:
%       ~/Online_Stats/onlineOmight_interface.m

% File History:
% 2019-05-14 Adapted from onlineOmight_interface.m
% TODO: Finish input parser

%% Hard-coded parameters
validSelectionMethods = {'notNaN', 'maxRange2Mean'};

% Function-specific parameters
valueStr = 'oscDurationSec';
valueLabel = 'Oscillation duration (s)';
avgValueLabel = 'Average oscillation duration of last 5 sweeps (s)';
valueFunc = @compute_oscillation_duration;
xLimits = [];               % auto
yLimits = [];               % auto
xLabel = 'File Name';
yLabel = 'Oscillation Duration (s)';
xTickAngle = 75;
figTitle = 'Oscillation Duration for each file';

% Plot parameters
valueColor = 'b';
marker = 'o';
markerSize = 8;
lineWidth = 2;
maxXTicks = 20;             % maximum number of X ticks

% Analysis parameters
pauseTime = 1;              % 1 second
minFileSizeBytes = 10000;   % 10000 bytes
nLast = 10;
nToAverage = 5;

% Save parameters
sheetType = 'csv';
abfPathStr = 'abfPaths';

%% Default values for optional arguments
avgFuncDefault = [];            % set later
selectionMethodDefault = 'maxRange2Mean';   
                                % select using maxRange2Mean by default
maxRange2MeanDefault = 40;      % range is not more than 40% of mean by default
saveSheetFlagDefault = true;    % save values in a spreadsheet by default
saveMatFlagDefault = false;      % don't save values in a .mat file by default

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

% Add optional inputs to the Input Parser
addOptional(iP, 'avgFunc', avgFuncDefault, ...
    @(x) validateattributes(x, {'function_handle'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SelectionMethod', selectionMethodDefault, ...
    @(x) any(validatestring(x, validSelectionMethods)));
addParameter(iP, 'MaxRange2Mean', maxRange2MeanDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'SaveSheetFlag', saveSheetFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveMatFlag', saveMatFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, valueFunc, varargin{:});
avgFunc = iP.Results.avgFunc;
selectionMethod = validatestring(iP.Results.SelectionMethod, ...
                                    validSelectionMethods);
maxRange2Mean = iP.Results.MaxRange2Mean;
saveSheetFlag = iP.Results.SaveSheetFlag;
saveMatFlag = iP.Results.SaveMatFlag;

% Keep unmatched arguments for the plot() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Set default averaging function
if isempty(avgFunc)
    avgFunc = @(x) compute_phase_average(x, 'NToAverage', nToAverage, ...
                                    'SelectionMethod', selectionMethod, ...
                                    'MaxRange2Mean', maxRange2Mean);
end

%% Select a directory that contains .abf files
abfDir = uigetdir();

%% Initialize variables
% The .abf files already analyzed
abfPathsAnalyzed = {};

% The computed values
values = [];

% The .abf file path previously under acquisition
pathUnderAcquis = '';

% The last file number analyzed
iFile = 0;

%% Create a figure
fig = figure(1);
clf(fig);
title(figTitle);

xlim('manual');
if ~isempty(xLimits)
    xlim(xLimits);
end

if ~isempty(yLimits)
    ylim(yLimits); ylim('manual');
end

if ~isempty(xLabel)
    xlabel(xLabel);
end

if ~isempty(yLabel)
    ylabel(yLabel);
end

if ~isempty(xTickAngle)
    xtickangle(xTickAngle);
end

% Hold on
hold on;

% Initialize lines
lines = gobjects(3, 1);

%% Analyze as long as the figure is open
while ishandle(fig)
    % Get all .abf files from this directory
    [abfFiles, abfPaths] = ...
        all_files('Directory', abfDir, 'Extension', '.abf');

    % Determine which files have not been analyzed
    [abfPathsNotAnalyzed, iPathsNotAnalyzed] = ...
        setdiff(abfPaths, abfPathsAnalyzed);

    % Count the number of files not analyzed
    nFilesNotAnalyzed = numel(abfPathsNotAnalyzed);

    % If there are no files not analyzed, pause and continue
    if nFilesNotAnalyzed == 0
        pause(pauseTime); continue
    end

    % If there is only one file not analyzed, 
    %   assume that it is under acquisition
    if nFilesNotAnalyzed == 1
        % Update the path under acquisition
        pathUnderAcquis = abfPathsNotAnalyzed{1};

        % Pause and continue
        pause(pauseTime); continue
    end

    % Otherwise, there should be at least two files not analyzed
    %   Retrieve the corresponding file structures
    abfFilesNotAnalyzed = abfFiles(iPathsNotAnalyzed);

    % Determine whether there is the file previously under acquisition
    idxToAnalyze = find_in_strings(pathUnderAcquis, abfPathsNotAnalyzed, ...
                                    'ReturnNan', false);

    % If the file is not found, assume that it has been deleted
    %   and choose a new file to be under acquisition
    if isempty(idxToAnalyze)
        % Update the path under acquisition
        pathUnderAcquis = abfPathsNotAnalyzed{1};

        % Pause and continue
        pause(pauseTime); continue
    end

    % Otherwise, analyze the file that was previously under acquisition
    abfFileToAnalyze = abfFilesNotAnalyzed(idxToAnalyze);
    abfPathToAnalyze = abfPathsNotAnalyzed{idxToAnalyze};
    iFile = iFile + 1;

    % Update the files not analyzed
    abfPathsNotAnalyzed = setdiff(abfPathsNotAnalyzed, {pathUnderAcquis});

    % Choose a new file to be under acquisition
    pathUnderAcquis = abfPathsNotAnalyzed{1};

    % Extract the file size
    abfFileSize = abfFileToAnalyze.bytes;
    
    % Only apply the analysis function if a file size threshold is met
    if abfFileSize >= minFileSizeBytes
        % Apply the analysis function
        newValue = valueFunc(abfPathToAnalyze);
    else
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
    newAvg = avgFunc(lastValues);

    % Compute a new limits
    upperLimit = newAvg * (1 + maxRange2Mean / 200);
    lowerLimit = newAvg * (1 - maxRange2Mean / 200);

    % Print new value and average value
    fprintf("%s: %g\n", valueLabel, newValue);
    fprintf("%s: %g\n", avgValueLabel, newAvg);

    % Plot new value
    plot(iFile, newValue, 'Color', valueColor, 'Marker', marker, ...
            'MarkerSize', markerSize, 'LineWidth', lineWidth, ...
            otherArguments{:});

    % Delete previous horizontal line(s) if any
    if any(ishandle(lines))
        delete(lines)
    end

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

    % Plot horizontal line(s)
    %   Note: Do this after x limits is updated
    lines(1, 1) = ...
        plot_horizontal_line(newAvg, 'Color', 'g', 'LineStyle', '--', ...
                                'LineWidth', 2);
    lines(2:3, 1) = ...
        plot_horizontal_line([upperLimit, lowerLimit], ...
                                'Color', 'r', 'LineStyle', '--', ...
                                'LineWidth', 2);

    % Expand y limits
    if ~isempty(yLimits)
        ylim(compute_axis_limits(values, 'y'));
    end

    % Update plot
    drawnow;
end

%% Save active data
if saveSheetFlag || saveMatFlag
    % Create a time stamp
    timeStamp = create_time_stamp('FormatOut', 'yyyymmddTHHMM');

    % Create the file base for output files
    valuesOutFileBase = ...
        fullfile(abfDir, [timeStamp, '_online_values_', valueStr]);
end

if saveSheetFlag
    % Create a spreadsheet path
    valuesSheetPath = [valuesOutFileBase, '.', sheetType];

    % Create a table
    valueTable = table(abfPathsAnalyzed, values, ...
                        'VariableNames', {abfPathStr, valueStr});

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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%