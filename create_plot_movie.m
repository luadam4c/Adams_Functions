function plotFrames = create_plot_movie (figHandle, varargin)
%% Creates a plot movie from a given figure with traces
% Usage: plotFrames = create_plot_movie (figHandle, fiSeconds (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       plotFrames  - frames for the plot movie, with fields
%                       cdata    - RGB intensity data
%                       colormap - color map used
%                       time     - time of the frame in seconds
%                       duration - duration of frame in seconds
%                   specified as a structure array
%
% Arguments:
%       figHandle   - figure handle
%                   must be a numeric array or a cell array of numeric vectors
%       fiSeconds   - (opt) frame intervals in seconds
%                   must be a positive vector
%                   default == 1 s (1 Hz)
%       varargin    - 'FrameTimes': frame times for each sample
%                   must be a numeric vector
%                   default == constructed from fiSeconds
%                   - 'FileBase': file base for the movie
%                   must be a string scalar or a character vector
%                   default == set in write_frames.m
%                   - 'PlotLeadPoints': whether to plot leading points
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'AlignToSelected': whether to align to selected points
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for create_empty_frames()
%
% Requires:
%       cd/align_vectors_by_index.m
%       cd/argfun.m
%       cd/create_empty_frames.m
%       cd/create_error_for_nargin.m
%       cd/extract_fields.m
%       cd/find_custom.m
%       cd/find_first_match.m
%       cd/hold_on.m
%       cd/hold_off.m
%       cd/ismatch.m
%       cd/match_time_info.m
%       cd/set_default_flag.m
%       cd/write_frames.m
%
% Used by:
%       cd/create_plot_movie.m

% File History:
% 2020-05-04 Moved from create_plot_movie.m
% 2020-05-04 Added 'FileBase' as an optional argument
% 2020-05-04 Added 'PlotLeadPoints' as an optional argument
% 2020-05-06 Added 'AlignToSelected' as an optional argument

%% Hard-coded constants
MS_PER_S = 1000;
markerExpansionRatio = 2;

%% Hard-coded parameters

%% Default values for optional arguments
fiSecondsDefault = [];      % set later
frameTimesDefault = [];     % set later
fileBaseDefault = '';       % don't save by default
plotLeadPointsDefault = true;   % plot leading points by default
alignToSelectedDefault = false; % don't align to selected by default

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
addRequired(iP, 'figHandle', ...
    @(x) validateattributes(x, {'matlab.ui.Figure'}, {'2d'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'fiSeconds', fiSecondsDefault, ...
    @(x) assert(isempty(x) || ispositivevector(x), ...
                'fiSeconds must be either empty or a positive vector!'));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FrameTimes', frameTimesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'FileBase', fileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PlotLeadPoints', plotLeadPointsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'AlignToSelected', alignToSelectedDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, figHandle, varargin{:});
fiSeconds = iP.Results.fiSeconds;
frameTimes = iP.Results.FrameTimes;
fileBase = iP.Results.FileBase;
plotLeadPoints = iP.Results.PlotLeadPoints;
alignToSelected = iP.Results.AlignToSelected;

% Keep unmatched arguments for the create_empty_frames() function
otherArguments = iP.Unmatched;

%% Preparation
% Make current figure
figure(figHandle);

% Get the figure position
figPosition = get(figHandle, 'Position');

% Get the figure width and height, rounded to the nearest integer
figWidth = round(figPosition(3));
figHeight = round(figPosition(4));

% Extract the line handles
lineHandles = findobj(figHandle, 'Type', 'Line');

% Deal with selected points
[lineHandles, selectedHandles] = ...
    deal_with_selected(lineHandles, alignToSelected);

% Extract all y data
yData = extract_fields(lineHandles, 'YData', 'UniformOutput', false);

% Find the number of samples for each line object
nSamples = cellfun(@numel, yData);

% Set the number of plot frames to be the maximum number of samples
nPlotFrames = max(nSamples) + 1;

% Set default frame interval
if isempty(fiSeconds) && isempty(frameTimes)
    fiSeconds = 1;
end

% Convert the frame interval to milliseconds
fiMs = fiSeconds * MS_PER_S;

% Match the frame times with the frame interval
[frameTimes, fiMs] = ...
    match_time_info(frameTimes, fiMs, nPlotFrames, 'TimeUnits', 's');

% Convert the matched frame interval to seconds
fiSeconds = fiMs / MS_PER_S;

%% Do the job
% Initialize plot movie frames
plotFrames = create_empty_frames(figHeight, figWidth, [nPlotFrames, 1], ...
                                'Duration', fiSeconds, otherArguments);

% Initialize lead point handles
leadPointHandles = plot([]);

% Find all legends
legends = findobj(gcf, 'type', 'Legend');

% Set auto update to be off
if ~isempty(legends)
    set(legends, 'AutoUpdate', 'off');
end

% Loop through all frame times in reverse
for iPlotFrame = nPlotFrames:-1:1
    % Get the current plot frame time
    plotFrameTimeThis = frameTimes(iPlotFrame);

    % Update plot
    drawnow;

    % Capture this plot frame
    plotFrameThis = getframe(figHandle);

    % Store info in plotFrames array and add time and duration
    plotFrames(iPlotFrame, 1).cdata = plotFrameThis.cdata;
    plotFrames(iPlotFrame, 1).colormap = plotFrameThis.colormap;
    plotFrames(iPlotFrame, 1).time = plotFrameTimeThis;

    % Only remove the last point if not the last frame
    removeLastPoints = set_default_flag([], iPlotFrame < nPlotFrames);

    % Generate the plot for the previous frame
    leadPointHandles = ...
        generate_previous_frame(lineHandles, leadPointHandles, ...
                                selectedHandles, removeLastPoints, ...
                                plotLeadPoints, markerExpansionRatio);
end

% Write frames to a movie file if requested
if ~isempty(fileBase)
    write_frames(plotFrames, 'FileBase', fileBase);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lineHandles, selectedHandles] = ...
                deal_with_selected (lineHandlesOrig, alignToSelected)

% Extract all y data
yData = extract_fields(lineHandlesOrig, 'YData', 'UniformOutput', false);

% Find the number of samples for each object
nSamples = cellfun(@numel, yData);

% If there are any objects with only one data point, make it a selected handle
%   otherwise, make it a line handle
selectedHandles = lineHandlesOrig(nSamples == 1);
lineHandles = lineHandlesOrig(nSamples > 1);

% Align line objects to selected indices by padding data with NaNs
if alignToSelected
    % Match each line handle with one of the selected handles
    %   and find the index of the selected data in the line data
    [lineHandles, idxInLine] = ...
        match_each_line_with_a_selected(lineHandles, selectedHandles);

    % Align line objects to selected indices
    align_lines_to_index(lineHandles, idxInLine);
end

% Match each selected handle with one of the line handles 
%   and pad the selected object
arrayfun(@(s) pad_selected_to_a_line(s, lineHandles), selectedHandles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lineHandles, idxInLine, indMatchedSelected] = ...
                match_each_line_with_a_selected (lineHandles, selectedHandles)
%% Matches each line object with a selected object

% Extract x and y data from line objects
[xLines, yLines] = ...
    extract_fields(lineHandles, {'XData', 'YData'}, 'UniformOutput', false);

% Extract x and y data from selected objects
[xSelected, ySelected] = ...
    extract_fields(selectedHandles, {'XData', 'YData'}, 'UniformOutput', true);

% For each line object, find the first matching selected object by data overlap 
indMatchedSelectedByData = ...
    cellfun(@(x, y) find_custom(ismember(xSelected, x) & ...
                                    ismember(ySelected, y), ...
                                1, 'first', 'ReturnNan', true), ...
            xLines, yLines, 'UniformOutput', true);

% Find the index of the matching selected point in the line data
idxInLine = ...
    cellfun(@(x, y, i) find_selected_index (x, y, xSelected, ySelected, i), ...
            xLines, yLines, num2cell(indMatchedSelectedByData), ...
            'UniformOutput', true);

% Extract color from line and selected objects
[colorLines, colorSelected] = ...
    argfun(@(x) extract_fields(x, 'Color', 'UniformOutput', false), ...
            lineHandles, selectedHandles);

% For each line object, find the first matching selected object by color
indMatchedSelectedByColor = find_first_match(colorLines, colorSelected);

% Replace NaN values by data overlap with the match by color
toReplace = isnan(indMatchedSelectedByData);
indMatchedSelected = indMatchedSelectedByData;
indMatchedSelected(toReplace) = indMatchedSelectedByColor(toReplace);

% Try to replace NaN values by data overlap with the match by color
idxInLine = update_nan_indices(idxInLine, indMatchedSelected);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function idxInLine = find_selected_index (xLine, yLine, ...
                                    xSelectedAll, ySelectedAll, idxSelected)

if isnan(idxSelected)
    idxInLine = NaN;
else
    idxInLine = find(xLine == xSelectedAll(idxSelected) & ...
                    yLine == ySelectedAll(idxSelected), 1, 'first');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function idxInLine = update_nan_indices (idxInLine, indMatchedSelected)

% Find all rows that are NaNs
rowsNans = find(isnan(idxInLine));

% Return if there are none
if isempty(rowsNans)
    return
end

% For each row that is NaN, find the first row 
%   with the same matched selected object that already has idxInLine
rowsSameSelected = ...
    arrayfun(@(x) find_custom(indMatchedSelected == indMatchedSelected(x) & ...
                                ~isnan(idxInLine), 1, 'first', ...
                                'ReturnNan', true), ...
            rowsNans);

% Update the rows that were NaNs
toUpdate = ~isnan(rowsSameSelected);
idxInLine(rowsNans(toUpdate)) = idxInLine(rowsSameSelected(toUpdate));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pad_selected_to_a_line (selectedHandle, lineHandles)
%% Pads a selected object to match the dimensions of a matching line object

% Extract x and y data from the selected object
xSelected = selectedHandle.XData;
ySelected = selectedHandle.YData;

% Extract x and y data from line objects
[xLines, yLines] = ...
    extract_fields(lineHandles, {'XData', 'YData'}, 'UniformOutput', false);

% Find the first matching line object by data overlap
idxMatchedLine = find(ismatch(xLines, xSelected, 'MatchMode', 'parts') & ...
                ismatch(yLines, ySelected, 'MatchMode', 'parts'), 1, 'first');

% If doesn't exist, do nothing
if isempty(idxMatchedLine)
    return
end

% Extract the data for this line object
xThisLine = xLines{idxMatchedLine};
yThisLine = yLines{idxMatchedLine};

% Find the number of samples
nSamplesThisLine = numel(yThisLine);

% Find the index of the matching selected point in the line data
idxSelected = find(xThisLine == xSelected & yThisLine == ySelected, 1, 'first');

% Pad NaNs for selected object
[selectedHandle.XData, selectedHandle.YData] = ...
    argfun(@(x) [nan(1, idxSelected - 1), x, ...
                nan(1, nSamplesThisLine - idxSelected)], ...
            xSelected, ySelected);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function align_lines_to_index(lineHandles, idxSelected)
%% Aligns data in line objects to given indices for each line

% Extract x and y data from line objects
[xDataOrig, yDataOrig] = ...
    extract_fields(lineHandles, {'XData', 'YData'}, 'UniformOutput', false);

% Align vectors to a particular index
[xDataNew, yDataNew] = ...
    argfun(@(x) align_vectors_by_index(x, idxSelected), xDataOrig, yDataOrig);

% Set as new data
cellfun(@(a, b) set(a, 'XData', b), num2cell(lineHandles), xDataNew);
cellfun(@(a, b) set(a, 'YData', b), num2cell(lineHandles), yDataNew);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function leadPointHandles = generate_previous_frame (lineHandles, ...
                                        leadPointHandles, selectedHandles, ...
                                        removeLastPoints, plotLeadPoints, ...
                                        markerExpansionRatio)
%% Modifies the plot for the previous frame

% Remove the last data point of each line and selected object
if removeLastPoints
    arrayfun(@remove_last_data_point, lineHandles);
    arrayfun(@remove_last_data_point, selectedHandles);
end

% Return the last data point remaining in each Line object
[xLast, yLast] = arrayfun(@return_last_data_point, lineHandles, ...
                            'UniformOutput', true);

% Plot over the last data points
if plotLeadPoints
    % Extract the color for each Line object
    lineColors = extract_fields(lineHandles, 'Color', 'UniformOutput', false);

    % Extract the parent axes each Line object
    axHandles = extract_fields(lineHandles, 'Parent', 'UniformOutput', false);

    % Extract the marker sizes for each Line object
    markerSizeOrig = extract_fields(lineHandles, 'MarkerSize', ...
                                    'UniformOutput', true);

    % Decide on new marker size
    markerSize = num2cell(markerSizeOrig .* markerExpansionRatio);

    % Plot or update each of the last data points as large dots
    if isempty(leadPointHandles)
        leadPointHandles = cellfun(@plot_large_dots, axHandles, ...
                                    num2cell(xLast), num2cell(yLast), ...
                                    lineColors, markerSize);
    else
        arrayfun(@update_large_dot, leadPointHandles, xLast, yLast);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function remove_last_data_point (lineHandle)
%% Removes the last data point from the plot and returns the coordinates

% Extract x and y data
xDataOld = lineHandle.XData;
yDataOld = lineHandle.YData;

% Remove and return last data point
if ~isempty(xDataOld) && ~isempty(yDataOld)
    % Remove last data point
    lineHandle.XData = xDataOld(1:end-1);
    lineHandle.YData = yDataOld(1:end-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLast, yLast] = return_last_data_point (lineHandle)
%% Return the coordinates of the last data point from the plot and returns

% Extract x and y data
xDataOld = lineHandle.XData;
yDataOld = lineHandle.YData;

% Return last data point
if ~isempty(xDataOld) && ~isempty(yDataOld)
    xLast = xDataOld(end);
    yLast = yDataOld(end);
else
    xLast = NaN;
    yLast = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lineObject = plot_large_dots (axHandle, x, y, colorMap, ...
                                        markerSize, varargin)

% Get current axes
set(gcf, 'CurrentAxes', axHandle);

% Hold on
wasHold = hold_on;

% Plot the dot
lineObject = plot(axHandle, x, y, 'Color', colorMap, ...
                'LineStyle', 'none', 'Marker', '.', ...
                'MarkerSize', markerSize, varargin{:});

% Hold off
hold_off(wasHold);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_large_dot (leadPointHandle, xLast, yLast)

leadPointHandle.XData = xLast;
leadPointHandle.YData = yLast;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
