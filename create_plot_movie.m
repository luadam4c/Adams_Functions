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
%                   - Any other parameter-value pair for plot_traces()
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/extract_fields.m
%       cd/hold_on.m
%       cd/hold_off.m
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

%% Hard-coded constants
MS_PER_S = 1000;
markerExpansionRatio = 2;

%% Hard-coded parameters
validPlotModes = {'overlapped', 'parallel', 'staggered'};

%% Default values for optional arguments
fiSecondsDefault = [];      % set later
frameTimesDefault = [];     % set later
fileBaseDefault = '';       % don't save by default
plotLeadPointsDefault = true;   % plot leading points by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
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

% Read from the Input Parser
parse(iP, figHandle, varargin{:});
fiSeconds = iP.Results.fiSeconds;
frameTimes = iP.Results.FrameTimes;
fileBase = iP.Results.FileBase;
plotLeadPoints = iP.Results.PlotLeadPoints;

% Keep unmatched arguments for the plot_traces() function
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
                                'Duration', fiSeconds);

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
                        removeLastPoints, plotLeadPoints, markerExpansionRatio);
end

% Write frames to a movie file if requested
if ~isempty(fileBase)
    write_frames(plotFrames, 'FileBase', fileBase);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function leadPointHandles = generate_previous_frame (lineHandles, ...
                                leadPointHandles, removeLastPoints, ...
                                plotLeadPoints, markerExpansionRatio)
%% Modifies the plot for the previous frame

% Remove the last data point of each Line object
if removeLastPoints
    arrayfun(@remove_last_data_point, lineHandles);
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
