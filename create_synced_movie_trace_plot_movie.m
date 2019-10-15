function [plotFrames, handles] = create_synced_movie_trace_plot_movie (frames, data, varargin)
%% Creates a plot movie showing a movie and a trace in synchrony
% Usage: [plotFrames, handles] = create_synced_movie_trace_plot_movie (frames, data, varargin)
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
%       frames      - frames structure, with fields:
%                       cdata    - RGB intensity data
%                       colormap - color map used
%                       time     - time of the frame in seconds
%                       duration - duration of frame in seconds
%                   must be a structure array
%       data        - trace data
%                   must be a numeric array or a cell array of numeric vectors
%       siMs        - (opt) sampling interval in ms
%                   must be a positive vector
%                   default == 5 ms (200 Hz)
%       varargin    - 'TimeVec': original time vector
%                   must be a numeric vector
%                   default == constructed from siMs
%                   - 'TraceLabels': labels for the traces, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == {'Trace #1', 'Trace #2', ...}
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/compute_axis_limits.m
%       cd/create_subplots.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/create_time_vectors.m
%       cd/extract_fields.m
%       cd/force_matrix.m
%       cd/struct2arglist.m
%       cd/match_time_info.m
%       cd/plot_frame.m
%       cd/plot_vertical_line.m
%
% Used by:
%       cd/create_pleth_EEG_movies.m

% File History:
% 2019-09-05 Created by Adam Lu
% 2019-09-06 Added 'TraceLabels' as an optional argument
% 

%% Hard-coded parameters
axesCoverage = 90;

% TODO: Make optional arguments
viewWindowSec = 5;          % width of viewing window in seconds
plotFrameRate = [];         % frame rate for playing the plot movie
yLimits = [];               % y limits for trace subplot
timeLabel = 'Time (s)';

%% Default values for optional arguments
siMsDefault = [];
tVecDefault = [];           % set later
traceLabelsDefault = '';    % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'frames', ...
    @(x) validateattributes(x, {'struct'}, {'vector'}));
addRequired(iP, 'data', ...
    @(x) validateattributes(x, {'cell', 'numeric'}, {'2d'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'siMs', siMsDefault, ...
    @(x) assert(isempty(x) || ispositivevector(x), ...
                'siMs must be either empty or a positive vector!'));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TimeVec', tVecDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'TraceLabels', traceLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));

% Read from the Input Parser
parse(iP, frames, data, varargin{:});
siMs = iP.Results.siMs;
tVec = iP.Results.TimeVec;
traceLabels = iP.Results.TraceLabels;

% Keep unmatched arguments for the TODO() function
% otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Force as a matrix
%   Note: data vectors are assumed to have the same length
data = force_matrix(data);

% Compute the number of samples
nSamples = size(data, 1);

% Compute the number of traces
nTraces = size(data, 2);

% Compute the number of subplot rows
nRowsSubplots = nTraces * 2;

% Compute the subplot grid positions
subPlotGridPositions = [{1:nTraces}, num2cell(nTraces + 1:nRowsSubplots)];

% Create default trace labels
if isempty(traceLabels)
    traceLabels = create_labels_from_numbers(1:nTraces, 'Prefix', 'Trace #');
end

% Construct the time vector if necessary
tVec = match_time_info(tVec, siMs, nSamples);

% Decide on the x tick locations (every second)
firstSecond = round(tVec(1));
lastSecond = round(tVec(end));
xTicks = firstSecond:lastSecond;

% Extract the frame start times
frameTimes = extract_fields(frames, 'time', 'UniformOutput', true);

% Extract the first frame's data
firstFrameData = frames(1).cdata;

% Extract the width and height of the first frame
frameHeight = size(firstFrameData, 1);
frameWidth = size(firstFrameData, 2);

% Decide on the figure position
figPosition = get(0, 'defaultfigureposition');
figPosition(3) = round(frameWidth / (axesCoverage / 100));
figPosition(4) = round(frameHeight * 2 / (axesCoverage / 100));

% % Count the number of frames
% nFrames = numel(frames);

% Compute default frame rate for the plot movie
if isempty(plotFrameRate)
    % Extract the frame durations
    frameDurations = extract_fields(frames, 'duration', 'UniformOutput', true);

    % Get the average frame duration in seconds
    meanFrameDuration = mean(frameDurations);

    % Use the average frame rate of the movie
    plotFrameRate = round(1/meanFrameDuration);
end

% Compute the appropriate y limits
if isempty(yLimits)
    yLimits = arrayfun(@(x) compute_axis_limits(data(:, x), 'y', ...
                                        'Coverage', 90, 'AutoZoom', true), ...
                        transpose(1:nTraces), 'UniformOutput', false);
end

% Compute the plot frame duration in seconds
plotFrameDuration = 1 / plotFrameRate;

% Compute the total plot duration
totalDuration =  tVec(end) - tVec(1);

% Count the number of plot frames
nPlotFrames = floor(totalDuration / plotFrameDuration);

% Construct the time points for each frame in seconds
plotFrameTimes = create_time_vectors(nPlotFrames, 'TimeStart', tVec(1), ...
                    'SamplingIntervalSeconds', plotFrameDuration, ...
                    'TimeUnits', 's', 'BoundaryMode', 'leftadjust');

%% Do the job
% Create subplots
[fig, ax] = create_subplots(nRowsSubplots, 1, subPlotGridPositions, ...
                            'FigPosition', figPosition);

% Get the figure position
figPosition = get(fig, 'Position');

% Get the figure width and height
figWidth = figPosition(3);
figHeight = figPosition(4);

% Initialize plot movie frames
plotFrames = create_empty_frames(figHeight, figWidth, [nPlotFrames, 1], ...
                                'Duration', plotFrameDuration);
                
% Rename the axes
% Plot the trace on the bottom
movieSubPlot = ax(1);
traceSubPlots = ax(2:end);

% Get the initial plot frame time
plotFrameTimeThis = plotFrameTimes(1);

% Compute the first x limits
xLimitsFirst = plotFrameTimeThis + viewWindowSec * 0.5 * [-1, 1];

% TODO: Update plot_traces to use 'AxesHandle' and use it
% Plot the traces
for iTrace = 1:nTraces
    traceLine = plot(traceSubPlots(iTrace), tVec, data(:, iTrace));
end

% Zoom in to the first window, set the axis limits, and remove time ticks
for iTrace = 1:nTraces
    set(traceSubPlots(iTrace), 'XLim', xLimitsFirst, 'XTick', [], ...
        'YLim', yLimits{iTrace}, 'TickLength', [0, 0]);
end

% Link the time axis limits
linkaxes(traceSubPlots, 'x');

% Add time tick locations for the bottom-most plot
set(traceSubPlots(end), 'XTick', xTicks);

% Create time label for the bottom-most plot
traceSubPlots(end).XLabel.String = timeLabel;

% Create data labels
for iTrace = 1:nTraces
    traceSubPlots(iTrace).YLabel.String = traceLabels{iTrace};
end

% Plot a vertical line through the plots
vertLines = gobjects(nTraces, 1);
for iTrace = 1:nTraces
    vertLines(iTrace) = ...
        plot_vertical_line(plotFrameTimeThis, 'LineStyle', '--', ...
                                    'LineWidth', 2, 'ColorMap', 'r', ...
                                    'AxesHandle', traceSubPlots(iTrace));
end

% Look for the corresponding movie frame
[iFrameThis, frameTimeThis] = find_nearest_frame(plotFrameTimeThis, frameTimes);

% Extract this frame
frameThis = frames(iFrameThis);

% Plot the frame on the top
handles = plot_frame(frameThis, 'AxesHandle', movieSubPlot, ...
                                'AxesCoverage', axesCoverage);
im = handles.im;

% Capture this plot frame
plotFrameThis = getframe(fig);

% Store info in plotFrames array and add time and duration
plotFrames(1, 1).cdata = plotFrameThis.cdata;
plotFrames(1, 1).colormap = plotFrameThis.colormap;
plotFrames(1, 1).time = plotFrameTimeThis;

% Loop through all frame times
if nPlotFrames > 1
    % Prevent y ticks from automatically updating
    % TODO: How to do this?
    % NOT this: yticklabels('manual');
	% NOT this: set(traceSubPlot, 'YTickLabelMode', 'manual');
	% NOT this: set(traceSubPlot, 'DataAspectRatioMode', 'manual');
    
    for iPlotFrame = 2:nPlotFrames
        % Get the current plot frame time
        plotFrameTimeThis = plotFrameTimes(iPlotFrame);

        % Update the zoom window of the first subplot
        %   Note: the other trace subplots are linked so will be updated too
        set(traceSubPlots(1), 'XLim', plotFrameTimeThis + viewWindowSec * 0.5 * [-1, 1]);

        % Update the vertical line position
        for iTrace = 1:nTraces
            set(vertLines(iTrace), 'XData', [plotFrameTimeThis, plotFrameTimeThis]);
        end
        
        % Update the corresponding movie frame if necessary
        if plotFrameTimeThis < frameTimeThis || ...
                plotFrameTimeThis >= frameTimes(iFrameThis + 1)
            % Look for the corresponding movie frame
            [iFrameThis, frameTimeThis] = ...
                find_nearest_frame(plotFrameTimeThis, frameTimes);

            % Update the frame data
            set(im, 'CData', frames(iFrameThis).cdata);
        end
        
        % Capture this plot frame
        plotFrameThis = getframe(fig);

        % Store info in plotFrames array and add time and duration
        plotFrames(iPlotFrame, 1).cdata = plotFrameThis.cdata;
        plotFrames(iPlotFrame, 1).colormap = plotFrameThis.colormap;
        plotFrames(iPlotFrame, 1).time = plotFrameTimeThis;
    end
end

%% Output results
% Save handles
handles.fig = fig;
handles.ax = ax;
handles.traceLine = traceLine;
handles.vertLines = vertLines;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [iFrame, startTime] = find_nearest_frame (t, allStartTimes)

% Compute all time differences
timeDiffs = t - allStartTimes;

% Ignore all negative entries
timeDiffs(timeDiffs < 0) = NaN;

% Take the minimum over all time differences
[startTime, iFrame] = min(timeDiffs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

plotFrameTimes = tVec(1) + transpose(0:(nPlotFrames - 1)) * plotFrameDuration;

% Display data aspect ratio
disp(get(traceSubPlot, 'PlotBoxAspectRatio'));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
