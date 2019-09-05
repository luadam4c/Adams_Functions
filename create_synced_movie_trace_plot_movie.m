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
%                   must be a numeric vector
%       siMs        - (opt) sampling interval in ms
%                   must be a positive vector
%                   default == 5 ms (200 Hz)
%       varargin    - 'TimeVec': original time vector
%                   must be a numeric vector
%                   default == constructed from siMs
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/compute_axis_limits.m
%       cd/create_subplots.m
%       cd/create_error_for_nargin.m
%       cd/create_time_vectors.m
%       cd/extract_fields.m
%       cd/struct2arglist.m
%       cd/match_time_info.m
%       cd/plot_frame.m
%       cd/plot_vertical_line.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-05 Created by Adam Lu
% 

%% Hard-coded parameters
% TODO: Make optional arguments
viewWindowSec = 5;          % width of viewing window in seconds
plotFrameRate = [];         % frame rate for playing the plot movie
yLimits = [];               % y limits for trace subplot
timeLabel = 'Time (s)';
dataLabel = 'EEG amplitude (uV)';

%% Default values for optional arguments
siMsDefault = [];
tVecDefault = [];           % set later

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
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'siMs', siMsDefault, ...
    @(x) assert(isempty(x) || ispositivevector(x), ...
                'siMs must be either empty or a positive vector!'));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TimeVec', tVecDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Read from the Input Parser
parse(iP, frames, data, varargin{:});
siMs = iP.Results.siMs;
tVec = iP.Results.TimeVec;

% Keep unmatched arguments for the TODO() function
% otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Compute the number of samples
nSamples = numel(data);

% Construct the time vector if necessary
tVec = match_time_info(tVec, siMs, nSamples);

% Extract the frame start times
frameTimes = extract_fields(frames, 'time', 'UniformOutput', true);

% Extract the first frame's data
firstFrameData = frames(1).cdata;

% Extract the width and height of the first frame
frameHeight = size(firstFrameData, 1);
frameWidth = size(firstFrameData, 2);

% Decide on the first subplot position
firstSubplotPosition = get(0, 'defaultfigureposition');
firstSubplotPosition(3) = frameWidth;
firstSubplotPosition(4) = frameHeight;

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
    yLimits = compute_axis_limits(data, 'y', 'Coverage', 90, 'AutoZoom', true);
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
[fig, ax] = create_subplots(2, 1, 'CenterPosition', firstSubplotPosition);

% Get the figure position
figPosition = get(fig, 'Position');

% Get the figure width and height
figWidth = figPosition(3);
figHeight = figPosition(4);

% Initialize plot movie frames
plotFrames = create_empty_frames(figHeight, figWidth, [nPlotFrames, 1]);
                
% Rename the axes
% Plot the trace on the bottom
movieSubPlot = ax(1);
traceSubPlot = ax(2);

% Get the initial plot frame time
plotFrameTimeThis = plotFrameTimes(1);

% Compute the first x limits
xLimitsFirst = plotFrameTimeThis + viewWindowSec * 0.5 * [-1, 1];

% TODO: Update plot_traces to use 'AxesHandle' and use it
% Plot the trace
traceLine = plot(traceSubPlot, tVec, data);

% Zoom in to the first window and set the y limits
set(traceSubPlot, 'XLim', xLimitsFirst, 'YLim', yLimits);

% Create time label
xlabel(timeLabel);

% Create data label
ylabel(dataLabel);

% Plot a vertical line through the plot
vertLine = plot_vertical_line(plotFrameTimeThis, 'LineStyle', '--', ...
                                'LineWidth', 2, 'ColorMap', 'r');

% Look for the corresponding movie frame
[iFrameThis, frameTimeThis] = find_nearest_frame(plotFrameTimeThis, frameTimes);

% Extract this frame
frameThis = frames(iFrameThis);

% Plot the frame on the top
handles = plot_frame(frameThis, 'AxesHandle', movieSubPlot);
im = handles.im;

% Capture this plot frame
plotFrameThis = getframe(fig);

% Store info in plotFrames array and add time and duration
plotFrames(1, 1).cdata = plotFrameThis.cdata;
plotFrames(1, 1).colormap = plotFrameThis.colormap;
plotFrames(1, 1).time = plotFrameTimeThis;
plotFrames(1, 1).duration = plotFrameDuration;

% Loop through all frame times
if nPlotFrames > 1
    for iPlotFrame = 2:nPlotFrames
        % Get the current plot frame time
        plotFrameTimeThis = plotFrameTimes(iPlotFrame);

        % Update the zoom window
        set(traceSubPlot, 'XLim', plotFrameTimeThis + viewWindowSec * 0.5 * [-1, 1]);

        % Update the vertical line position
        set(vertLine, 'XData', [plotFrameTimeThis, plotFrameTimeThis]);

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
        plotFrames(iPlotFrame, 1).duration = plotFrameDuration;
    end
end

%% Output results
% Save handles
handles.fig = fig;
handles.ax = ax;
handles.traceLine = traceLine;
handles.vertLine = vertLine;

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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
