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
%       siSeconds   - (opt) sampling interval in seconds
%                   must be a positive vector
%                   default == 0.005 s (200 Hz)
%       varargin    - 'TimeVec': original time vector in seconds
%                   must be a numeric vector
%                   default == constructed from siSeconds
%                   - 'TraceLabels': labels for the traces, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == {'Trace #1', 'Trace #2', ...}
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/compute_axis_limits.m
%       cd/compute_spectrogram.m
%       cd/create_subplots.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/create_time_vectors.m
%       cd/extract_fields.m
%       cd/find_in_strings.m
%       cd/force_column_cell.m
%       cd/force_matrix.m
%       cd/struct2arglist.m
%       cd/match_time_info.m
%       cd/plot_frame.m
%       cd/plot_spectrogram.m
%       cd/plot_vertical_line.m
%
% Used by:
%       cd/create_pleth_EEG_movies.m

% File History:
% 2019-09-05 Created by Adam Lu
% 2019-09-06 Added 'TraceLabels' as an optional argument
% 2019-10-15 Added plotSpectrogram
% 

%% Hard-coded constants
MS_PER_S = 1000;

%% Hard-coded parameters
axesCoverage = 90;

% TODO: Make optional arguments
viewWindowSec = 5;          % width of viewing window in seconds
plotFrameRate = [];         % frame rate for playing the plot movie
yLimits = [];               % y limits for trace subplot
yLimitsSpect = [1, 50];
timeLabel = 'Time (s)';
plotSpectrogram = [];       % set later

%% Default values for optional arguments
siSecondsDefault = 0.005;
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
addOptional(iP, 'siSeconds', siSecondsDefault, ...
    @(x) assert(isempty(x) || ispositivevector(x), ...
                'siSeconds must be either empty or a positive vector!'));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TimeVec', tVecDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'TraceLabels', traceLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));

% Read from the Input Parser
parse(iP, frames, data, varargin{:});
siSeconds = iP.Results.siSeconds;
tVec = iP.Results.TimeVec;
traceLabels = iP.Results.TraceLabels;

% Keep unmatched arguments for the TODO() function
% otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Force as a matrix
%   Note: data vectors are assumed to have the same length
data = force_matrix(data);

% Determine the index for the first EEG trace
%   Note: this will be used to plot a spectrogram
if ~isempty(traceLabels)
    idxEEG = find_in_strings('EEG', traceLabels, 'MaxNum', 1);
else
    idxEEG = [];
end

% Determine whether to plot spectrogram
if isempty(idxEEG)
    plotSpectrogram = false;
elseif isempty(plotSpectrogram)
    plotSpectrogram = true;
end

% Compute the number of samples
nSamples = size(data, 1);

% Compute the number of traces
nTraces = size(data, 2);

% Compute the number of plots
if plotSpectrogram
    nPlots = nTraces + 1;
else
    nPlots = nTraces;
end

% Compute the number of subplot rows
nRowsSubplots = nPlots * 2;

% Compute the subplot grid positions
subPlotGridPositions = [{1:nPlots}, num2cell(nPlots + 1:nRowsSubplots)];

% Create default trace labels
if isempty(traceLabels)
    traceLabels = create_labels_from_numbers(1:nTraces, 'Prefix', 'Trace #');
else
    traceLabels = force_column_cell(traceLabels);
end

% Update trace labels if necessary
if plotSpectrogram && numel(traceLabels) < nPlots
    % Update trace labels
    traceLabels = [traceLabels; {'Spect'}];
end

% Construct the time vector if necessary
siMs = siSeconds * MS_PER_S;
[tVec, siMs] = match_time_info(tVec, siMs, nSamples, 'TimeUnits', 's');

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
if plotSpectrogram && numel(yLimits) < nPlots
    yLimits = [yLimits; {yLimitsSpect}];

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

% Only compute spectrogram if to be plotted
if plotSpectrogram
    % Get the entire EEG trace
    eegChannelValues = data(:, idxEEG);

    % Compute sampling interval in seconds
    siSeconds = siMs / MS_PER_S;

    % Compute the spectrogram
    [spectData, freqHz, timeInstantsSeconds] = ...
        compute_spectrogram(eegChannelValues, siSeconds);
end

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

% TODO: Update plot_traces.m with plotSpectrogramFlag?
% Plot the traces
traceLines = gobjects(nTraces, 1);
spect = gobjects;
for iPlot = 1:nPlots
    % Plot the appropriate trace or map
    if plotSpectrogram && iPlot == nPlots
        spect = plot_spectrogram(spectData, timeInstantsSeconds, freqHz, ...
                            'AxesHandle', traceSubPlots(iPlot));
    else
        traceLines(iPlot) = plot(traceSubPlots(iPlot), tVec, data(:, iPlot));
    end
end

% Zoom in to the first window, set the axis limits, and remove time ticks
for iPlot = 1:nPlots
    set(traceSubPlots(iPlot), 'XLim', xLimitsFirst, 'XTick', [], ...
        'YLim', yLimits{iPlot}, 'TickLength', [0, 0]);
end

% Link the time axis limits
linkaxes(traceSubPlots, 'x');

% Add time tick locations for the bottom-most plot
set(traceSubPlots(end), 'XTick', xTicks);

% Create time label for the bottom-most plot
traceSubPlots(end).XLabel.String = timeLabel;

% Create data labels
for iPlot = 1:nPlots
    traceSubPlots(iPlot).YLabel.String = traceLabels{iPlot};
end

% Plot a vertical line through the plots
vertLines = gobjects(nPlots, 1);
for iPlot = 1:nPlots
    vertLines(iPlot) = ...
        plot_vertical_line(plotFrameTimeThis, 'LineStyle', '--', ...
                                    'LineWidth', 2, 'ColorMap', 'r', ...
                                    'AxesHandle', traceSubPlots(iPlot));
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
        set(traceSubPlots(1), 'XLim', ...
                            plotFrameTimeThis + viewWindowSec * 0.5 * [-1, 1]);

        % Update the vertical line position
        arrayfun(@(x) set(x, 'XData', plotFrameTimeThis * [1, 1]), vertLines);
        
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
handles.traceLines = traceLines;
handles.spect = spect;
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
