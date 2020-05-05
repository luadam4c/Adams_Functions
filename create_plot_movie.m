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
%                   - Any other parameter-value pair for plot_traces()
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/extract_fields.m
%       cd/match_time_info.m
%       cd/write_frames.m
%
% Used by:
%       cd/create_plot_movie.m

% File History:
% 2020-05-04 Moved from create_plot_movie.m
% 2020-05-04 Added 'FileBase' as an optional argument

%% Hard-coded constants
MS_PER_S = 1000;

%% Hard-coded parameters
validPlotModes = {'overlapped', 'parallel', 'staggered'};

%% Default values for optional arguments
fiSecondsDefault = [];      % set later
frameTimesDefault = [];     % set later
fileBaseDefault = '';       % don't save by default

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

% Read from the Input Parser
parse(iP, figHandle, varargin{:});
fiSeconds = iP.Results.fiSeconds;
frameTimes = iP.Results.FrameTimes;
fileBase = iP.Results.FileBase;

% Keep unmatched arguments for the plot_traces() function
otherArguments = iP.Unmatched;

%% Preparation
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
nPlotFrames = max(nSamples);

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

% Loop through all frame times in reverse
for iPlotFrame = nPlotFrames:-1:1
    % Update plot
    drawnow;

    % Get the current plot frame time
    plotFrameTimeThis = frameTimes(iPlotFrame);

    % Capture this plot frame
    plotFrameThis = getframe(figHandle);

    % Store info in plotFrames array and add time and duration
    plotFrames(iPlotFrame, 1).cdata = plotFrameThis.cdata;
    plotFrames(iPlotFrame, 1).colormap = plotFrameThis.colormap;
    plotFrames(iPlotFrame, 1).time = plotFrameTimeThis;

    % Remove last data point
    arrayfun(@remove_last_data_point, lineHandles, 'UniformOutput', false);
end

% Write frames to a movie file if requested
if ~isempty(fileBase)
    write_frames(plotFrames, 'FileBase', fileBase);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lineHandle = remove_last_data_point (lineHandle)
%% Removes the last data point from the plot

% Extract x and y data
xDataOld = lineHandle.XData;
yDataOld = lineHandle.YData;

% Remove last data point
if ~isempty(xDataOld) && ~isempty(yDataOld)
    lineHandle.XData = xDataOld(1:end-1);
    lineHandle.YData = yDataOld(1:end-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
