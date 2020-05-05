function [plotFrames, handles] = create_trace_plot_movie (xData, yData, varargin)
%% Creates a plot movie for traces
% Usage: [plotFrames, handles] = create_trace_plot_movie (xData, yData, fiSeconds (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [plotFrames, handles] = create_trace_plot_movie(1:3, magic(3), 'FileBase', 'magic3');
%       dataWithNaN = magic(7);
%       dataWithNaN([1, 4, 7], [2, 4]) = NaN;
%       [plotFrames, handles] = create_trace_plot_movie(1:7, dataWithNaN, 'PlotMode', 'overlapped');
%       write_frames(plotFrames, 'FileBase', 'magic7withNaN');
%       [plotFrames, handles] = create_trace_plot_movie(1:100, randn(100, 3), 0.1);
%       write_frames(plotFrames, 'FileBase', 'randn100');
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
%       xData       - trace x data
%                   must be a numeric array or a cell array of numeric vectors
%       yData       - trace y data
%                   must be a numeric array or a cell array of numeric vectors
%       fiSeconds   - (opt) frame intervals in seconds
%                   must be a positive vector
%                   default == 1 s (1 Hz)
%       varargin    - 'FrameTimes': frame times for each sample
%                   must be a numeric vector
%                   default == constructed from fiSeconds
%                   - 'PlotMode': plotting mode for multiple traces
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'parallel'      - in parallel in subPlots
%                       'overlapped'    - overlapped in a single plot
%                       'staggered'     - staggered in a single plot
%                   default == 'parallel'
%                   - 'FileBase': file base for the movie
%                   must be a string scalar or a character vector
%                   default == set in write_frames.m
%                   - Any other parameter-value pair for plot_traces()
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/create_plot_movie.m
%       cd/force_data_as_matrix.m
%       cd/match_format_vector_sets.m
%       cd/plot_traces.m
%
% Used by:

% File History:
% 2020-05-03 Modified from create_synced_movie_trace_plot_movie.m
% 2020-05-04 Added 'FileBase' as an optional argument

%% Hard-coded constants
MS_PER_S = 1000;

%% Hard-coded parameters
validPlotModes = {'overlapped', 'parallel', 'staggered'};

%% Default values for optional arguments
fiSecondsDefault = [];      % set later
frameTimesDefault = [];     % set later
plotModeDefault = 'parallel'; % plot traces in parallel by default
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
addRequired(iP, 'xData', ...
    @(x) validateattributes(x, {'cell', 'numeric'}, {'2d'}));
addRequired(iP, 'yData', ...
    @(x) validateattributes(x, {'cell', 'numeric'}, {'2d'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'fiSeconds', fiSecondsDefault, ...
    @(x) assert(isempty(x) || ispositivevector(x), ...
                'fiSeconds must be either empty or a positive vector!'));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FrameTimes', frameTimesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'PlotMode', plotModeDefault, ...
    @(x) any(validatestring(x, validPlotModes)));
addParameter(iP, 'FileBase', fileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, xData, yData, varargin{:});
fiSeconds = iP.Results.fiSeconds;
frameTimes = iP.Results.FrameTimes;
plotMode = validatestring(iP.Results.PlotMode, validPlotModes);
fileBase = iP.Results.FileBase;

% Keep unmatched arguments for the plot_traces() function
otherArguments = iP.Unmatched;

%% Preparation
% Match vectors
[xData, yData] = match_format_vector_sets(xData, yData, 'MatchVectors', true);

% Force as a matrix
[xData, yData] = ...
    argfun(@(a) force_data_as_matrix(a, 'ForceVectorAsColumns', true), ...
            xData, yData);

%% Do the job
% Plot everything
handles = plot_traces(xData, yData, 'PlotMode', plotMode, otherArguments);

% Extract figure handle used
figHandle = handles.fig;

% Create the plot movie
plotFrames = create_plot_movie(figHandle, fiSeconds, 'FileBase', fileBase, ...
                                'FrameTimes', frameTimes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Extract the current x and y data, treating each sample point as
%   its own vector
[xDataThis, yDataThis] = ...
    argfun(@(a) num2cell(a(iPlotFrame, :)), xData, yData);

% Make sure to only update plots for all frames except the first
if iPlotFrame > 1
    otherArguments.PlotOnly = true;
end

% Update plots with new data point
handles = plot_traces(xDataThis, yDataThis, 'AxesHandles', axHandles, ...
                        'PlotMode', plotMode, otherArguments);

lineHandles = force_column_cell(handles.plotsData);
cellfun(@remove_last_data_point, lineHandles, 'UniformOutput', false);

% Create subplots if not provided
if isempty(axHandles)
    [figHandle, axHandles] = create_subplots(nTraces, 1, 'AlwaysNew', true);
else
    if numel(axHandles) ~= nTraces
        error('Number of axes provided must match the number of traces!');
    end
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
