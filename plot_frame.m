function handles = plot_frame (frame, varargin)
%% Plots a specific movie frame
% Usage: handles = plot_frame (frame, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       xylo = read_frames('xylophone.mp4');
%       figure(1); plot_frame(xylo(4))
%       xylo4 = read_frames('xylophone.mp4', 'FrameIndex', 4);
%       figure(2); plot_frame(xylo4)
%
% Outputs:
%       handles     - handles, with fields:
%                       fig - Figure object that frame was plotted on
%                       ax  - Axes object that frame was plotted on
%                       im  - Image object for the frame
%                   specified as a structure
%
% Arguments:
%       frame       - MATLAB movie frame structure, with fields:
%                       cdata    - RGB intensity data
%                       colormap - color map used
%                   specified as a scalar structure
%       plotFunc    - (opt) plotting function used
%                   must be a function handle
%                   default == @imshow
%       varargin    - 'FigHandle': figure handle for created figure
%                   must be a empty or a figure object handle
%                   default == []
%                   - 'AxesHandle': axes handle for created axes
%                   must be a empty or a axes object handle
%                   default == []
%                   - 'VideoObject': VideoReader object for the video
%                   must be a VideoReader object
%                   default == VideoReader.empty
%                   - Any other parameter-value pair for plotFunc()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/set_axes_properties.m
%       cd/set_figure_properties.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/create_synced_movie_trace_plot_movie.m

% File History:
% 2019-09-04 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
plotFuncDefault = function_handle.empty;
figHandleDefault = [];          % no existing figure by default
axHandleDefault = [];           % no existing axes by default
videoObjectDefault = VideoReader.empty;

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
addRequired(iP, 'frame', ...
    @(x) validateattributes(x, {'struct'}, {'2d'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'plotFunc', plotFuncDefault, ...
    @(x) validateattributes(x, {'function_handle'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FigHandle', figHandleDefault);
addParameter(iP, 'AxesHandle', axHandleDefault);
addParameter(iP, 'VideoObject', videoObjectDefault, ...
    @(x) validateattributes(x, {'VideoReader'}, {'2d'}));

% Read from the Input Parser
parse(iP, frame, varargin{:});
figHandle = iP.Results.FigHandle;
axHandle = iP.Results.AxesHandle;
plotFunc = iP.Results.plotFunc;
vidObj = iP.Results.VideoObject;

% Keep unmatched arguments for the plotFunc() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Set default plotting function

% Grab the image data
imData = frame.cdata;

% Decide on the height and width
if ~isempty(vidObj)
    % Read the height and width of the video object
    vidHeight = vidObj.Height;
    vidWidth = vidObj.Width;
else
    % Assume each element is a pixel
    vidHeight = size(imData, 1);
    vidWidth = size(imData, 2);
end

% Set figure properties and retrieve handles
if isempty(findall(0, 'type', 'figure'))
    % Create a figure and update the figure height and width
    fig = set_figure_properties('FigHandle', figHandle, ...
                            'Height', vidHeight, 'Width', vidWidth);
else
    % Decide on the figure
    fig = set_figure_properties('FigHandle', figHandle);
end

% Decide on the axes
ax = set_axes_properties('AxesHandle', axHandle);

%% Do the job
% Add other arguments
if isempty(plotFunc)
    plotFunc = @(x) imshow(x, 'Border', 'tight', 'Parent', ax, ...
                            'InitialMagnification', 'fit', otherArguments{:});
else
    plotFunc = @(x) plotFunc(x, 'Parent', ax, otherArguments{:});
end

% Plot the image
im = plotFunc(frame.cdata);

% Remove ticks
set(ax, 'XTick', []);
set(ax, 'YTick', []);

%% Output results
handles.fig = fig;
handles.ax = ax;
handles.im = im;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%