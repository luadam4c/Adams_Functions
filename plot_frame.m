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
%       varargin    - Any other parameter-value pair for plotFunc()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/set_axes_properties.m
%       cd/set_figure_properties.m
%       cd/struct2arglist.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-04 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
plotFuncDefault = function_handle.empty;

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
addOptional(iP, 'plotFunc', plotFuncDefault);
    @(x) validateattributes(x, {'function_handle'}, {'2d'});

% Read from the Input Parser
parse(iP, frame, varargin{:});
plotFunc = iP.Results.plotFunc;

% Keep unmatched arguments for the plotFunc() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Set default plotting function
if isempty(plotFunc)
    plotFunc = @imshow;
end

% Set figure properties and retrieve the figure handle
if ~isempty(vidObj)
    % Read the height and width of the video object
    vidHeight = vidObj.Height;
    vidWidth = vidObj.Width;

    fig = set_figure_properties('Height', vidHeight, 'Width', vidWidth]);
else
    fig = set_figure_properties;
end

%% Do the job
% Add other arguments
plotFuncWithOther = @(x) plotFunc(x, otherArguments{:});

% Plot the image
im = plotFuncWithOther(frame.cdata);

% Retrieve the axes handle
ax = set_axes_properties;

%% Output results
handles.fig = fig;
handles.ax = ax;
handles.im = im;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%