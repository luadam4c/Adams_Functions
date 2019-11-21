function handles = plot_spectrogram (spectData, timeInstants, freqHz, varargin)
%% Plots the spectrogram data computed by compute_spectrogram.m
% Usage: handles = plot_spectrogram (spectData, timeInstants, freqHz, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [data, freq, time] = compute_spectrogram(rand(1000, 1), 0.005);
%       plot_spectrogram(data, time, freq);
%
% Outputs:
%       handles     - handles
%                       im
%                       ax
%                   specified as a structure
%
% Arguments:
%       spectData       - TODO: Description of spectData
%                       must be a TODO
%       timeInstants    - TODO: Description of timeInstants
%                       must be a TODO
%       freqHz      - TODO: Description of freqHz
%                   must be a TODO
%       varargin    - 'AxesHandle': axes handle for created axes
%                   must be a empty or a axes object handle
%                   default == []
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to empty
%                   must be a 2-element increasing numeric vector
%                   default == uses compute_axis_limits.m
%                   - 'ColorMap': a color map
%                   must be a numeric array with 3 columns
%                   default == '/media/adamX/Settings_Matlab/spectrogram_colormap.mat'
%                   - Any other parameter-value pair for imagesc()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/decide_on_colormap.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/plot_traces_spike2_mat.m

% File History:
% 2019-10-15 Moved from plot_traces_spike2_mat.m
% 

%% Hard-coded parameters
spectColorMapFile = '/media/adamX/Settings_Matlab/spectrogram_colormap.mat';

%% Default values for optional arguments
axHandleDefault = [];           % gca by default
yLimitsDefault = [1, 50];       % Look at 1-50 Hz by default
colorMapDefault = [];           % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
% TODO: validation functions
addRequired(iP, 'spectData');
addRequired(iP, 'timeInstants');
addRequired(iP, 'freqHz');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'AxesHandle', axHandleDefault);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'ColorMap', colorMapDefault);

% Read from the Input Parser
parse(iP, spectData, timeInstants, freqHz, varargin{:});
axHandle = iP.Results.AxesHandle;
yLimits = iP.Results.YLimits;
colorMap = iP.Results.ColorMap;

% Keep unmatched arguments for the imagesc() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Decide on the axes
ax = set_axes_properties('AxesHandle', axHandle);

% Decide on the color map
if isempty(colorMap)
    colorMapFile = matfile(spectColorMapFile);
    colorMap = colorMapFile.colorMap;
else
    colorMap = decide_on_colormap(colorMap);
end

%% Do the job
% Set the colormap
colormap(ax, colorMap);

% Plot the spectrogram
im = imagesc(ax, timeInstants, freqHz, abs(spectData), otherArguments{:});

% Flip the Y Axis so lower frequencies are at the bottom
set(ax, 'YDir', 'normal');

% Restrict to certain y axis limits
if ~isempty(yLimits)
    set(ax, 'YLim', yLimits);
end

%% Outputs
handles.im = im;
handles.ax = ax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%