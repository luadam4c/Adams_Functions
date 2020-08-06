function handles = plot_spectrogram (spectData, timeInstants, freqHz, varargin)
%% Plots the magnitude of the spectrogram data computed by compute_spectrogram.m
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
%                   default == [1, 50]
%                   - 'CLimits': limits of color axis, 
%                               suppress by setting value to empty
%                   must be a 2-element increasing numeric vector
%                   default == set by built-in imagesc()
%                   - 'XLabel': label for the time axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == 'Time (s)'
%                   - 'YLabel': label(s) for the y axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == 'Frequency (Hz)'
%                   - 'ColorMap': a color map
%                   must be a numeric array with 3 columns
%                   default == '/media/adamX/Settings_Matlab/spectrogram_colormap.mat'
%                   - 'PlotOnly': whether to plot spectrogram only
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotColorBar': whether to plot color bar
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for imagesc()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/decide_on_colormap.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/plot_spectrogram_multiunit.m
%       cd/plot_traces_spike2_mat.m

% File History:
% 2019-10-15 Moved from plot_traces_spike2_mat.m
% 2019-12-11 Added color axis limits
% 2020-08-05 Added 'PlotColorBar' as an optional argument
% TODO: Default x label, y label, title and modify effects of PlotOnly flag
% TODO: Allow option to plot the following:
%           (1) FFT amplitude
%           (2) FFT power
%           (3) power spectral density (make this default)
%           (4) decibel

%% Hard-coded parameters
spectColorMapFile = '/media/adamX/Settings_Matlab/spectrogram_colormap.mat';

%% Default values for optional arguments
axHandleDefault = [];           % gca by default
yLimitsDefault = [1, 50];       % Look at 1-50 Hz by default
cLimitsDefault = [];            % use imagesc() defaults
xLabelDefault = '';             % set later
yLabelDefault = '';             % set later
colorMapDefault = [];           % set later
plotOnlyDefault = false;        % plots labels and annotations by default
plotColorBarDefault = true;     % plots color bar by default

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
    @(x) isempty(x) || isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'CLimits', cLimitsDefault, ...
    @(x) isempty(x) || isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'XLabel', xLabelDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ColorMap', colorMapDefault);
addParameter(iP, 'PlotOnly', plotOnlyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotColorBar', plotColorBarDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, spectData, timeInstants, freqHz, varargin{:});
axHandle = iP.Results.AxesHandle;
yLimits = iP.Results.YLimits;
cLimits = iP.Results.CLimits;
xLabel = iP.Results.XLabel;
yLabel = iP.Results.YLabel;
colorMap = iP.Results.ColorMap;
plotOnly = iP.Results.PlotOnly;
plotColorBar = iP.Results.PlotColorBar;

% Keep unmatched arguments for the imagesc() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Modify flags if necessary
if plotOnly
    plotColorBar = false;
    xLabel = 'suppress';
    yLabel = 'suppress';
end

% Set the default x-axis label
if isempty(xLabel)
    xLabel = 'Time (s)';
end

% Set the default y-axis label
if isempty(yLabel)
    yLabel = 'Frequency (Hz)';
end

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
% Hold on
wasHold = hold_on;

% Set the colormap
colormap(ax, colorMap);

% Plot the spectrogram (magnitude of FFT)
if ~isempty(cLimits)
    im = imagesc(ax, timeInstants, freqHz, abs(spectData), ...
                cLimits, otherArguments{:});
else
    im = imagesc(ax, timeInstants, freqHz, abs(spectData), otherArguments{:});
end

% Flip the Y Axis so lower frequencies are at the bottom
set(ax, 'YDir', 'normal');

% Restrict to certain y axis limits
if ~isempty(yLimits)
    set(ax, 'YLim', yLimits);
end

% Generate an x-axis label
if ~strcmpi(xLabel, 'suppress')
    xlabel(xLabel);
end

% Generate a y-axis label
if ~strcmpi(yLabel, 'suppress')
    ylabel(yLabel);
end

% Plot a color bar
if plotColorBar
    colorBar = colorbar;
    colorBar.Label.String = '(V \cdot s)';
end

% Hold off
hold_off(wasHold);

%% Outputs
handles.im = im;
handles.ax = ax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%