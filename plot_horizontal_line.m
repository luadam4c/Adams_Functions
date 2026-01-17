function h = plot_horizontal_line (yValue, varargin)
%% Plots horizontal line(s)
% Usage: h = plot_horizontal_line (yValue, varargin)
% Explanation:
%       TODO
%       Note: This function plots horizontal lines with the same y value
%               as the same color
%
% Example(s):
%       h = plot_horizontal_line(yValue)
%       h = plot_horizontal_line(yValue, 'XLimits', xLimits)
%       h = plot_horizontal_line(3)
%       h = plot_horizontal_line(3, 'XLimits', [])
%       h = plot_horizontal_line(3, 'XLimits', [0, 0])
%       h = plot_horizontal_line(3, 'XLimits', [1, 2])
%       h = plot_horizontal_line(3, 'XLimits', [1, 2, 4, 5])
%       h = plot_horizontal_line([3, 4, 5])
%       h = plot_horizontal_line([3 4], 'XLimits', {[2 4], [1 2 4 5]})
%       h = plot_horizontal_line([3 4], 'XLimits', {[2 4], [1 2 4 5]})
%       h = plot_horizontal_line([3, 4, 5], 'Color', 'r')
%       TODO: h = plot_horizontal_line([3, 4, 5], 'ColorMapFunc', 'jet')
%
% Outputs:
%       h           - handle to the line object(s) created
%                   specified as a primitive line object handle array
% Arguments:
%       yValue      - the y value(s) for the horizontal line(s)
%                   must be a numeric, datetime or duration array
%       varargin    - 'XLimits': x value limits for the line(s)
%                   must be empty or a numeric vector of 2 elements
%                       or an array of 2 rows
%                   default == get(gca, 'XLim')
%                   - 'ColorMap' - color map used
%                   must be a 2-D numeric array with 3 columns
%                   default == decide_on_colormap([], nLines)
%                   - Any other parameter-value pair for the plot_vertical_line() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/plot_vertical_line.m
%
% Used by:
%       cd/m3ha_simulate_population.m
%       cd/parse_multiunit.m
%       cd/plot_autocorrelogram.m
%       cd/plot_bar.m
%       cd/plot_error_bar.m
%       cd/plot_raw_multiunit.m
%       cd/plot_pulse_response_with_stimulus.m
%       cd/plot_spike_histogram.m
%       cd/plot_struct.m
%       cd/plot_tuning_curve.m
%       cd/plot_window_boundaries.m

% File History:
% 2018-12-19 Created by Adam Lu
% 2018-12-27 Now allows yValue to be an array
% 2018-12-27 Now accepts datetime and duration arrays
% 2019-01-24 Now accepts multiple x limits
% 2019-03-17 Allow each x limits to be of length > 2 and break them up
%               into pairs
% 2019-08-22 Added 'ColorMap' as an optional argument
% 2019-08-23 Made sure that each y value is a single color
% 2019-08-30 Now uses plot_vertical_line.m with option 'HorizontalInstead'

%% Hard-coded parameters

%% Default values for optional arguments
xLimitsDefault = [];
colorMapDefault = [];               % set later

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
addRequired(iP, 'yValue', ...
    @(x) validateattributes(x, {'numeric', 'datetime', 'duration'}, {'3d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'XLimits', xLimitsDefault);
addParameter(iP, 'ColorMap', colorMapDefault);

% Read from the Input Parser
parse(iP, yValue, varargin{:});
xLimits = iP.Results.XLimits;
colorMap = iP.Results.ColorMap;

% Keep unmatched arguments for the plot_vertical_line() function
otherArguments = iP.Unmatched;

%% Do the job
% Use plot_vertical_line
%   Note: 'ColorMap' still needs to be parsed and passed on
%           so that 'Color' will be accepted as 'ColorMap'
h = plot_vertical_line(yValue, 'YLimits', xLimits, 'ColorMap', colorMap, ...
                    'HorizontalInstead', true, otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
