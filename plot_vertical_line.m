function h = plot_vertical_line (xValue, varargin)
%% Plots vertical line(s)
% Usage: h = plot_vertical_line (xValue, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       h = plot_vertical_line(xValue)
%       h = plot_vertical_line(xValue, 'YLimits', yLimits)
%       h = plot_vertical_line(3)
%       h = plot_vertical_line(3, 'YLimits', [])
%       h = plot_vertical_line(3, 'YLimits', [0, 0])
%       h = plot_vertical_line(3, 'YLimits', [1, 2])
%       h = plot_vertical_line(3, 'YLimits', [1, 2, 4, 5])
%       h = plot_vertical_line([3, 4, 5])
%       h = plot_vertical_line([3 4], 'YLimits', {[2 4], [1 2 4 5]})
%       h = plot_vertical_line([3 4], 'YLimits', {[2 4], [1 2 4 5]})
%       h = plot_vertical_line([3, 4, 5], 'Color', 'r')
%       h = plot_vertical_line(3, 'HorizontalInstead', true)
%
% Outputs:
%       h           - handle to the line object(s) created
%                   specified as a primitive line object handle
% Arguments:
%       xValue      - the x value(s) for the vertical line(s)
%                   must be a numeric, datetime or duration array
%       varargin    - 'HorizontalInstead': whether to plot a horizontal shade
%                                               instead
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'YLimits': y value limits for the line(s)
%                   must be empty or a numeric vector of 2 elements
%                       or an array of 2 rows
%                   default == get(gca, 'YLim')
%                   - 'ColorMap' - color map used
%                   must be a 2-D numeric array with 3 columns
%                   default == decide_on_colormap([], nLines)
%                   - Any other parameter-value pair for the line() function
%
% Requires:
%       cd/apply_over_cells.m
%       cd/create_error_for_nargin.m
%       cd/decide_on_colormap.m
%       cd/force_column_cell.m
%       cd/match_format_vector_sets.m
%
% Used by:
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/plot_bar.m
%       cd/plot_error_bar.m
%       cd/plot_horizontal_line.m
%       cd/plot_pulse_response_with_stimulus.m
%       cd/plot_tuning_curve.m
%       cd/plot_swd_histogram.m
%       cd/plot_window_boundaries.m

% File History:
% 2018-12-19 Created by Adam Lu
% 2018-12-27 Now allows xValue to be an array
% 2018-12-27 Now accepts datetime and duration arrays
% 2019-01-24 Now accepts multiple y limits
% 2019-08-30 Added 'ColorMap' as an optional argument
% TODO: Finish up 'HorizontalInstead' and use this in plot_horizontal_line
% TODO: Allow 2-D arrays for x values
% 

%% Hard-coded parameters

%% Default values for optional arguments
yLimitsDefault = [];
colorMapDefault = [];               % set later
horizontalInsteadDefault = false;

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
addRequired(iP, 'xValue', ...
    @(x) validateattributes(x, {'numeric', 'datetime', 'duration'}, {'3d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'YLimits', yLimitsDefault);
addParameter(iP, 'ColorMap', colorMapDefault);
addParameter(iP, 'HorizontalInstead', horizontalInsteadDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, xValue, varargin{:});
yLimits = iP.Results.YLimits;
colorMap = iP.Results.ColorMap;
horizontalInstead = iP.Results.HorizontalInstead;

% Keep unmatched arguments for the line() function
otherArguments = iP.Unmatched;

%% Preparation
% Set default y value limits
if isempty(yLimits)
    yLimits = get(gca, 'YLim');
end

% Force as a cell array of column vectors and match vectors
[xValue, yLimits] = match_format_vector_sets(num2cell(xValue), yLimits);

% Place in column cell arrays and 
%   expand x limits if there are more than 2 values
[yLimitsCell, xValueCell] = cellfun(@(x, y) expand_limits(x, y), yLimits, xValue, ...
                                    'UniformOutput', false);

% Vertically concatenate all column cell arrays
xValueAll = apply_over_cells(@vertcat, xValueCell);
yLimitsAll = apply_over_cells(@vertcat, yLimitsCell);

% Count the number of y values
nYValues = numel(xValue);

% Count the number of lines for each y value
nLinesEachX = count_vectors(xValueCell);

% Count the number of lines to plot
nLines = numel(xValueAll);

% Decide on the number of colors to plot
nColors = nYValues;

% Set default color map
colorMap = decide_on_colormap(colorMap, nColors);

% Expand to nLines
% TODO: Add an option to decide_on_colormap.m to do this 'ExpandBy'
colorMapCell = arrayfun(@(x) repmat(colorMap(x, :), nLinesEachX(x), 1), ...
                    1:nYValues, 'UniformOutput', false);
colorMapExpanded = vertcat(colorMapCell{:});

%% Do the job
% Hold on if plotting more than one line
if nLines > 1
    hold on;
end

% Plot all lines
if horizontalInstead
    h = cellfun(@(y, x, z) line(x, repmat(y, size(x)), ...
                                'Color', colorMapExpanded(z, :), ...
                                otherArguments), ...
                xValueAll, yLimitsAll, num2cell(transpose(1:nLines)));
else
    h = cellfun(@(x, y, z) line(repmat(x, size(y)), y, ...
                                'Color', colorMapExpanded(z, :), ...
                                otherArguments), ...
                xValueAll, yLimitsAll, num2cell(transpose(1:nLines)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [yLimitsCell, xValueCell] = expand_limits(yLimits, xValue)

nXEndpoints = numel(yLimits);

if mod(nXEndpoints, 2) ~= 0
    error('Number of x endpoints must be even!');
end

% Actual number of lines to plot
nLines = nXEndpoints / 2;

% Reshape as two rows
yLimits = reshape(yLimits, 2, nLines);

% Force as column cell array of column vectors
yLimitsCell = force_column_cell(yLimits);

% Expand y value accordingly
xValueCell = repmat({xValue}, nLines, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%