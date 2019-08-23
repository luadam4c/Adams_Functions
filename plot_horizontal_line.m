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
%                   - Any other parameter-value pair for the line() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/decide_on_colormap.m
%       cd/force_column_cell.m
%       cd/isnum.m
%       cd/match_format_vector_sets.m
%       cd/match_row_count.m
%
% Used by:
%       cd/parse_multiunit.m
%       cd/plot_bar.m
%       cd/plot_error_bar.m
%       cd/plot_pulse_response_with_stimulus.m
%       cd/plot_struct.m
%       cd/plot_tuning_curve.m

% File History:
% 2018-12-19 Created by Adam Lu
% 2018-12-27 Now allows yValue to be an array
% 2018-12-27 Now accepts datetime and duration arrays
% 2019-01-24 Now accepts multiple x limits
% 2019-03-17 Allow each x limits to be of length > 2 and break them up
%               into pairs
% 2019-08-22 Added 'ColorMap' as an optional argument
% 2019-08-23 Made sure that each y value is a single color
% TODO: Allow 2-D arrays for y values

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

% Keep unmatched arguments for the line() function
otherArguments = iP.Unmatched;

%% Preparation
% Set default x value limits
if isempty(xLimits)
    xLimits = get(gca, 'XLim');
end

% Force as a cell array of column vectors and match vectors
[yValue, xLimits] = match_format_vector_sets(num2cell(yValue), xLimits);

% Place in column cell arrays and 
%   expand x limits if there are more than 2 values
[xLimitsCell, yValueCell] = cellfun(@(x, y) expand_limits(x, y), xLimits, yValue, ...
                                    'UniformOutput', false);

% Vertically concatenate all column cell arrays
yValueAll = apply_over_cells(@vertcat, yValueCell);
xLimitsAll = apply_over_cells(@vertcat, xLimitsCell);

% Count the number of y values
nYValues = numel(yValue);

% Count the number of lines for each y value
nLinesEachY = count_vectors(yValueCell);

% Count the number of lines to plot
nLines = numel(yValueAll);

% Decide on the number of colors to plot
nColors = nYValues;

% Set default color map
colorMap = decide_on_colormap(colorMap, nColors);

% Expand to nLines
% TODO: Add an option to decide_on_colormap.m to do this 'ExpandBy'
colorMapCell = arrayfun(@(x) repmat(colorMap(x, :), nLinesEachY(x), 1), ...
                    1:nYValues, 'UniformOutput', false);
colorMapExpanded = vertcat(colorMapCell{:});

%% Do the job
% Hold on if plotting more than one line
if nLines > 1
    hold on;
end

% Plot all lines
h = cellfun(@(y, x, z) line(x, repmat(y, size(x)), ...
                        'Color', colorMapExpanded(z, :), otherArguments), ...
            yValueAll, xLimitsAll, num2cell(transpose(1:nLines)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLimitsCell, yValueCell] = expand_limits(xLimits, yValue)

nXEndpoints = numel(xLimits);

if mod(nXEndpoints, 2) ~= 0
    error('Number of x endpoints must be even!');
end

% Actual number of lines to plot
nLines = nXEndpoints / 2;

% Reshape as two rows
xLimits = reshape(xLimits, 2, nLines);

% Force as column cell array of column vectors
xLimitsCell = force_column_cell(xLimits);

% Expand y value accordingly
yValueCell = repmat({yValue}, nLines, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

@(x) isempty(x) || isnumeric(x) && isvector(x) && length(x) == 2);

h = arrayfun(@(y) line(xLimits, repmat(y, size(xLimits)), otherArguments), ...
            yValue);

if isfield(otherArguments, 'Color')
    h = cellfun(@(y, x) line(x, repmat(y, size(x)), otherArguments), ...
                yValueAll, xLimitsAll);
else
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
