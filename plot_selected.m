function selected = plot_selected (xValues, yValues, indSelected, varargin)
%% Plots selected values
% Usage: selected = plot_selected (xValues, yValues, indSelected, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       xValues = transpose(1:10);
%       yValues = randi(numel(xValues), size(xValues));
%       indSelected = [2, 6];
%       plot_tuning_curve(xValues, yValues)
%       hold on
%       plot_selected(xValues, yValues, indSelected)
%
% Outputs:
%       selected  	- TODO: Description of handles
%                   specified as a TODO
%
% Arguments:
%       xValues     - x values
%                   must be a TODO
%       yValues     - y values
%                   must be a TODO
%       indSelected - indices selected
%                   must be a TODO
%       varargin    - 'Marker': type of markers
%                   default == 'o'
%                   - 'LineStyle': line style of connecting the markers
%                   must be an unambiguous, case-insensitive match to one of: 
%                       '-'     - solid line
%                       '--'    - dashed line
%                       ':'     - dotted line
%                       '-.'    - dash-dotted line
%                       'none'  - no line
%                   default == 'none'
%                   - 'LineWidth': line width of markers
%                   must be empty or a positive scalar
%                   default == 3
%                   - 'ColorMap': color map passed in
%                   must be empty or a string/character vector
%                       or an n-by-3 numeric array
%                   default == 'r'
%                   - 'AxesHandle': axes handle to plot on
%                   must be a empty or a axes object handle
%                   default == gca
%                   - Any other parameter-value pair for plot()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/decide_on_colormap.m
%       cd/force_column_cell.m
%       cd/hold_on.m
%       cd/hold_off.m
%       cd/match_row_count.m
%
% Used by:
%       cd/m3ha_plot_simulated_traces.m
%       cd/plot_bar.m
%       cd/plot_tuning_curve.m

% File History:
% 2019-11-24 Created by Adam Lu
% TODO: Allow inputs to be cell arrays

%% Hard-coded parameters

%% Default values for optional arguments
markerDefault = 'x';            % plot crosses by default
lineStyleDefault = 'none';      % no line connecting markers by default
lineWidthDefault = 3;           % line width of 3 points by default
colorMapDefault = 'r';          % red markers by default
axHandleDefault = [];           % gca by default

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
%TODO: Add requirements
addRequired(iP, 'xValues');
addRequired(iP, 'yValues');
addRequired(iP, 'indSelected');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Marker', markerDefault);
addParameter(iP, 'LineStyle', lineStyleDefault, ...
    @(x) all(islinestyle(x, 'ValidateMode', true)));
addParameter(iP, 'LineWidth', lineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'ColorMap', colorMapDefault);
addParameter(iP, 'AxesHandle', axHandleDefault);

% Read from the Input Parser
parse(iP, xValues, yValues, indSelected, varargin{:});
marker = iP.Results.Marker;
[~, lineStyle] = islinestyle(iP.Results.LineStyle, 'ValidateMode', true);
lineWidth = iP.Results.LineWidth;
colorMap = iP.Results.ColorMap;
axHandle = iP.Results.AxesHandle;

% Keep unmatched arguments for the plot() function
otherArguments = iP.Unmatched;

%% Preparation
% Force as column cell arrays
[xValues, yValues, indSelected] = ...
    argfun(@force_column_cell, xValues, yValues, indSelected);

% Match row count 
nRows = max([numel(xValues), numel(yValues), numel(indSelected)]);

% Match row count 
[xValues, yValues, indSelected] = ...
    argfun(@(x) match_row_count(x, nRows), xValues, yValues, indSelected);

% Convert color map to a cell array or row vectors
colorMapCell = decide_on_colormap(colorMap, nRows, 'ForceCellOutput', true);

%% Do the job
% Hold on
wasHold = hold_on(axHandle);

% Plot selected
selected = cellfun(@(a, b, c, d) plot_selected_helper(a, b, c, d, ...
                           	marker, lineStyle, lineWidth, axHandle, otherArguments), ...
                    xValues, yValues, indSelected, colorMapCell, ...
                    'UniformOutput', false);

% Hold off
hold_off(wasHold, axHandle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function selected = plot_selected_helper(xValues, yValues, ...
                                            indSelected, colorMap, ...
                                            marker, lineStyle, lineWidth, ...
                                            axHandle, otherArguments)
% Plots a set of selected values

%% Preparation
% Remove NaN values from indices
indSelected = indSelected(~isnan(indSelected));

% If no indices remaining, return
if isempty(indSelected) || isempty(xValues) || isempty(yValues)
    selected = gobjects;
    return
end

%% Do the job
% Selected x locations
xLocsSelected = xValues(indSelected);

% Selected y locations
yLocsSelected = yValues(indSelected, :);

% Plot values
selected = plot(axHandle, xLocsSelected, yLocsSelected, ...
                'LineStyle', lineStyle, 'Marker', marker, ...
                'Color', colorMap, 'LineWidth', lineWidth, otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%