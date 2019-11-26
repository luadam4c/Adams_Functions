function handles = plot_selected (xValues, yValues, indSelected, varargin)
%% Plots selected values
% Usage: handles = plot_selected (xValues, yValues, indSelected, varargin)
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
%       handles     - TODO: Description of handles
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
%                   - Any other parameter-value pair for plot()
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/plot_bar.m
%       cd/plot_tuning_curve.m

% File History:
% 2019-11-24 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
markerDefault = 'x';            % plot crosses by default
lineStyleDefault = 'none';      % no line connecting markers by default
lineWidthDefault = 3;           % line width of 3 points by default
colorMapDefault = 'r';          % red markers by default

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

% Read from the Input Parser
parse(iP, xValues, yValues, indSelected, varargin{:});
marker = iP.Results.Marker;
[~, lineStyle] = islinestyle(iP.Results.LineStyle, 'ValidateMode', true);
lineWidth = iP.Results.LineWidth;
colorMap = iP.Results.ColorMap;

% Keep unmatched arguments for the plot() function
otherArguments = iP.Unmatched;

%% Preparation
% Remove NaN values from indices
indSelected = indSelected(~isnan(indSelected));

% If no indices remaining, return
if isempty(indSelected)
    handles = gobjects;
    return
end

%% Do the job
% Selected x locations
xLocsSelected = xValues(indSelected);

% Selected y locations
yLocsSelected = yValues(indSelected, :);

% Plot values
selected = plot(xLocsSelected, yLocsSelected, ...
                'LineStyle', lineStyle, 'Marker', marker, ...
                'Color', colorMap, 'LineWidth', lineWidth, otherArguments);

%% Output results
handles = selected;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%