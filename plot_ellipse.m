function [h, xValues, yValues] = plot_ellipse (center, halflengths, theta0, varargin)
%% Plot an ellipse that may be oblique
% Usage: [h, xValues, yValues] = plot_ellipse (center, halflengths, theta0, varargin)
% Explanation:
%       Plots an ellipse with a given center, half-axis lengths 
%           and rotation angle (radians). 
%
% Example:
%               plot_ellipse([2, 3], [3, 2], pi/6);
%
% Outputs:
%       h           - the ellipse
%                   specified as a chart line object
%
% Side Effects:
%       Plots an ellipse
%
% Arguments:    
%       center      - center of ellipse
%                   must be a 2-element numeric vector
%       halflengths - halflengths of ellipse
%                   must be a 2-element numeric, positive vector
%       theta0      - angle between x axis and 
%                       the first axis of ellipse (radians)
%                   must be a numeric scalar
%       varargin    - 'ToPlot': whether to plot
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'NPoints': number of points to plot
%                   must be a positive integer scalar
%                   default == 1000
%                   - Any other parameter-value pair for the plot() function
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/plot_grouped_scatter.m
%       cd/ZG_plot_grouped_scatter.m

% File History:
% 2017-12-15 Created by Adam Lu
% 2018-05-16 Now uses islinestyle.m
% 2018-12-18 Now uses iP.Unmatched
% 

%% Default values for optional arguments
toPlotDefault = true;                   % whether to plot
nPointsDefault = 1000;                  % default number of points to plot
lineColorDefault = 'r';                 % default line color of ellipse
lineStyleDefault = '-';                 % default line style of ellipse
lineWidthDefault = 1;                   % default line width of ellipse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'plot_ellipse';
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'center', ...                   % center of ellipse
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
addRequired(iP, 'halflengths', ...              % halflengths of ellipse
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive', 'numel', 2}));
addRequired(iP, 'theta0', ...                   % angle of ellipse
    @(x) validateattributes(x, {'numeric'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ToPlot', toPlotDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'NPoints', nPointsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Read from the Input Parser
parse(iP, center, halflengths, theta0, varargin{:});
nPoints = iP.Results.NPoints;
toPlot = iP.Results.ToPlot;

% Keep unmatched arguments for the plot() function
otherArguments = iP.Unmatched;

%% Prepare ellipse
% Parametric variable
t = linspace(0, 2*pi, nPoints);     % this is a row vector

% Extract parameters from arguments
a = halflengths(1);
b = halflengths(2);

% Obtain unrotated and unshifted x and y values using the parametric equation
canonical = [a * cos(t); ...
             b * sin(t)];

% Prepare rotation matrix
R = [cos(theta0), -sin(theta0); ...
     sin(theta0), cos(theta0)];

% Obtain rotated x and y values
rotated = R * canonical;

% Obtain shifted x and y values
center = center(:);         % make sure it's a column vector
shifted = repmat(center, 1, nPoints) + rotated;

%% Plot ellipse
xValues = shifted(1, :);
yValues = shifted(2, :);
if toPlot
    h = plot(xValues, yValues, ...
            'LineStyle', lineStyleDefault, 'Color', lineColorDefault, ...
            'LineWidth', lineWidthDefault, otherArguments);
else 
    h = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%% Hard-coded parameters
validLineStyles = {'-', '--', ':', '-.', 'none'};

addParameter(iP, 'LineStyle', lineStyleDefault, ...
    @(x) any(validatestring(x, validLineStyles)));

% Display warning message if some inputs are unmatched
if ~isempty(fieldnames(iP.Unmatched))
    fprintf('WARNING: The following name-value pairs could not be parsed: \n');
    disp(iP.Unmatched);
end

% Display warning message if some inputs are unmatched
if ~isempty(fieldnames(iP.Unmatched))
    fprintf('WARNING: The following name-value pairs could not be parsed: \n');
    disp(iP.Unmatched);
end

addParameter(iP, 'LineColor', lineColorDefault);
addParameter(iP, 'LineStyle', lineStyleDefault, ...
    @(x) all(islinestyle(x, 'ValidateMode', true)));
addParameter(iP, 'LineWidth', lineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
lineColor = iP.Results.LineColor;
[~, lineStyle] = islinestyle(iP.Results.LineStyle, 'ValidateMode', true);
lineWidth = iP.Results.LineWidth;
%                   - 'LineColor': color of ellipse
%                   must be recognized by the plot() function
%                   default == 'r'
%                   - 'LineStyle': line style of ellipse
%                   must be an unambiguous, case-insensitive match to one of: 
%                       '-'     - solid line
%                       '--'    - dashed line
%                       ':'     - dotted line
%                       '-.'    - dash-dotted line
%                       'none'  - no line
%                   default == '-'
%                   - 'LineWidth': line width of ellipse
%                   must be a positive scalar
%                   default == 1
%       cd/islinestyle.m

    h = plot(xValues, yValues, ...
            'LineStyle', lineStyle, 'Color', lineColor, ...
            'LineWidth', lineWidth);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
