function [h, xValues, yValues] = plot_ellipse (center, halflengths, theta0, varargin)
%% Plot an ellipse that may be oblique
% Usage: [h, xValues, yValues] = plot_ellipse (center, halflengths, theta0, varargin)
% Explanation:
%       Plots an ellipse with a given center, half-axis lengths 
%           and rotation angle (radians). Example:
%               plot_ellipse([2, 3], [3, 2], pi/6);
% Outputs:
%       h           - the ellipse
%                   specified as a chart line object
% Side Effects:
%       Plots an ellipse
% Arguments:    
%       center      - center of ellipse
%                   must be a 2-element numeric vector
%       halflengths - halflengths of ellipse
%                   must be a 2-element numeric, positive vector
%       theta0      - angle between x axis and 
%                       the first axis of ellipse (radians)
%                   must be a numeric scalar
%       varargin    - 'Color': color of ellipse
%                   must be recognized by the plot() function
%                   default == 'r'
%                   - 'NPoints': number of points to plot
%                   must be a positive integer scalar
%                   default == 1000
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
%                   - 'ToPlot': whether to plot
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%
% Used by:    
%       /home/Matlab/Adams_Functions/plot_grouped_scatter.m
%
% File History:
% 2017-12-15 Created by Adam Lu
% 

%% Default values for optional arguments
nPointsDefault = 1000;                  % default number of points to plot
colorDefault = 'r';                     % default color of ellipse
lineStyleDefault = '-';                 % default line style of ellipse
lineWidthDefault = 1;                   % default line width of ellipse
toPlotDefault = true;                   % whether to plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(['Not enough input arguments, ', ...
            'type ''help plot_ellipse'' for usage']);
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
addParameter(iP, 'Color', colorDefault);
addParameter(iP, 'NPoints', nPointsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'LineStyle', lineStyleDefault, ...
    @(x) any(validatestring(x, {'-', '--', ':', '-.', 'none'})));
addParameter(iP, 'LineWidth', lineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'ToPlot', toPlotDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, center, halflengths, theta0, varargin{:});
nPoints = iP.Results.NPoints;
color = iP.Results.Color;
lineStyle = iP.Results.LineStyle;
lineWidth = iP.Results.LineWidth;
toPlot = iP.Results.ToPlot;

% Display warning message if some inputs are unmatched
if ~isempty(fieldnames(iP.Unmatched))
    fprintf('WARNING: The following name-value pairs could not be parsed: \n');
    disp(iP.Unmatched);
end

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
            'LineStyle', lineStyle, 'Color', color, ...
            'LineWidth', lineWidth);
else h = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
