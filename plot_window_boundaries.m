function h = plot_window_boundaries (win, varargin)
%% Plots window boundaries as vertical lines
% Usage: h = plot_window_boundaries (win, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       h           - handles to each line object (left, right)
%                   specified as a 2-element column array 
%                       of primitive line object handles
% Arguments:
%       win         - window to plot boundaries for
%                   must be a 2-element numeric vector
%       varargin    - 'YLimits': y value limits for the line
%                   must be empty or a numeric vector of 2 elements
%                   default == get(gca, 'YLim')
%                   - 'LineColor': color of boundaries
%                   must be recognized by the plot() function
%                   default == 'r'
%                   - 'LineStyle': line style of boundaries
%                   must be an unambiguous, case-insensitive match to one of: 
%                       '-'     - solid line
%                       '--'    - dashed line
%                       ':'     - dotted line
%                       '-.'    - dash-dotted line
%                       'none'  - no line
%                   default == '-'
%                   - Any other parameter-value pair for the line() function
% 
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_numeric.m
%       cd/islinestyle.m
%       cd/plot_vertical_line.m
%
% Used by:    
%       cd/m3ha_plot_individual_traces.m

% File History:
% 2018-10-29 Created by Adam Lu
% 2018-12-19 Now passes done extra arguments
% 2018-12-19 Now returns object handles
% 

%% Default values for optional arguments
yLimitsDefault = [];
lineStyleDefault = '--';    % boundaries are dotted lines by default
lineColorDefault = 'g';     % boundaries are green by default

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
addRequired(iP, 'win', ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'LineColor', lineColorDefault);
addParameter(iP, 'LineStyle', lineStyleDefault, ...
    @(x) all(islinestyle(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, win, varargin{:});
yLimits = iP.Results.YLimits;
lineColor = iP.Results.LineColor;
[~, lineStyle] = islinestyle(iP.Results.LineStyle, 'ValidateMode', true);

% Keep unmatched arguments for the line() function
otherArguments = iP.Unmatched;

%% Preparation
% Set default y value limits
if isempty(yLimits)
    yLimits = get(gca, 'YLim');
end

% Force as a column
win = force_column_numeric(win);

%% Do the job
% Plot lines
h = arrayfun(@(x) plot_vertical_line(x, 'YLimits', yLimits, ...
                'Color', lineColor, 'LineStyle', lineStyle, otherArguments), ...
            win);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

line(win(1) * ones(1, 2), yLimits, ...
    'Color', lineColor, 'LineStyle', lineStyle);
line(win(2) * ones(1, 2), yLimits, ...
    'Color', lineColor, 'LineStyle', lineStyle);

% Initialize a graphics object handle array
h = gobjects(2, 1);

h(1) = plot_vertical_line(win(1), 'YLimits', yLimits, ...
                'Color', lineColor, 'LineStyle', lineStyle, otherArguments);
h(2) = plot_vertical_line(win(2), 'YLimits', yLimits, ...
                'Color', lineColor, 'LineStyle', lineStyle, otherArguments);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%