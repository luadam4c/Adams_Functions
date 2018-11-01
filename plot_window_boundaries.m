function plot_window_boundaries (win, varargin)
%% TODO: A summary of what the function does (must be a single unbreaked line)
% Usage: plot_window_boundaries (win, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Arguments:
%       win         - window to plot boundaries for
%                   must be a 2-element numeric vector
%       varargin    - 'LineColor': color of boundaries
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
%
% Requires:
%       cd/islinestyle.m
%
% Used by:    
%       cd/m3ha_plot_individual_traces.m

% File History:
% 2018-10-29 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
lineColorDefault = 'g';     % boundaries are green by default
lineStyleDefault = '--';    % boundaries are dotted lines by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'win', ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'LineColor', lineColorDefault);
addParameter(iP, 'LineStyle', lineStyleDefault, ...
    @(x) all(islinestyle(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, win, varargin{:});
lineColor = iP.Results.LineColor;
[~, lineStyle] = islinestyle(iP.Results.LineStyle, 'ValidateMode', true);

%% Preparation
% Get the y-axis limits
yLimits = get(gca, 'YLim');

%% Do the job
% Plot lines spanning the y-axis
line(win(1) * ones(1, 2), yLimits, ...
    'Color', lineColor, 'LineStyle', lineStyle);
line(win(2) * ones(1, 2), yLimits, ...
    'Color', lineColor, 'LineStyle', lineStyle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%