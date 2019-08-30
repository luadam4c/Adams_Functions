function h = plot_horizontal_shade (varargin)
%% Plots a shaded area at specific y values, either between specific x values or extend to the current x-axis limits
% Usage: h = plot_horizontal_shade (y, xLow, xHigh, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       plot_horizontal_shade([10, 20])
%       plot_horizontal_shade(1:5, rand(5, 1), rand(5, 1) + 2)
%       plot_horizontal_shade([1, 2], 'Color', 'Blue')
%       plot_horizontal_shade({[10, 20], [30, 40]})
%
% Outputs:
%       h           - handle to the shade
%                   specified as a handle to a Patch object
%
% Arguments:
%       y           - (opt) y value(s)
%                   must be a numeric vector
%       xLow        - (opt) low x value(s)
%                   must be a numeric vector
%       xHigh       - (opt) high x value(s)
%                   must be a numeric vector
%       varargin    - 'LineStyle': line style of boundaries
%                   must be an unambiguous, case-insensitive match to one of: 
%                       '-'     - solid line
%                       '--'    - dashed line
%                       ':'     - dotted line
%                       '-.'    - dash-dotted line
%                       'none'  - no line
%                   default == 'none'
%                   - 'ColorMap': color map passed in
%                   must be empty or a string/character vector
%                       or an n-by-3 numeric array
%                   default == LightGray
%                   - Any other parameter-value pair for the fill() function
%
% Requires:
%       cd/plot_vertical_shade.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/plot_window_boundaries.m

% File History:
% 2019-08-27 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
xDefault = [];
yLowDefault = [];
yHighDefault = [];
horizontalInsteadDefault = false;
lineStyleDefault = 'none';
colorMapDefault = 'LightGray';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add optional inputs to the Input Parser
addOptional(iP, 'y', xDefault);
addOptional(iP, 'xLow', yLowDefault);
addOptional(iP, 'xHigh', yHighDefault);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'LineStyle', lineStyleDefault, ...
    @(y) all(islinestyle(y, 'ValidateMode', true)));
addParameter(iP, 'ColorMap', colorMapDefault);

% Read from the Input Parser
parse(iP, varargin{:});
y = iP.Results.y;
xLow = iP.Results.xLow;
xHigh = iP.Results.xHigh;
[~, lineStyle] = islinestyle(iP.Results.LineStyle, 'ValidateMode', true);
colorMap = iP.Results.ColorMap;

% Keep unmatched arguments for the fill() function
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
h = plot_vertical_shade(y, xLow, xHigh, 'HorizontalInstead', true, ...
                        'LineStyle', lineStyle, 'ColorMap', colorMap, ...
                        otherArguments{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%