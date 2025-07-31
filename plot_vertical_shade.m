function h = plot_vertical_shade (varargin)
%% Plots a shaded area at specific x values, either between specific y values or extend to the current y-axis limits
% Usage: h = plot_vertical_shade (x, yLow, yHigh, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       plot_vertical_shade([10, 20])
%       plot_vertical_shade(1:5, rand(5, 1), rand(5, 1) + 2)
%       plot_vertical_shade([1, 2], 'Color', 'Blue')
%       plot_vertical_shade({[10, 20], [30, 40]})
%       plot_vertical_shade([1, 2], 'HorizontalInstead', true)
%
% Outputs:
%       h           - handle to the shade
%                   specified as a handle to a Patch object
%
% Arguments:
%       x           - (opt) x value(s)
%                   must be a numeric vector
%       yLow        - (opt) low y value(s)
%                   must be a numeric vector
%       yHigh       - (opt) high y value(s)
%                   must be a numeric vector
%       varargin    - 'HorizontalInstead': whether to plot a horizontal shade
%                                               instead
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'LineStyle': line style of boundaries
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
%       cd/argfun.m
%       cd/decide_on_colormap.m
%       cd/force_column_cell.m
%       cd/islinestyle.m
%       cd/match_format_vector_sets.m
%       cd/set_default_flag.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/plot_horizontal_shade.m
%       cd/plot_psth.m
%       cd/plot_relative_events.m
%       cd/plot_window_boundaries.m
%       cd/virt_golomb_generate_output.m

% File History:
% 2019-08-27 Created by Adam Lu
% 2019-08-30 Now accepts cell arrays as inputs
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
addOptional(iP, 'x', xDefault);
addOptional(iP, 'yLow', yLowDefault);
addOptional(iP, 'yHigh', yHighDefault);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'HorizontalInstead', horizontalInsteadDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'LineStyle', lineStyleDefault, ...
    @(x) all(islinestyle(x, 'ValidateMode', true)));
addParameter(iP, 'ColorMap', colorMapDefault);

% Read from the Input Parser
parse(iP, varargin{:});
x = iP.Results.x;
yLow = iP.Results.yLow;
yHigh = iP.Results.yHigh;
horizontalInstead = iP.Results.HorizontalInstead;
[~, lineStyle] = islinestyle(iP.Results.LineStyle, 'ValidateMode', true);
colorMap = iP.Results.ColorMap;

% Keep unmatched arguments for the fill() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Initialize handle
h = gobjects(1);

% Match the vectors if any input is a cell array
if iscell(x) || iscell(yLow) || iscell(yHigh)
    % Force as column cell arrays
    x = force_column_cell(x);

    % Match formats of vectors
    [x, yLow] = match_format_vector_sets(x, yLow);
    [x, yHigh] = match_format_vector_sets(x, yHigh);

    % Plot many vertical shades
    h = cellfun(@(x, y, z) plot_vertical_shade_helper(x, y, z, ...
                            horizontalInstead, lineStyle, colorMap, ...
                            otherArguments), ...
                x, yLow, yHigh);
else
    % Plot one vertical shade
    h = plot_vertical_shade_helper(x, yLow, yHigh, ...
        horizontalInstead, lineStyle, colorMap, otherArguments);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = plot_vertical_shade_helper(x, yLow, yHigh, ...
                        horizontalInstead, lineStyle, colorMap, otherArguments)
%% Plots one vertical shade

% Force as column vectors
[x, yLow, yHigh] = argfun(@force_column_vector, x, yLow, yHigh);

% Compute number of each type of inputs
[nXs, nYLows, nYHighs] = argfun(@numel, x, yLow, yHigh);

% Make sure yLow and yHigh are of equal length
if nYLows > 1 && nYHighs > 1 && nYLows ~= nYHighs
    disp('yLow and yHigh must be of equal length if any of them are not 1!');
    return
end

% Compute number of x values to plot
nXToPlot = max([nXs, nYLows, nYHighs]);

% Set default x values
if isempty(x)
    if nXToPlot > 1
        x = transpose(1:nXToPlot);
    else
        if horizontalInstead
            x = get(gca, 'XLim');
        else
            x = get(gca, 'YLim');
        end
    end
end

% Get current y limits
if isempty(yLow) || isempty(yHigh)
    if horizontalInstead
        yLimits = get(gca, 'XLim');
    else
        yLimits = get(gca, 'YLim');
    end
end

% Set default y lows
if isempty(yLow)
    yLow = repmat(yLimits(1), nXToPlot, 1);
end

% Set default y highs
if isempty(yHigh)
    yHigh = repmat(yLimits(2), nXToPlot, 1);
end

% Decide on the color map
colorMap = decide_on_colormap(colorMap, 1);

%% Do the job
% Save whether was hold
wasHold = set_default_flag([], ishold);

% Hold on
hold on

% The x and y values for the confidence intervals
xValues = [x; flipud(x)];
yValues = [yHigh; flipud(yLow)];

% Fill the area between xValues and yValues
if horizontalInstead
    h = fill_helper(yValues, xValues, colorMap, lineStyle, otherArguments);
else
    h = fill_helper(xValues, yValues, colorMap, lineStyle, otherArguments);
end

% Hold off
if ~wasHold
    hold off
end

% Move it to the bottom of the figure stack
uistack(h, 'bottom');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = fill_helper (x, y, colorMap, lineStyle, otherArguments)

h = fill(x, y, colorMap, 'LineStyle', lineStyle, otherArguments{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%