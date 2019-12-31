function handles = plot_scale_bar (varargin)
%% Plots a scale bar for x and/or y axis based on the current axis limits
% Usage: handles = plot_scale_bar (axisType (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       figure; plot(1:10, randi(10, 1, 10));
%       plot_scale_bar;
%       plot_scale_bar('x', 'BarLength', 1, 'BarUnits', 'cm');
%       plot_scale_bar('xy', 'XBarUnits', 'ms', 'YBarUnits', 'mV');
%
% Outputs:
%       handles     - TODO: Description of handles
%                   specified as a TODO
%
% Arguments:
%       axisType    - (opt) axis type
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'x'     - x axis
%                       'y'     - y axis
%                       'xy'    - both x and y axis
%                   default == 'xy'
%       varargin    - 'BarLength': bar length in axis units
%                   must be a numeric scalar
%                   default == TODO
%                   - 'XBarLength': x bar length in axis units
%                   must be a numeric scalar
%                   default == TODO
%                   - 'YBarLength': y bar length in axis units
%                   must be a numeric scalar
%                   default == TODO
%                   - 'BarUnits': scale bar units
%                   must be a string scalar or a character vector
%                   default == 'units'
%                   - 'XBarUnits': x scale bar units
%                   must be a string scalar or a character vector
%                   default == same as barUnits
%                   - 'YBarUnits': y scale bar units
%                   must be a string scalar or a character vector
%                   default == same as barUnits
%                   - 'XPosNormalized': x position in normalized units
%                   must be a numeric scalar
%                   default == 0.1
%                   - 'YPosNormalized': y position in normalized units
%                   must be a numeric scalar
%                   default == 0.8
%                   - 'Bar2TextOffsetNormalized': bar to text offset
%                                                   in normalized units
%                   must be a numeric scalar
%                   default == 0.05
%                   - Any other parameter-value pair 
%                       for plot_window_boundaries()
%
% Requires:
%       cd/plot_window_boundaries.m
%
% Used by:
%       cd/m3ha_oscillations_analyze.m
%       cd/m3ha_plot_figure02.m
%       cd/m3ha_plot_figure03.m
%       cd/m3ha_plot_figure05.m

% File History:
% 2019-12-01 Created by Adam Lu
% 

%% Hard-coded parameters
validAxisTypes = {'x', 'y', 'xy'};

% TODO: Make optional arguments
barColor = 'k';
barLineWidth = 1;
relativeBarLength = 0.2;
xPosNormalized = 0.1;
yPosNormalized = 0.8;
bar2TextOffsetNormalized = 0.05;

%% Default values for optional arguments
axisTypeDefault = 'xy';
barLengthDefault = [];              % set later
xBarLengthDefault = [];             % set later
yBarLengthDefault = [];             % set later
barUnitsDefault = 'units';
xBarUnitsDefault = '';              % set later
yBarUnitsDefault = '';              % set later
xPosNormalizedDefault = 0.1;
yPosNormalizedDefault = 0.8;
bar2TextOffsetNormalizedDefault = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add optional inputs to the Input Parser
addOptional(iP, 'axisType', axisTypeDefault, ...
    @(x) any(validatestring(x, validAxisTypes)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BarLength', barLengthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'XBarLength', xBarLengthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'YBarLength', yBarLengthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'BarUnits', barUnitsDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'XBarUnits', xBarUnitsDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'YBarUnits', yBarUnitsDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'XPosNormalized', xPosNormalizedDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'YPosNormalized', yPosNormalizedDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'Bar2TextOffsetNormalized', bar2TextOffsetNormalizedDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));

% Read from the Input Parser
parse(iP, varargin{:});
axisType = iP.Results.axisType;
barLength = iP.Results.BarLength;
xBarLength = iP.Results.XBarLength;
yBarLength = iP.Results.YBarLength;
barUnits = iP.Results.BarUnits;
xBarUnits = iP.Results.XBarUnits;
yBarUnits = iP.Results.YBarUnits;
xPosNormalized = iP.Results.XPosNormalized;
yPosNormalized = iP.Results.YPosNormalized;
bar2TextOffsetNormalized = iP.Results.Bar2TextOffsetNormalized;

% Keep unmatched arguments for the plot_window_boundaries() function
otherArguments = iP.Unmatched;

% Validate axisType
axisType = validatestring(axisType, validAxisTypes);

%% Preparation
% Decide on the scale bars to plot
switch axisType
    case 'x'
        plotXScaleBar = true;
        plotYScaleBar = false;
    case 'y'
        plotXScaleBar = false;
        plotYScaleBar = true;
    case 'xy'
        plotXScaleBar = true;
        plotYScaleBar = true;
    otherwise
        error('axisType unrecognized!');
end

% Decide on the x an y bar units
if isempty(xBarUnits)
    xBarUnits = barUnits;
end
if isempty(yBarUnits)
    yBarUnits = barUnits;
end

%% Do the job
% Get current x and y limits
xLimits = get(gca, 'XLim');
yLimits = get(gca, 'YLim');

% Plot x scale bar
if plotXScaleBar
    % Decide on the appropriate bar length
    if isempty(xBarLength)
        if ~isempty(barLength)
            xBarLength = barLength;
        else
            xBarLength = compute_best_bar_length(xLimits, relativeBarLength);
        end
    end

    % Create a label for the x scale bar
    xBarLabel = sprintf('%g %s', xBarLength, xBarUnits);

    % Decide on the x and y position for the x scale bar
    xPosXBar = xLimits(1) + xPosNormalized * diff(xLimits);
    yPosXBar = yLimits(1) + yPosNormalized * diff(yLimits);

    % Decide on the x bar limits
    xLimitsXBar = xPosXBar + [0, xBarLength];

    % Decide on the x and y position for the x bar text
    xPosXText = xPosXBar;
    yPosXText = yPosXBar - bar2TextOffsetNormalized * diff(yLimits);

    % Plot the x scale bar
    xBar = plot_window_boundaries(xLimitsXBar, 'BarValue', yPosXBar, ...
                        'BoundaryType', 'horizontalBars', ...
                        'Color', barColor, 'LineWidth', barLineWidth, ...
                        otherArguments);

    % Plot the x bar text
    xBarText = text(xPosXText, yPosXText, xBarLabel);
else
    xBarLength = 0;
end

% Plot y scale bar
if plotYScaleBar
    % Decide on the appropriate bar length
    if isempty(yBarLength)
        if ~isempty(barLength)
            yBarLength = barLength;
        else
            yBarLength = compute_best_bar_length(yLimits, relativeBarLength);
        end
    end

    % Create a label for the y scale bar
    yBarLabel = sprintf('%g %s', yBarLength, yBarUnits);

    % Decide on the x and y position for the y scale bar
    xPosYBar = xLimits(1) + xPosNormalized * diff(xLimits) + xBarLength;
    yPosYBar = yLimits(1) + yPosNormalized * diff(yLimits);

    % Decide on the y bar limits
    yLimitsYBar = yPosYBar + [0, yBarLength];

    % Decide on the x and y position for the y bar text
    xPosYText = xPosYBar + bar2TextOffsetNormalized * diff(xLimits);
    yPosYText = yPosYBar + yBarLength / 2;

    % Plot the y scale bar
    yBar = plot_window_boundaries(yLimitsYBar, 'BarValue', xPosYBar, ...
                        'BoundaryType', 'verticalBars', ...
                        'Color', barColor, 'LineWidth', barLineWidth, ...
                        otherArguments);

    % Plot the y bar text
    yBarText = text(xPosYText, yPosYText, yBarLabel);
end

%% Output results
if plotXScaleBar
    handles.xBar = xBar;
    handles.xBarText = xBarText;
end
if plotYScaleBar
    handles.yBar = yBar;
    handles.yBarText = yBarText;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bestLength = compute_best_bar_length(axisLimits, relativeBarLength)
%% Computes the "best" scale bar length

% Compute the approximate bar length
approxBarLength = diff(axisLimits) .* relativeBarLength;

% Round to just one significant digit
bestLength = round(approxBarLength, 1, 'significant');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
