function [hRegression, hText, textStr, results] = plot_regression_line (varargin)
%% Plots a linear regression line for x and y data
% Usage: [hRegression, hText, textStr, results] = plot_regression_line (varargin)
% Explanation:
%       This function fits a first-degree polynomial (a straight line) to
%       the provided X and Y data and plots it on a specified or current axes.
%       It can automatically extract data from line objects already on the axes.
%
% Example(s):
%       xData = (1:10);
%       yData = xData * 1.5 + 2 * randn(1, 10);
%       figure; plot(xData, yData, 'o');
%       plot_regression_line;
%       figure; plot(xData, yData, 'ks');
%       plot_regression_line('Color', 'r', 'LineWidth', 2, 'LineStyle', '-');
%       plot_regression_line('ShowRSquared', false, 'TextLocation', 'topleft', 'TextMargin', 0.1);
%       figure; plot(xData, yData, 'bo');
%       [hLine, hEq] = plot_regression_line('ShowEquation', true);
%       [hLine, hEq] = plot_regression_line('ShowEquation', true, 'ThroughOrigin', true);
%       [~, ~, ~, res] = plot_regression_line('ShowCorrCoeff', true, 'Delimiter', '  ');
%       disp(res);
%
% Outputs:
%       hRegression     - handle to the regression line
%                       specified as a Line object
%       hText           - handle to the equation/statistics text object
%                       specified as a Text object
%       textStr         - equation and/or R-squared of the regression line
%                       specified as a character array
%       results         - a structure with regression results:
%                           .slope: slope of the line
%                           .intercept: y-intercept of the line
%                           .ssResiduals: sum of squared residuals
%                           .ssTotal: total sum of squares
%                           .rSquared: R-squared value
%                           .isSignificant: whether correlation is significant
%                           .corrCoeff: Pearson correlation coefficient
%                           .pValue: p-value of the correlation
%                       specified as a structure
%
% Arguments:
%       varargin    - 'XData': x data values for the regression
%                   must be a a numeric vector
%                   default == detected from current axes
%                   - 'YData': y data values for the regression
%                   must be a a numeric vector
%                   default == detected from current axes
%                   - 'ShowEquation': whether to show the equation
%                   must be a logical scalar
%                   default == true
%                   - 'ShowRSquared': whether to show the R-squared value
%                   must be a logical scalar
%                   default == true
%                   - 'ShowCorrCoeff': whether to show the correlation coefficient
%                   must be a logical scalar
%                   default == false
%                   - 'ThroughOrigin': whether to force the intercept to be zero
%                   must be a logical scalar
%                   default == false
%                   - 'Delimiter': Delimiter for joining text strings
%                   must be a character vector or string scalar
%                   default == '\n' (newline)
%                   - 'Units': coordinate system for text positioning
%                   must be an unambiguous, case-insensitive match to one of:
%                       'normalized'  - relative to the axes (default)
%                       'data'        - uses the data coordinates of the plot
%                       'pixels'      - uses pixel units
%                       'inches'      - uses inches
%                       'centimeters' - uses centimeters
%                       'points'      - uses points
%                   default == 'normalized'
%                   - 'TextLocation': location for the text
%                   must be an unambiguous, case-insensitive match to one of:
%                       'topleft'     - places text in the top-left corner
%                       'topright'    - places text in the top-right corner
%                       'bottomleft'  - places text in the bottom-left corner
%                       'bottomright' - places text in the bottom-right corner
%                   default == 'bottomright'
%                   - 'TextMargin': textMargin for the text location
%                   must be a numeric scalar
%                   default == 0.05
%                   - 'HorizontalAlignment': Horizontal alignment of the text
%                   must be an unambiguous, case-insensitive match to one of:
%                       'left'      - Aligns the left edge of the text with the position
%                       'center'    - Centers the text at the position
%                       'right'     - Aligns the right edge of the text with the position
%                   default == set based on 'TextLocation'
%                   - 'VerticalAlignment': Vertical alignment of the text
%                   must be an unambiguous, case-insensitive match to one of:
%                       'top'       - Aligns the top edge of the text with the position
%                       'cap'       - Aligns with characters that have capital letters
%                       'middle'    - Centers the text vertically at the position
%                       'baseline'  - Aligns with the text baseline
%                       'bottom'    - Aligns the bottom edge of the text with the position
%                   default == set based on 'TextLocation'
%                   - 'Color': color for the line and text
%                   must be a valid color specification (e.g., 'r', [1 0 0])
%                   default == 'r' for significant, 'k' for non-significant
%                   - 'LineStyle': style for the regression line
%                   must be an unambiguous, case-insensitive match to one of:
%                       '-'     - solid line
%                       '--'    - dashed line
%                       ':'     - dotted line
%                       '-.'    - dash-dotted line
%                       'none'  - no line
%                   default == '--'
%                   - 'AxesHandle': axes handle to plot on
%                   must be a empty or an axes object handle
%                   default == gca
%                   - Any other parameter-value pair for plot()
%
% Requires:
%       cd/extract_data_from_lines.m
%       cd/fit_regression_line.m
%       cd/hold_off.m
%       cd/hold_on.m
%       cd/islinestyle.m
%       cd/plot_text.m
%       cd/struct2arglist.m
%       cd/set_axes_properties.m
%
% Used by:
%       cd/plot_correlation_coefficient.m
%       cd/virt_plot_amplitude_correlation.m
%       cd/virt_plot_phase_response.m

% File History:
% 2025-09-17 - Extracted from plot_correlation_coefficient.m
% 2025-09-17 - Added 'ShowEquation' argument and hText output.
% 2025-09-18 - Added 'ShowRSquared', 'ShowCorrCoeff', 'ThroughOrigin' arguments.
% 2025-09-18 - R-squared and equation text are now combined into one object.
% 2025-09-18 - Added 'TextLocation', 'Color', 'LineStyle' args & updated results.
% 2025-09-18 - Added 'TextMargin' as an optional argument.
% 2025-09-18 - Now uses plot_text.m for displaying text.
% 2025-09-18 - Added more plot_text passthrough arguments.
% 2025-09-18 - Now uses fit_regression_line.m for fitting.

%% Hard-coded parameters
validTextLocations = {'topleft', 'topright', 'bottomleft', 'bottomright'};
validUnits = {'normalized', 'data', 'pixels', 'inches', 'centimeters', 'points'};

%% Default values for optional arguments
xDataDefault = [];
yDataDefault = [];
showEquationDefault = true;
showRSquaredDefault = true;
showCorrCoeffDefault = false;
throughOriginDefault = false;
delimiterDefault = '\n';
unitsDefault = 'normalized';
textLocationDefault = 'bottomright';
textMarginDefault = 0.05;
hAlignDefault = '';
vAlignDefault = '';
colorDefault = [];              % determined by significance by default
lineStyleDefault = '--';
axHandleDefault = [];           % gca by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'XData', xDataDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'YData', yDataDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'ShowEquation', showEquationDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'scalar'}));
addParameter(iP, 'ShowRSquared', showRSquaredDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'scalar'}));
addParameter(iP, 'ShowCorrCoeff', showCorrCoeffDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'scalar'}));
addParameter(iP, 'ThroughOrigin', throughOriginDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'scalar'}));
addParameter(iP, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Units', unitsDefault, ...
    @(x) any(validatestring(x, validUnits)));
addParameter(iP, 'TextLocation', textLocationDefault, ...
    @(x) any(validatestring(x, validTextLocations)));
addParameter(iP, 'TextMargin', textMarginDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'HorizontalAlignment', hAlignDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'VerticalAlignment', vAlignDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Color', colorDefault);
addParameter(iP, 'LineStyle', lineStyleDefault, ...
    @(x) all(islinestyle(x, 'ValidateMode', true)));
addParameter(iP, 'AxesHandle', axHandleDefault);

% Read from the Input Parser
parse(iP, varargin{:});
xData = iP.Results.XData;
yData = iP.Results.YData;
showEquation = iP.Results.ShowEquation;
showRSquared = iP.Results.ShowRSquared;
showCorrCoeff = iP.Results.ShowCorrCoeff;
throughOrigin = iP.Results.ThroughOrigin;
delimiter = iP.Results.Delimiter;
units = validatestring(iP.Results.Units, validUnits);
textLocation = validatestring(iP.Results.TextLocation, validTextLocations);
textMargin = iP.Results.TextMargin;
hAlign = iP.Results.HorizontalAlignment;
vAlign = iP.Results.VerticalAlignment;
colorUser = iP.Results.Color;
[~, lineStyle] = islinestyle(iP.Results.LineStyle, 'ValidateMode', true);
axHandle = iP.Results.AxesHandle;

% Keep unmatched arguments for the plot() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Decide on the axes to plot on
axHandle = set_axes_properties('AxesHandle', axHandle);

% Get data for regression. This function handles extracting from the plot,
% overriding with any user-provided data, and removing NaNs.
[xData, yData] = extract_data_from_lines('AxesHandle', axHandle, ...
                                         'XData', xData, ...
                                         'YData', yData, ...
                                         'RemoveNaN', true);

% Initialize outputs
hText = gobjects;
textStr = '';

% Don't plot if there is no data
if isempty(xData)
    hRegression = gobjects;
    results = struct;
    return;
end

%% Fitting
% Fit a line to the data and get all relevant statistics
results = fit_regression_line(xData, yData, ...
                              'ThroughOrigin', throughOrigin, ...
                              'ComputeCorrCoeff', showCorrCoeff);

%% Plotting
% Extract from results structure
polyCoeffs = results.polyCoeffs;
slope = results.slope;
intercept = results.intercept;
rSquared = results.rSquared;
corrText = results.corrText;
corrCoeff = results.corrCoeff;
isSignificant = results.isSignificant;

% Decide on the color to use
if ~isempty(colorUser)
    colorMap = colorUser;
elseif ~isnan(isSignificant) && isSignificant && ...
        ~isnan(corrCoeff) && abs(corrCoeff) ~= 1
    colorMap = 'r';
else
    colorMap = 'k';
end

% Hold on
wasHold = hold_on(axHandle);

% Get the x-axis limits to define the line's extent
xLimits = get(axHandle, 'XLim');

% Evaluate the polynomial at the limits to get the y-values for the line
yFitLine = polyval(polyCoeffs, xLimits);

% Plot the regression line
hRegression = plot(axHandle, xLimits, yFitLine, ...
                   'Color', colorMap, 'LineStyle', lineStyle, otherArguments{:});

% Prepare text to display on the plot
textStrings = {};

% Create the equation string if requested
if showEquation
    if throughOrigin
        eqStr = sprintf('y = %.2fx', slope);
    elseif intercept >= 0
        eqStr = sprintf('y = %.2fx + %.2f', slope, intercept);
    else
        eqStr = sprintf('y = %.2fx - %.2f', slope, abs(intercept));
    end
    textStrings{end+1} = eqStr;
end

% Create the R-squared string if requested
if showRSquared
    if throughOrigin
        rSquaredStr = sprintf('Uncentered R^2 = %.2f', rSquared);
    else
        rSquaredStr = sprintf('R^2 = %.2f', rSquared);
    end
    textStrings{end+1} = rSquaredStr;
end

% Show the correlation coefficient if requested
if showCorrCoeff && ~isempty(corrText)
    % Add to text strings
    textStrings{end+1} = corrText;
end

% Plot the text string if it's not empty
if ~isempty(textStrings)
    hText = plot_text(textStrings, 'AxesHandle', axHandle, ...
                      'Delimiter', delimiter, ...
                      'Units', units, ...
                      'TextLocation', textLocation, ...
                      'TextMargin', textMargin, ...
                      'HorizontalAlignment', hAlign, ...
                      'VerticalAlignment', vAlign, ...
                      'Color', colorMap);
end

% Hold off
hold_off(wasHold, axHandle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%