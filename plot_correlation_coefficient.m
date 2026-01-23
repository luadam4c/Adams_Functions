function [textObjectOrString, isSignificant, corrCoeff, pValue, hRegression] = plot_correlation_coefficient (varargin)
%% Computes and plots the correlation coefficient and p-value for X-Y data
% Usage: [textObjectOrString, isSignificant, corrCoeff, pValue, hRegression] = plot_correlation_coefficient (varargin)
% Explanation:
%       Computes and plots the Pearson correlation coefficient and its p-value
%       for a given set of X-Y data using the corr2 function (from the Image Processing Toolbox). 
%       The p-value is computed using a t-test with n-2 degrees of freedom, 
%       where the t-statistic is defined as t = r * sqrt((n-2)/(1-r^2)) 
%       (Cohen, Mike X. Analyzing Neural Time Series Data, Eq 34.5).
%       Note that this is the exact same behavior as built in corrcoef()
%
%       Data can be provided directly or extracted from an existing plot. 
%       Optionally, it can also plot the linear regression line.
%
% Example(s):
%       xData = (1:10) + 3 * randn(1, 10);
%       yData = (2:2:20) + 3 * randn(1, 10);
%       figure; plot(xData, yData, 'o');
%       plot_correlation_coefficient;
%       plot_correlation_coefficient('PlotRegression', true, 'Color', 'blue');
%       plot_correlation_coefficient('TextLocation', 'bottomright', 'TextMargin', 0.1);
%       corrText = plot_correlation_coefficient('ComputeOnly', true);
%
% Outputs:
%       textObjectOrString - handle to text object or the text string itself
%                           specified as a text object handle or a character vector
%       isSignificant   - whether correlation is deemed significant
%                       specified as a logical scalar
%       corrCoeff       - correlation coefficient value
%                       specified as a numeric scalar
%       pValue          - p value
%                       specified as a numeric scalar
%       hRegression     - handle to the regression line
%                       specified as Line object
%
% Arguments:
%       varargin    - 'XData': x data values
%                   must be a a numeric vector
%                   default == detected from current axes
%                   - 'YData': ydata values
%                   must be a a numeric vector
%                   default == detected from current axes
%                   - 'PlotRegression': whether to plot a regression line
%                   must be a logical scalar
%                   default == false
%                   - 'ShowPValue': whether to show the p-value
%                   must be a logical scalar
%                   default == true
%                   - 'ComputeOnly': whether to compute only and not plot
%                   must be a logical scalar
%                   default == false
%                   - 'TextLocation': location for the text
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'topleft'     - places text in the top-left corner
%                       'topright'    - places text in the top-right corner
%                       'bottomleft'  - places text in the bottom-left corner
%                       'bottomright' - places text in the bottom-right corner
%                   default == 'topleft'
%                   - 'TextMargin': textMargin for the text location
%                   must be a numeric scalar
%                   default == 0.05
%                   - 'Color': text color
%                   must be a valid color specification (e.g., 'r', [1 0 0])
%                   default == 'r' for significant, 'k' for non-significant
%                   - 'AxesHandle': axes handle to plot on
%                   must be a empty or an axes object handle
%                   default == set in set_axes_properties.m
%                   - Any other parameter-value pair for plot_text() or text()
%
% Requires:
%       cd/extract_data_from_lines.m
%       cd/plot_regression_line.m
%       cd/plot_text.m
%       cd/hold_off.m
%       cd/hold_on.m
%       cd/struct2arglist.m
%       cd/set_axes_properties.m
%
% Used by:
%       cd/m3ha_plot_grouped_scatter.m
%       cd/fit_regression_line.m
%       cd/plot_regression_line.m
%       cd/plot_relative_events.m
%       cd/virt_analyze_sniff_whisk.m

% File History:
% 2020-08-04 Adapted from m3ha_simulation_population.m
% 2020-08-19 Now computes and shows p values
% 2025-09-13 Now ignores NaN values
% 2025-09-17 Added 'AxesHandle' as an optional argument
% 2025-09-17 Added 'PlotRegression' as an optional argument
% 2025-09-17 Updated to use extract_data_from_lines.m and plot_regression_line.m
% 2025-09-18 Combined text outputs; added 'ComputeOnly' and 'TextLocation' args
% 2025-09-18 Added 'ShowPValue' as an optional argument
% 2025-09-18 Added 'TextMargin' and 'Color' as optional arguments
% 2025-09-18 Now uses plot_text.m for displaying text.
% TODO: If Image Processing Toolbox not installed, use corrcoef() instead and display warning

%% Hard-coded parameters
validTextLocations = {'topleft', 'topright', 'bottomleft', 'bottomright'};
sigLevel = 0.05;                % significance level

%% Default values for optional arguments
xDataDefault = [];
yDataDefault = [];
plotRegressionDefault = false;
showPValueDefault = true;
computeOnlyDefault = false;
textLocationDefault = 'topleft';
textMarginDefault = 0.05;
colorDefault = [];              % determined by significance by default
axHandleDefault = [];           % gca by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

addParameter(iP, 'XData', xDataDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'YData', yDataDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'PlotRegression', plotRegressionDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'scalar'}));
addParameter(iP, 'ShowPValue', showPValueDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'scalar'}));
addParameter(iP, 'ComputeOnly', computeOnlyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'scalar'}));
addParameter(iP, 'TextLocation', textLocationDefault, ...
    @(x) any(validatestring(x, validTextLocations)));
addParameter(iP, 'TextMargin', textMarginDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Color', colorDefault);
addParameter(iP, 'AxesHandle', axHandleDefault);

% Read from the Input Parser
parse(iP, varargin{:});
xData = iP.Results.XData;
yData = iP.Results.YData;
plotRegression = iP.Results.PlotRegression;
showPValue = iP.Results.ShowPValue;
computeOnly = iP.Results.ComputeOnly;
textLocation = validatestring(iP.Results.TextLocation, validTextLocations);
textMargin = iP.Results.TextMargin;
colorUser = iP.Results.Color;
axHandle = iP.Results.AxesHandle;

% Keep unmatched arguments for the text() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Decide on the axes to plot on
if ~computeOnly
    axHandle = set_axes_properties('AxesHandle', axHandle);
end

% Get data for correlation. This handles extraction from the plot,
% overrides with user-provided data, and removes NaNs.
[xData, yData] = extract_data_from_lines('AxesHandle', axHandle, ...
                                         'XData', xData, ...
                                         'YData', yData, ...
                                         'RemoveNaN', true);

%% Compute
% Return empty if no data
if isempty(xData)
    if computeOnly
        textObjectOrString = '';
    else
        textObjectOrString = gobjects;
    end
    hRegression = gobjects; 
    isSignificant = false; corrCoeff = NaN; pValue = NaN;
    return;
end

% Compute the 2D correlation coefficient
corrCoeff = corr2(xData, yData);

% Count the number of data points
nPoints = numel(xData);

% Test the significance
[isSignificant, pValue] = test_corr_significance(corrCoeff, nPoints, sigLevel);

% Test using corrcoef, should be the same for vector data
% [corrMatrix, pValues] = corrcoef(xData, yData);
% corrValue2 = corrMatrix(2, 1);
% pValue2 = pValues(2, 1);

% Create the text string for the correlation coefficient and p value
textCell = {sprintf('Correlation coefficient: %.2f', corrCoeff)};
if showPValue
    textCell{end+1} = sprintf('p value = %.2g', pValue);
end
textStr = strjoin(textCell, '\n');

% Return text if compute only
if computeOnly
    textObjectOrString = textStr;
    hRegression = gobjects; % Ensure graphics handle is empty
    return;
end

%% Plot
% Initialize regression line handle
hRegression = gobjects;

% Decide on the text color
if ~isempty(colorUser)
    colorMap = colorUser;
elseif isSignificant && abs(corrCoeff) ~= 1
    colorMap = 'r';
else
    colorMap = 'k';
end

% Hold on
wasHold = hold_on(axHandle);

% Plot the correlation text using plot_text
textObjectOrString = ...
    plot_text(textCell, 'AxesHandle', axHandle, ...
            'TextLocation', textLocation, 'TextMargin', textMargin, ...
            'Color', colorMap, otherArguments{:});

% Plot the regression line if requested
if plotRegression
    hRegression = plot_regression_line('XData', xData, 'YData', yData, ...
                                   'AxesHandle', axHandle, 'Color', colorMap);
end

% Hold off
hold_off(wasHold, axHandle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [isSignificant, pValue] = test_corr_significance (corrCoeff, nPoints, sigLevel)
%% Tests whether a correlation coefficient is significant

% Compute the degree of freedom
degreeFreedom = (nPoints - 2);

% Compute the t statistic
%   Taken from Equation 34.5 of Analyzing Neural Time Series Data
%       by Mike X Cohen
tStatistic = corrCoeff * sqrt(degreeFreedom / (1 - corrCoeff^2));

% Compute a p value
%   Note: 1 - tcdf(x, nu) is inaccurate when tcdf(x, nu) is too close to 1
pValue = tcdf(abs(tStatistic), degreeFreedom, 'upper') * 2;

% Decide whether the correlation coefficient is significant
isSignificant = pValue < sigLevel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%