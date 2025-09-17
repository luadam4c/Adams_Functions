function [textObjects, isSignificant, corrValue, pValue] = plot_correlation_coefficient (varargin)
%% Plots the correlation coefficient between x and y data
% Usage: [textObjects, isSignificant, corrValue, pValue] = plot_correlation_coefficient (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       xData = (1:10) + 3 * randn(1, 10);
%       yData = (2:2:20) + 3 * randn(1, 10);
%       figure; plot(xData, yData, 'o');
%       plot_correlation_coefficient;
%
% Outputs:
%       textObjects     - correlation coefficient Text object returned
%                       specified as an array of Text objects
%       isSignificant   - whether correlation is deemed significant
%                       specified as a logical scalar
%       corrValue       - correlation coefficient value
%                       specified as a numeric scalar
%       pValue          - p value
%                       specified as a numeric scalar
%
% Arguments:
%       varargin    - 'XData': x data values
%                   must be a a numeric vector
%                   default == detected from current axes
%                   - 'YData': ydata values
%                   must be a a numeric vector
%                   default == detected from current axes
%                   - 'AxesHandle': axes handle to plot on
%                   must be a empty or an axes object handle
%                   default == set in set_axes_properties.m
%                   - Any other parameter-value pair for text()
%
% Requires:
%       cd/apply_over_cells.m
%       cd/extract_fields.m
%       cd/force_column_vector.m
%       cd/hold_off.m
%       cd/hold_on.m
%       cd/struct2arglist.m
%       cd/set_axes_properties.m
%
% Used by:
%       cd/m3ha_plot_grouped_scatter.m
%       cd/plot_relative_events.m
%       cd/virt_analyze_sniff_whisk.m
%       \Shared\Code\vIRt\virt_moore.m

% File History:
% 2020-08-04 Adapted from m3ha_simulation_population.m
% 2020-08-19 Now computes and shows p values
% 2025-09-13 Now ignores NaN values
% 2025-09-17 Added 'AxesHandle' as an optional argument
% TODO: Add regression line

%% Hard-coded parameters

%% Default values for optional arguments
xDataDefault = [];
yDataDefault = [];
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
addParameter(iP, 'AxesHandle', axHandleDefault);

% Read from the Input Parser
parse(iP, varargin{:});
xData = iP.Results.XData;
yData = iP.Results.YData;
axHandle = iP.Results.AxesHandle;

% Keep unmatched arguments for the text() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Find all line objects
if isempty(xData) || isempty(yData)
    lines = findobj(gca, 'Type', 'Line');
end

% Extract x data values
if isempty(xData) 
    xData = extract_fields(lines, 'XData');
    xData = force_column_vector(xData, 'ToLinearize', true);
    xData = apply_over_cells(@vertcat, xData);
end

% Extract y data values
if isempty(yData) 
    yData = extract_fields(lines, 'YData');
    yData = force_column_vector(yData, 'ToLinearize', true);
    yData = apply_over_cells(@vertcat, yData);
end

% If there are NaN values, remove from each
if any(isnan(xData)) || any(isnan(yData))
    toRemove = isnan(xData) | isnan(yData);
    xData(toRemove) = [];
    yData(toRemove) = [];
end

% Decide on the axes to plot on
axHandle = set_axes_properties('AxesHandle', axHandle);

%% Compute
% Compute the 2D correlation coefficient
corrValue = corr2(xData, yData);

% Count the number of data points
nPoints = numel(xData);

% Test the significance
[isSignificant, pValue] = test_corr_significance(corrValue, nPoints);

% Test using corrcoef, should be the same for vector data
% [corrMatrix, pValues] = corrcoef(xData, yData);
% corrValue2 = corrMatrix(2, 1);
% pValue2 = pValues(2, 1);

%% Plot
% Decide on the text color
if isSignificant && abs(corrValue) ~= 1
    textColor = 'r';
else
    textColor = 'k';
end

% Hold on
wasHold = hold_on;

% Plot the correlation coefficient
textObjects(1, 1) = ...
    text(axHandle, 0.1, 0.95, sprintf('Correlation coefficient: %.2f', corrValue), ...
        'Units', 'normalized', 'Color', textColor, otherArguments{:}); 

% Plot the p value
textObjects(2, 1) = ...
    text(axHandle, 0.1, 0.9, sprintf('p value = %.2g', pValue), ...
        'Units', 'normalized', 'Color', textColor, otherArguments{:}); 

% Hold off
hold_off(wasHold);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [isSignificant, pValue] = test_corr_significance (corrValue, nPoints)
%% Tests whether a correlation coefficient is significant

%% Hard-coded parameters
sigLevel = 0.05;

% Compute the degree of freedom
degreeFreedom = (nPoints - 2);

% Compute the t statistic
%   Taken from Equation 34.5 of Analyzing Neural Time Series Data
%       by Mike X Cohen
tStatistic = corrValue * sqrt(degreeFreedom / (1 - corrValue^2));

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