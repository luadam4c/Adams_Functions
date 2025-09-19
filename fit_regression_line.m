function [results] = fit_regression_line (xData, yData, varargin)
%% Fits a linear regression line to x and y data
% Usage: [results] = fit_regression_line (xData, yData, varargin)
% Explanation:
%       This function fits a first-degree polynomial (a straight line) to
%       the provided X and Y data and returns the fitting results.
%
% Example(s):
%       xData = (1:10);
%       yData = xData * 1.5 + 2 * randn(1, 10);
%       results = fit_regression_line(xData, yData, 'ComputeCorrCoeff', true);
%       disp(results);
%
% Outputs:
%       results         - a structure with regression results:
%                           .polyCoeffs: polynomial coefficients [slope, intercept]
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
%       xData       - x data values for the regression
%                   must be a numeric vector
%       yData       - y data values for the regression
%                   must be a numeric vector
%       varargin    - 'ComputeCorrCoeff': whether to compute correlation
%                   must be a logical scalar
%                   default == false
%                   - 'ThroughOrigin': whether to force the intercept to be zero
%                   must be a logical scalar
%                   default == false
%
% Requires:
%       cd/plot_correlation_coefficient.m
%
% Used by:
%       cd/plot_regression_line.m

% File History:
% 2025-09-18 - Created by extracting from plot_regression_line.m

%% Default values for optional arguments
computeCorrCoeffDefault = false;
throughOriginDefault = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'xData', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'yData', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ComputeCorrCoeff', computeCorrCoeffDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'scalar'}));
addParameter(iP, 'ThroughOrigin', throughOriginDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'scalar'}));

% Read from the Input Parser
parse(iP, xData, yData, varargin{:});
computeCorrCoeff = iP.Results.ComputeCorrCoeff;
throughOrigin = iP.Results.ThroughOrigin;

%% Preparation
% Give a warning if user wants correlation for a through-origin regression
if throughOrigin && computeCorrCoeff
    warning(['Computing the standard Pearson correlation coefficient ', ...
             'for a regression forced through the origin can be misleading.']);
end

%% Do the job
% Fit a line to the data
if throughOrigin
    % Forcing the intercept to be zero (y = mx)
    % The slope m is calculated as sum(x*y) / sum(x^2)
    m = xData(:) \ yData(:);
    polyCoeffs = [m, 0];
else
    % Fit a first-degree polynomial (y = mx + b)
    polyCoeffs = polyfit(xData, yData, 1);
end

% Get the slope and intercepts
m = polyCoeffs(1);
b = polyCoeffs(2);

% Evaluate the polynomial at the actual data points
yFitData = polyval(polyCoeffs, xData);

% Calculate the sum of squares
ssResiduals = sum((yData - yFitData).^2);
if throughOrigin
    % R-squared for through-origin is relative to y = 0
    ssTotal = sum(yData.^2);
else
    % Standard R-squared is relative to the mean of y
    ssTotal = sum((yData - mean(yData)).^2);
end

% Calculate R-squared, handle case where ssTotal is 0
if ssTotal == 0
    rSquared = 1;
else
    rSquared = 1 - (ssResiduals / ssTotal);
end

% Save in results
results.polyCoeffs = polyCoeffs;
results.slope = m;
results.intercept = b;
results.ssResiduals = ssResiduals;
results.ssTotal = ssTotal;
results.rSquared = rSquared;

% Compute the correlation coefficient if requested
if computeCorrCoeff
    % Compute without plotting
    [corrText, isSignificant, corrCoeff, pValue] = ...
        plot_correlation_coefficient('XData', xData, 'YData', yData, ...
                                     'ComputeOnly', true);
else
    corrText = ''; isSignificant = NaN; corrCoeff = NaN; pValue = NaN;
end

% Add to results
results.corrText = corrText;
results.isSignificant = isSignificant;
results.corrCoeff = corrCoeff;
results.pValue = pValue;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
