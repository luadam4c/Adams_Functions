function lines = plot_error_bar (xValue, yLow, yHigh, varargin)
%% Plots error bar(s)
% Usage: lines = plot_error_bar (xValue, yLow, yHigh, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       lines       - Line objects plotted
%                   specified as a Line object handle array
% Arguments:
%       xValue      - x value of the error bar(s)
%                   must be a numeric vector
%       yLow        - lower y value of the error bar(s)
%                   must be a numeric vector
%       yHigh       - upper y value of the error bar(s)
%                   must be a numeric vector
%       varargin    - 'BarWidth': bar width(s)
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the line() function
%
% Requires:
%       ~/Adams_Functions/create_error_for_nargin.m
%       ~/Adams_Functions/struct2arglist.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-01-24 Created by Adam and Katerina
% 

%% Hard-coded parameters
relativeBarWidth = 0.25;

%% Default values for optional arguments
barWidthDefault = [];             % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'xValue');
addRequired(iP, 'yLow');
addRequired(iP, 'yHigh');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BarWidth', barWidthDefault);

% Read from the Input Parser
parse(iP, xValue, yLow, yHigh, varargin{:});
barWidth = iP.Results.BarWidth;

% Keep unmatched arguments for the line() function
otherArguments = struct2arglist(iP.Unmatched);

% Check relationships between arguments
if numel(xValue) ~= numel(yLow)
    disp("STOP NOW because number x values does not equal number of y Lows")
    return
elseif numel(xValue) ~= numel(yHigh)
    disp("STOP NOW because number x values does not equal number of y Highs")
    return
end

%% Preparation
% Set default bar width
if isempty(barWidth)
    % Set default
    barWidth = relativeBarWidth * min(diff(xValue));
else
    % Don't overwrite default
end

% Count the number of lines
nLines = numel(xValue);

% Apply force_column_vector to each input argument
[xValue, yLow, yHigh] = argfun(@force_column_vector, xValue, yLow, yHigh);

% Compute the lower and upper x limits for the bars
xLow = xValue - 0.5 * barWidth;
xHigh = xValue + 0.5 * barWidth;

% Set y limits
yLimits = transpose([yLow, yHigh]);

% Set x limits
xLimits = transpose([xLow, xHigh]);

% Preallocate an array for Line objects
lines = gobjects(nLines, 3);

%% Do the job
% Plot the vertical line
lines(:, 1) = plot_vertical_line(xValue, 'YLimits', yLimits, otherArguments{:});

% Plot the upper limit horizontal line
lines(:, 2) = plot_horizontal_line(yHigh, 'XLimits', xLimits, otherArguments{:});

% Plot the lower limit horizontal line
lines(:, 3) = plot_horizontal_line(yLow, 'XLimits', xLimits, otherArguments{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

xValue = force_column_vector(xValue)
yLow = force_column_vector(yLow);
yHigh = force_column_vector(yHigh);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%