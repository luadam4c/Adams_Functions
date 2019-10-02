function [handles1, handles2] = plot_chevron (data, varargin)
%% Plots a Chevron (paired comparison) plot from data
% Usage: [handles1, handles2] = plot_chevron (data, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       randVec1 = randi(10, 10, 1);
%       randVec2 = randi(10, 10, 1) + 2;
%       data = [randVec1, randVec2];
%       plot_chevron(data)
%
% Outputs:
%       handles1    - TODO: Description of handles1
%                   specified as a TODO
%       handles2    - TODO: Description of handles1
%                   specified as a TODO
%
% Arguments:
%       data        - data table or data vectors
%                   must be a table or a numeric array
%                       or a cell array of numeric vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for plot_tuning_curve()
%
% Requires:
%       ~/Adams_Functions/create_error_for_nargin.m
%       ~/Adams_Functions/struct2arglist.m
%       /TODO:dir/TODO:file
%       cd/argfun.m
%       cd/compute_stats.m
%       cd/create_labels_from_numbers.m
%       cd/force_matrix.m
%       cd/plot_tuning_curve.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-10-01 Created by Adam Lu
% 

%% Hard-coded parameters
readoutLabel = 'SWD count';
lineWidth = 1;
markerColor = [0, 0, 0];
markerSize = 4;
colorMap = [0, 0, 0];
meanLineWidth = 2;
meanMarkSize = 6;
meanColorMap = 'r';
plotStar = false;
figTitle = '';

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'data', ...
    @(x) validateattributes(x, {'numeric', 'cell', 'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
    % TODO: validation function %);

% Read from the Input Parser
parse(iP, data, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the plot_tuning_curve() function
otherArguments = iP.Unmatched;

%% Preparation
% Decide on data values
if istable(data)
    % Extract values
    dataValues = table2array(data);
else
    % Force as a matrix
    dataValues = force_matrix(data, 'AlignMethod', 'leftadjustpad');
end

% Count the number of conditions
nConds = size(dataValues, 2);

% Count the number of samples
nSamples = size(dataValues, 1);

% Compute parameter values for the plot
pValues = transpose(1:nConds);

% Compute parameter axis limits
pLimits = compute_paxis_limits_chevron(pValues);

% Decide on p tick labels
if istable(data)
    % Extract variable names if any
    pTickLabels = data.Properties.VariableNames;
else
    % Create labels
    pTickLabels = create_labels_from_numbers(pValues, 'Prefix', 'time');
end

% Decide on column labels
if istable(data) && ~isempty(data.Properties.RowNames)
    % Extract row names if any
    columnLabels = data.Properties.RowNames;
else
    % Create labels
    columnLabels = create_labels_from_numbers(1:nSamples, 'Prefix', 'data');
end

% Compute means and confidence intervals
if nConds == 2
    % Compute the mean of the baseline values
    baseMean = compute_stats(dataValues(:, 1), 'mean');

    % Compute the differences
    diffValues = dataValues(:, 2) - dataValues(:, 1);

    % Compute the mean and confidence intervals of the differences
    [diffMean, diffLower95, diffUpper95] = ...
        argfun(@(x) compute_stats(diffValues, x), 'mean', 'lower95', 'upper95');

    % Compute the values to plot
    [meanValues, lower95Values, upper95Values] = ...
        argfun(@(x) baseMean + [0; x], diffMean, diffLower95, diffUpper95);
end

%% Do the job
% Plot a tuning curve
handles1 = plot_tuning_curve(pValues, transpose(dataValues), ...
                    'PLimits', pLimits, 'PTicks', pValues, ...
                    'PTickLabels', pTickLabels, ...
                    'ColumnLabels', columnLabels, ...
                    'PLabel', 'suppress', 'ReadoutLabel', readoutLabel, ...
                    'FigTitle', figTitle, ...
                    'RunTTest', true, 'RunRankTest', true, ...
                    'LineWidth', lineWidth, 'ColorMap', colorMap, ...
                    'Marker', 'o', 'MarkerFaceColor', colorMap, ...
                    'MarkerSize', markerSize, ...
                    'LegendLocation', 'suppress', otherArguments);

% Hold on
if ~ishold
    wasHold = false;
    hold on
else
    wasHold = true;
end

% Plot the mean and confidence intervals of the differences
if nConds == 2
    handles2 = plot_tuning_curve(pValues, meanValues, ...
                    'LowerCI', lower95Values, 'UpperCI', upper95Values, ...
                    'PlotCurveOnly', true, ...
                    'LineWidth', meanLineWidth, 'ColorMap', meanColorMap, ...
                    'Marker', 'o', 'MarkerFaceColor', meanColorMap, ...
                    'MarkerSize', meanMarkSize);
end

% TODO: Place in plot_tuning_curve 'PlotStar'
% Plot a star if significant
if plotStar
    starXPos = mean(pValues);
    yLimits = get(gca, 'YLim');
    starYPos = yLimits(1) + 0.8 * (yLimits(2) - yLimits(1));
    plot(starXPos, starYPos, 'k*');
end

% Hold off
if ~wasHold
    hold off
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xLimits = compute_paxis_limits_chevron (tickValues)
%% Computes x axis limits for Chevron plots

% Hard-coded parameters
marginPercentage = 25;          % 25% margins

% Compute the mean difference between ticks
meanTickInt = mean(diff(tickValues));

% Compute x axis limits
xLimits = [tickValues(1), tickValues(end)] + ...
            (marginPercentage / 100) * meanTickInt * [-1, 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Compute the mean, lower and upper confidence interval bounds
[means, lower95s, upper95s] = ...
    argfun(@(x) compute_stats(dataValues, x), 'mean', 'lower95', 'upper95');

% Plot the means
plot(pValues, means, 'r-o', ...
    'LineWidth', meanLineWidth, 'MarkerSize', meanMarkSize, ...
    'MarkerFaceColor', meanColorMap);

% Plot error bars
plot_error_bar(pValues, lower95s, upper95s, 'Color', meanColorMap, ...
                'LineWidth', meanLineWidth);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%