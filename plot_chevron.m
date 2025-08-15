function [handles, handlesMean] = plot_chevron (data, varargin)
%% Plots a Chevron (paired comparison) plot from data
% Usage: [handles, handlesMean] = plot_chevron (data, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       randVec1 = randi(10, 10, 1);
%       randVec2 = randi(10, 10, 1) + 10;
%       data = [randVec1, randVec2];
%       plot_chevron(data)
%
% Outputs:
%       handles     - handles to plotted objects for the data
%                   specified as a structure
%       handlesMean - handles to plotted objects for the mean difference
%                   specified as a structure
%
% Arguments:
%       data        - data vectors or data table or path to data table
%                   Note: The dimension with fewer elements is taken as 
%                           the parameter
%                   must be a table or a numeric array
%                       or a cell array of numeric vectors
%                       or a string scalar or a character vector
%       varargin    - 'IsLog2Data': whether data is log 2-scaled
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotMeanValues': whether to plot the mean values
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotMeanDifference': whether to plot the mean difference
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotErrorBars': whether to plot error bars 
%                                       even though it's wrong to do so
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RunTTest': whether to run paired t-test
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'RunRankTest': whether to run paired 
%                                       Wilcoxon signed-rank test
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PLimits': limits of parameter axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == set in compute_paxis_limits_chevron()
%                   - 'PTicks': x tick values for the parameter values
%                   must be a numeric vector
%                   default == all parameter values
%                   - 'PTickLabels': x tick labels in place of parameter values
%                   must be a cell array of character vectors/strings
%                   default == table column names or time1, time2, ...
%                   - 'PLabel': label for the parameter, 
%                               suppress by setting value to {'suppress'}
%                   must be a string scalar or a character vector
%                   default == 'suppress'
%                   - 'ColumnLabels': labels for the readout columns, 
%                               suppress by setting value to {'suppress'}
%                   must be a scalartext 
%                       or a cell array of strings or character vectors
%                   default == table row names or data1, data2, ...
%                   - 'ColorMap': color map used when nColumnsToPlot > 1
%                   must be a 2-D numeric array with 3 columns
%                   default == set in decide_on_colormap.m
%                   - 'MeanColorMap': color map for the mean
%                   must be a 2-D numeric array with 3 columns
%                   default == 'r' if colorMap is black; 'k' otherwise
%                   - 'Marker': marker type
%                   must be empty or a positive scalar
%                   default == 'o'
%                   - 'MarkerSize': marker size for individual data
%                   must be empty or a positive scalar
%                   default == 6
%                   - 'LineWidth': line width for individual data
%                   must be empty or a positive scalar
%                   default == 1
%                   - 'LegendLocation': location for legend
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'auto'      - use default
%                       'suppress'  - no legend
%                       anything else recognized by the legend() function
%                   default == 'eastoutside'
%                   - 'FigExpansion': expansion factor for figure position
%                   must be a must be a positive scalar or 2-element vector
%                   default == [1, 0.4]
%                   - 'AxesHandle': axes handle to plot
%                   must be a empty or a axes object handle
%                   default == set in set_axes_properties.m
%                   - 'AlwaysNew': whether to always create a new figure even if
%                                   figNumber is not passed in
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true if data is a path, false otherwise
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - Any other parameter-value pair for plot_tuning_curve()
%
% Requires:
%       cd/argfun.m
%       cd/compute_stats.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/decide_on_colormap.m
%       cd/extract_fileparts.m
%       cd/force_data_as_matrix.m
%       cd/hold_off.m
%       cd/hold_on.m
%       cd/plot_error_bar.m
%       cd/plot_tuning_curve.m
%       cd/save_all_figtypes.m
%       cd/set_axes_properties.m
%       cd/set_figure_properties.m
%
% Used by:
%       cd/m3ha_oscillations_analyze.m
%       cd/metabolismR01_plot_chevrons.m
%       cd/plot_measures.m
%       cd/plot_relative_events.m
%       cd/plot_small_chevrons.m

% File History:
% 2019-10-01 Created by Adam Lu
% 2019-10-03 Made many things optional arguments
% 2019-10-07 Now plots means with open circles
% 2019-10-08 Added 'IsLog2Data' as an optional argument
% 2019-11-27 Added 'PlotMeanValues' as an optional argument
% 2019-11-27 Now adds transparency to mean circles
% 2019-11-27 Added 'MeanColorMap' as an optional argument
% 2020-02-29 Now allows the first argument to be a string
% 2020-02-29 Added 'AlwaysNew' as an optional argument
% 2020-02-29 Added 'FigName' as an optional argument
% TODO: Use this in plot_table.m?

%% Hard-coded parameters

% TODO: Make optional arguments
meanLineWidthRatio = 1.5;
meanLineStyle = '-';
meanMarkerSizeRatio = 2;
meanColorMap = [];                  % set later
meanMarkerFaceColor = 'flat';
meanMarkerFaceAlpha = 0.25;

%% Default values for optional arguments
isLog2DataDefault = false;
plotMeanValuesDefault = true;
plotMeanDifferenceDefault = false;
plotErrorBarsDefault = false;
runTTestDefault = true;
runRankTestDefault = true;
pLimitsDefault = [];                % set later
pTicksDefault = [];                 % set later
pTickLabelsDefault = {};            % set later
pLabelDefault = 'suppress';
columnLabelsDefault = '';           % set later
colorMapDefault = [];               % set later
meanColorMapDefault = [];           % set later
markerDefault = 'o';
markerSizeDefault = 6;
lineWidthDefault = 1;
legendLocationDefault = 'eastoutside';
figExpansionDefault = [];
axHandleDefault = [];               % gca by default
alwaysNewDefault = [];              % set later
figNameDefault = '';                % don't save figure by default
figTypesDefault = {'png', 'epsc'};

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
    @(x) validateattributes(x, {'numeric', 'cell', 'table', ...
                            'char', 'string'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'IsLog2Data', isLog2DataDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotMeanValues', plotMeanValuesDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotMeanDifference', plotMeanDifferenceDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotErrorBars', plotErrorBarsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RunTTest', runTTestDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RunRankTest', runRankTestDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PLimits', pLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'PTicks', pTicksDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'PTickLabels', pTickLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'PLabel', pLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ColumnLabels', columnLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ColorMap', colorMapDefault);
addParameter(iP, 'MeanColorMap', meanColorMapDefault);
addParameter(iP, 'Marker', markerDefault);
addParameter(iP, 'MarkerSize', markerSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'LineWidth', lineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'LegendLocation', legendLocationDefault, ...
    @(x) all(islegendlocation(x, 'ValidateMode', true)));
addParameter(iP, 'FigExpansion', figExpansionDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive'}));
addParameter(iP, 'AxesHandle', axHandleDefault);
addParameter(iP, 'AlwaysNew', alwaysNewDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, data, varargin{:});
isLog2Data = iP.Results.IsLog2Data;
plotMeanValues = iP.Results.PlotMeanValues;
plotMeanDifference = iP.Results.PlotMeanDifference;
plotErrorBars = iP.Results.PlotErrorBars;
runTTest = iP.Results.RunTTest;
runRankTest = iP.Results.RunRankTest;
pLimits = iP.Results.PLimits;
pTicks = iP.Results.PTicks;
pTickLabels = iP.Results.PTickLabels;
pLabel = iP.Results.PLabel;
columnLabels = iP.Results.ColumnLabels;
colorMap = iP.Results.ColorMap;
meanColorMap = iP.Results.MeanColorMap;
marker = iP.Results.Marker;
markerSize = iP.Results.MarkerSize;
lineWidth = iP.Results.LineWidth;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);
figExpansion = iP.Results.FigExpansion;
axHandle = iP.Results.AxesHandle;
alwaysNew = iP.Results.AlwaysNew;
figName = iP.Results.FigName;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the plot_tuning_curve() function
otherArguments = iP.Unmatched;

%% Preparation
% Set default flags
if isempty(alwaysNew)
    if ischar(data) || isstring(data)
        alwaysNew = true;
        if isempty(figName)
            figName = extract_fileparts(data, 'pathbase');
        end
    else
        alwaysNew = false;
    end
end

% Force data values as a numeric matrix
%   where each group is a column and each row is a sample
[dataValues, groupLabelsAuto, sampleLabelsAuto] = force_data_as_matrix(data);

% Count the number of groups
nGroups = size(dataValues, 2);

% Count the number of samples
nSamples = size(dataValues, 1);

% Decide on a color map
colorMap = decide_on_colormap(colorMap, nSamples);

% Decide on the mean color map
if isempty(meanColorMap)
    if all(all(colorMap == 0))
        meanColorMap = decide_on_colormap('r', 1);
    else
        meanColorMap = [0, 0, 0];
    end
else
    meanColorMap = decide_on_colormap(meanColorMap, 1);
end

% Compute parameter values for the plot
pValues = transpose(1:nGroups);

% Compute parameter axis limits
if isempty(pLimits)
    pLimits = compute_paxis_limits_chevron(pValues);
end

% Decide on parameter tick values
if isempty(pTicks)
    pTicks = pValues;
end

% Decide on parameter tick labels
if isempty(pTickLabels)
    pTickLabels = groupLabelsAuto;
end

% Decide on column labels
if isempty(columnLabels)
    columnLabels = sampleLabelsAuto;
end

% Decide on the mean marker size
meanMarkerSize = ceil(meanMarkerSizeRatio * markerSize);
meanLineWidth = meanLineWidthRatio * lineWidth;

%% Do the job
% Decide on the figure
figHandle = set_figure_properties('FigExpansion', figExpansion, ...
                                    'AlwaysNew', alwaysNew);

% Decide on the axes
axHandle = set_axes_properties('AxesHandle', axHandle, 'FigHandle', figHandle);

% Plot a tuning curve
handles = plot_tuning_curve(pValues, transpose(dataValues), ...
                    'RunTTest', runTTest, 'RunRankTest', runRankTest, ...
                    'PLimits', pLimits, 'PTicks', pTicks, ...
                    'PTickLabels', pTickLabels, 'PLabel', pLabel, ...
                    'ColumnLabels', columnLabels, ...
                    'ColorMap', colorMap, ...
                    'LegendLocation', legendLocation, ...
                    'AxesHandle', axHandle, ...                    
                    'LineWidth', lineWidth, ...
                    'Marker', marker, 'MarkerSize', markerSize, ...
                    otherArguments);

% Plot the mean and confidence intervals of the differences
%   TODO: Does it ever make sense for more than 2 groups?
if plotMeanDifference && nGroups == 2
    % Compute the mean of the baseline values
    baseMean = compute_stats(dataValues(:, 1), 'mean');

    % Compute the differences
    diffValues = dataValues(:, 2) - dataValues(:, 1);

    % Compute the mean and confidence intervals of the differences
    [diffMean, diffLower95, diffUpper95] = ...
        argfun(@(x) compute_stats(diffValues, x, 'IgnoreNan', true), ...
                'mean', 'lower95', 'upper95');

    % Compute the values to plot
    [meanValues, lower95Values, upper95Values] = ...
        argfun(@(x) baseMean + [0; x], diffMean, diffLower95, diffUpper95);

    % Hold on
    wasHold = hold_on;

    % Plot mean and confidence intervals
    handlesMean = plot_tuning_curve(pValues, meanValues, 'PlotOnly', true, ...
                    'LowerCI', lower95Values, 'UpperCI', upper95Values, ...
                    'LineWidth', meanLineWidth, 'ColorMap', meanColorMap, ...
                    'LineStyle', meanLineStyle, ...
                    'Marker', 'none', 'AxesHandle', axHandle);

    % Hold off
    hold_off(wasHold);
else
    handlesMean = struct;
end

% Plot error bars if requested
% TODO: Plot mean circle and error bar without overlap
if plotErrorBars
    % Compute the mean, lower and upper confidence interval bounds
    [lower95s, upper95s] = ...
        argfun(@(x) compute_stats(dataValues, x, 'IgnoreNan', true), ...
                'lower95', 'upper95');

    % Hold on
    wasHold = hold_on;

    % Plot error bars
    plot_error_bar(pValues, lower95s, upper95s, 'Color', meanColorMap, ...
                    'LineWidth', meanLineWidth, 'AxesHandle', axHandle);

    % Hold off
    hold_off(wasHold);
end

% Plot mean circles if requested
if plotMeanValues
    % Compute the mean, lower and upper confidence interval bounds
    means = compute_stats(dataValues, 'mean', 'IgnoreNan', true);

    % Hold on
    wasHold = hold_on;

    % Plot the means with circles with transparency
    %   Note: Use scatter to set MarkerFaceAlpha
    scatter(axHandle, pValues, means, meanMarkerSize^2, meanColorMap, ...
            'Marker', 'o', 'LineWidth', meanLineWidth, ...
            'MarkerFaceColor', meanMarkerFaceColor, ...
            'MarkerFaceAlpha', meanMarkerFaceAlpha, ...
            'MarkerEdgeColor', 'none');

    % Plot the means with circles with transparency
    %   Note: Use scatter to set MarkerFaceAlpha
    plot(axHandle, pValues, means, 'Color', meanColorMap, ...
            'Marker', 'o', 'MarkerSize', meanMarkerSize, ...
            'LineWidth', meanLineWidth, ...
            'MarkerFaceColor', 'none');

    % Hold off
    hold_off(wasHold);
end

% Change the y tick locations and labels if data is log 2 ratio
% TODO: decide_on_tick_locations.m
if isLog2Data
    % Get the default y tick locations
    yTickLocs = get(axHandle, 'YTick');

    % Keep only the locations that are integers
    % yTickValuesOriginal = 2.^(yTickLocs);
    % newYTickLocs = yTickLocs(isaninteger(yTickValuesOriginal));
    newYTickLocs = yTickLocs(isaninteger(yTickLocs));

    % Make sure 0 is included
    newYTickLocs = unique([0, newYTickLocs], 'sorted');

    % Make sure there are at least two ticks
    if numel(newYTickLocs) < 2
        yLimits = get(axHandle, 'YLim');
        maxValue = max(max(dataValues));
        minValue = min(min(dataValues));

        [~, idx] = min(abs([minValue - floor(yLimits(1)), ...
                            maxValue - ceil(yLimits(2))]));

        newYLimits = yLimits;
        if idx == 1 && floor(yLimits(1)) < 0
            newYTickLocs = [floor(yLimits(1)), newYTickLocs];
            newYLimits(1) = floor(yLimits(1));
        else
            newYTickLocs = [newYTickLocs, ceil(yLimits(2))];
            newYLimits(2) = ceil(yLimits(2));
        end
        set(axHandle, 'YLim', newYLimits);
    end

    % Get the original tick values
    newYTickValuesOriginal = 2.^(newYTickLocs);

    % Create y tick labels
    newYTickLabels = create_labels_from_numbers(newYTickValuesOriginal);

    % Set new ticks and labels
    set(axHandle, 'YTick', newYTickLocs, 'YTickMode', 'manual');
    set(axHandle, 'YTickLabel', newYTickLabels, 'YTickLabelMode', 'manual');
end

% Save figure
if ~isempty(figName)
    save_all_figtypes(figHandle, figName, figTypes);
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

% Plot a star if significant
plotStar = false;
if plotStar
    starXPos = mean(pValues);
    yLimits = get(gca, 'YLim');
    starYPos = yLimits(1) + 0.8 * (yLimits(2) - yLimits(1));
    plot(starXPos, starYPos, 'k*');
end

colorMap = [0, 0, 0];

plot(axHandle, pValues, means, 'r-o', ...
    'LineWidth', meanLineWidth, 'MarkerSize', meanMarkerSize, ...
    'MarkerFaceColor', meanColorMap);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%