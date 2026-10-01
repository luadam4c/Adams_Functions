function [handles, statsOut] = plot_grouped_jitter (data, varargin)
%% Plots a jitter plot colored by group from data (either uses plotSpread or not)
% Usage: [handles, statsOut] = plot_grouped_jitter (data, grouping (opt), condition (opt), varargin)
% Explanation:
%       This function creates a grouped jitter plot (swarm plot) where data 
%       points are spread along the x-axis to avoid overlap. It supports two 
%       levels of grouping: 
%           1. 'Condition' (mapped to x-axis positions)
%           2. 'Grouping' (mapped to color)
%       It optionally overlays descriptive statistics (mean and error bars) 
%       and performs statistical tests. Depending on the number of conditions,
%       it performs either two-sample tests between adjacent conditions or 
%       single-sample tests against a null hypothesis. Tests can be pooled 
%       across groups or evaluated separately per group. The test statistic 
%       logic adapts automatically based on the 'DataType' parameter ('normal', 
%       'proportion', or 'ratio').
%
% Example(s):
%       % Example 1: Basic comparison of two distributions
%       data = [randn(50,1); randn(50,1)+2];
%       condition = [ones(50,1); 2*ones(50,1)];
%       plot_grouped_jitter(data, [], condition);
%
%       % Example 2: Grouped data (2 groups across 2 conditions)
%       % Data: G1-C1, G1-C2, G2-C1, G2-C2
%       data = [randn(50,1); randn(50,1)+2; randn(50,1)+1; randn(50,1)+3];
%       grouping = [ones(100,1); 2*ones(100,1)];       % Maps to Color
%       condition = repmat([ones(50,1); 2*ones(50,1)], 2, 1); % Maps to X-axis
%       plot_grouped_jitter(data, grouping, condition, ...
%           'GroupingLabels', {'Control', 'Treatment'}, ...
%           'XTickLabels', {'Baseline', 'Post-Stim'}, ...
%           'StatsMode', 'withinGroupAndCond');
%
%       % Example 3: Customizing limits and statistics
%       [handles, statsOut] = plot_grouped_jitter(data, grouping, condition, ...
%           'YLimits', [-5, 10], ...
%           'RunTTest', true, ...
%           'RunRankTest', false, ...
%           'StatsMode', 'pooled', ...
%           'PlotMeanValues', true);
%
%       % Example 4
%       data = [-1; -3; 10; 49; 71; 67]; grouping = {'p1'; 'p2'; 'p3'; 'p4'; 'p5'; 'p6'}; condition = {'second'; 'second'; 'second'; 'first'; 'first'; 'first'};
%       plot_grouped_jitter(data, grouping, condition, 'StatsMode', 'pooled');
%       plot_grouped_jitter(data, grouping, condition, 'StatsMode', 'pooled', 'UsePlotSpread', false)
%
%       % Example 5
%       randVec1 = randi(10, 10, 1);
%       randVec2 = randi(10, 10, 1) + 10;
%       data = [randVec1, randVec2];
%       plot_grouped_jitter(data)
%       plot_grouped_jitter(data, 'UsePlotSpread', true)
%       plot_grouped_jitter(data, 'UsePlotSpread', false)
%       plot_grouped_jitter(data, 'UsePlotSpread', false, 'JitterWidth', 0.5)
%
%       % Example 6
%       data = [randn(50,1);randn(50,1)+3.5]*[1 1];
%       grouping = [[ones(50,1);zeros(50,1)],[randi([0,1],[100,1])]];
%       plot_grouped_jitter(data, grouping, 'RunTTest', true)
%       plot_grouped_jitter(data, grouping, 'UsePlotSpread', false, 'RunRankTest', true)
%
%       % Example 7
%       data = [randn(50,1);randn(50,1)+5;randn(50,1)+10;randn(50,1)+15];
%       grouping = [ones(50,1);zeros(50,1);ones(50,1);zeros(50,1)];
%       condition = [2*ones(50,1);2*ones(50,1);zeros(50,1);zeros(50,1)];
%       plot_grouped_jitter(data, grouping, condition)
%       plot_grouped_jitter(data, grouping, condition, 'UsePlotSpread', false)
%       plot_grouped_jitter(data, grouping, condition, 'XTickLocs', 'suppress')
%       plot_grouped_jitter(data, grouping, condition, 'UsePlotSpread', false, 'XTickLocs', 'suppress')
%
% Outputs:
%       handles     - handles to plotted objects
%                       .ax - handle to the axes
%                       .distributions - handles to the distributions
%                       .stats - handles to the stats from plotSpread
%                       .eqText - handle to the mixed effect equation text
%                       .pText - handles to the p-value text objects
%                       .sigMarker - handles to the significance markers
%                   specified as a structure
%       statsOut    - computed statistical values (means, errors, lower95s, upper95s, p-values)
%                     for mixed effect models, also includes the model object
%                     and fitted equation string
%                   specified as a structure
%
% Arguments:
%       data        - cell array of distributions or an nDatapoints-by-mDistributions array, 
%                       or an array with data that is indexed by either
%                       condition (x-axis) or grouping (color), or both.
%                   Note: The dimension with fewer elements is taken as 
%                           the parameter
%                   must be a table or a numeric array
%                       or a cell array of numeric vectors
%       grouping    - (opt) group assignment for each data point
%                   must be an array of one the following types:
%                       'cell', 'string', numeric', 'logical', 
%                           'datetime', 'duration'
%                   default == the column number for a 2D array
%       condition - (opt) condition assignment for each data point
%                   must be an array of one the following types:
%                       'cell', 'string', numeric', 'logical', 
%                           'datetime', 'duration'
%                   default == the column number for a 2D array
%       varargin    - 'DataType': type of data to calculate statistics for
%                   must be an unambiguous, case-insensitive match to one of:
%                       'normal'        - expect normal distribution, 
%                                           compute standard error of the mean
%                                           perform T-tests or ranksum/signed rank tests
%                       'proportion'    - expect Bernoulli trials,
%                                           compute standard error of a proportion
%                                           perform z-tests or exact binomial/Fisher tests
%                       'ratio'         - expect ratios
%                                           compute tests in log-scale
%                                           perform log T-tests or rank tests against 1
%                   default == 'normal'
%                   - 'StatsMode': how to compute statistics for error bars and p-values
%                   must be an unambiguous, case-insensitive match to one of:
%                       'withinGroup'        - tests run for each group separately across conditions
%                       'withinGroupAndCond' - tests run for each group separately across conditions, AND across groups within each condition
%                       'pooled'             - groups are pooled together for one general test
%                       'mixedEffect'        - mixed effect model accounting for correlated 
%                                              observations within each group
%                   default == 'withinGroup'
%                   - 'ErrorBarType': type of statistic to calculate for error bars
%                   must be an unambiguous, case-insensitive match to one of:
%                       'err95'         - 95% confidence interval
%                       'stderr'        - standard error
%                   default == 'err95'
%                   - 'GroupingOrder': how to order the unique grouping labels internally
%                   must be an unambiguous, case-insensitive match to one of:
%                       'original'      - keeps the original data grouping order
%                       'bycondition'   - orders the groups strictly by their primary condition
%                   default == 'original'
%                   - 'Weights': optional weight vector for calculating 
%                       weighted standard errors
%                   must be empty or a numeric array of the same size as data
%                   default == []
%                   - 'Means': optional means to plot instead of calculating
%                   must be a numeric array 
%                   default == []
%                   - 'Errors': optional symmetric errors to plot instead of calculating
%                   must be a numeric array 
%                   default == []
%                   - 'Upper95s': optional upper bounds for asymmetric error bars
%                   must be a numeric array
%                   default == []
%                   - 'Lower95s': optional lower bounds for asymmetric error bars
%                   must be a numeric array
%                   default == []
%                   - 'UsePlotSpread': whether to use plotSpread.m
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == [] (true if plotSpread.m is found)
%                   - 'JitterWidth': width of the jitter
%                   must be a non-negative scalar
%                   default == [] (plotSpread default) or 0.3 (manual plot)
%                   - 'OffsetGroups': whether to offset groups along the x-axis
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true if StatsMode is 'withinGroup' or 'withinGroupAndCond', false otherwise
%                   - 'PlotMeanValues': whether to plot the mean values
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotErrorBars': whether to plot error bars 
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'RunTTest': whether to run unpaired t-test
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'RunRankTest': whether to run unpaired 
%                                       Wilcoxon rank-sum test
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == [0.5, nConditions + 0.5]
%                   - 'YLimits': limits of y axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == [] (auto-scale)
%                   - 'XTickLocs': x tick locations
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == unique(condition)
%                   - 'XTickLabels': x tick labels in place of parameter values
%                   must be a cell array of character vectors/strings
%                   default == {}
%                   - 'GroupingLabels': labels for the groupings, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == {'Group #1', 'Group #2', ...}
%                   - 'XTickAngle': angle for parameter tick labels
%                   must be a numeric scalar
%                   default == 0
%                   - 'XLabel': label for the x axis, 
%                   must be a string scalar or a character vector 
%                   default == none
%                   - 'YLabel': label for the y axis, 
%                   must be a string scalar or a character vector 
%                   default == none
%                   - 'Title': title for the plot, 
%                   must be a string scalar or a character vector 
%                   default == none
%                   - 'Marker': marker for the individual data points
%                   must be a string scalar or a character vector
%                   default == '.'
%                   - 'LineWidth': line width for the individual data points when marker is 'o' or 'x'
%                   must be a positive numeric scalar
%                   default == 1
%                   - 'MeanMarker': marker for the mean value points
%                   must be a string scalar or a character vector
%                   default == 'o'
%                   - 'MeanMarkerSize': size of the mean marker
%                   must be a positive numeric scalar
%                   default == 10
%                   - 'MeanLineWidth': line width of error bars and mean points
%                   must be a positive numeric scalar
%                   default == 1.5
%                   - 'ErrorBarCapSize': size of the error bar caps
%                   must be a non-negative numeric scalar
%                   default == 0
%                   - 'SigLevel': significance level for hypothesis tests
%                   must be a numeric scalar strictly between 0 and 1
%                   default == 0.05
%                   - 'ColorMap': a color map for each group
%                   must be a numeric array with 3 columns
%                   default == set in decide_on_colormap.m
%                   - 'LegendLocation': location for legend
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'auto'      - use default
%                       'suppress'  - no legend
%                       anything else recognized by the legend() function
%                   default == 'suppress' if nGroups == 1 
%                               'best' if nGroups is 2~9
%                               'eastoutside' if nGroups is 10+
%                   - 'AxesHandle': axes handle to plot on
%                   must be a empty or an axes object handle
%                   default == set in set_axes_properties.m
%                   - Any other parameter-value pair for plotSpread() or plot()
%
% Requires:
%       cd/addpath_custom.m
%       cd/compute_stats.m
%       cd/create_default_grouping.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/decide_on_colormap.m
%       cd/hold_off.m
%       cd/hold_on.m
%       cd/locate_functionsdir.m
%       cd/plot_test_result.m
%       cd/plot_text.m
%       cd/set_axes_properties.m
%       cd/struct2arglist.m
%       cd/test_normality.m
%       ~/Downloaded_Functions/plotSpread/plotSpread.m
%
% Used by:
%       cd/m3ha_plot_figure08.m
%       cd/m3ha_simulate_population.m
%       cd/virt_plot_jitter.m

% File History:
% 2020-04-13 Modified from plot_violin.m
% 2020-04-18 Now uses the categoryIdx option in plotSpread
% 2025-08-28 Now uses the distributionIdx option in plotSpread
% 2025-08-28 Added 'UsePlotSpread' as an optional argument
% 2025-08-28 Implemented stats plotting by Gemini
% 2025-08-29 Implemented normality testing from plot_tuning_curve.m by Gemini
% 2025-09-16 Now outputs distributions in handles for manual case
% 2025-09-17 Made 'JitterWidth' an optional argument by Gemini
% 2025-09-17 Added 'AxesHandle' as an optional argument
% 2026-01-23 Added 'StatisticsByGroup' optional argument by Gemini
% 2026-01-23 Added 'YLimits' optional argument and improved docs by Gemini
% 2026-03-04 Fixed behavior when 'XTickLocs' is set to 'suppress'
% 2026-03-04 Made marker, meanMarker, meanMarkerSize, meanLineWidth, errorBarCapSize, and sigLevel optional arguments
% 2026-03-04 Added 'LineWidth' as an optional argument
% 2026-03-12 Implemented single condition statistical testing against null hypothesis based on data type and normality
% 2026-03-12 Added support for 'ratio' Data Types and fixed corresponding logic for 'proportion' tests
% 2026-03-12 Handled iterating condition pairs sequentially when nConditions > 2 via arrayfun and subfunction packaging
% 2026-03-12 Changed default legend location to 'best' if nGroups is 2~9 
% 2026-03-17 Added 'ErrorBarType' argument and updated to use compute_stats.m
% 2026-03-18 Added 'GroupingOrder' argument to sort by condition by Gemini
% 2026-03-20 Changed 'StatisticsByGroup' to 'StatsMode' and implemented global mixed-effects model by Gemini
% 2026-03-20 Fixed unrecognized field name bug, dynamically sets p-value colors, and prevents text overlap by Gemini
% 2026-03-20 Now allows asymmetric error bars by Gemini
% 2026-03-20 Updated p-values to use normalized units so they update with y-axis limits, and removed pooledStatColor to rely on default test colors
% 2026-03-22 Aligned group points and means, separated by jitter width, and added comparison lines by Gemini
% 2026-03-22 Removed set_graphics_color and delegated color/position logic to plot_test_result.m
% 2026-03-22 Now outputs plot_test_result.m and plot_text.m graphic handles directly
% 2026-03-24 Added 'OffsetGroups' optional argument by Gemini
% 2026-09-30 Centered manual jitter and forced padding for single conditions to prevent off-center plots

%% Add withinGroupAndCond to valid stats modes
%% Define hard-coded parameters for the function
validDataTypes = {'normal', 'proportion', 'ratio'};
validErrorBarTypes = {'err95', 'stderr'};
validGroupingOrders = {'original', 'bycondition'};
validStatsModes = {'withinGroup', 'pooled', 'mixedEffect', 'withinGroupAndCond'};
defaultJitterWidthNotPlotSpread = 0.3;

%% TODO: Make the following optional arguments
maxNGroupsForInnerLegend = 8;
maxNGroupsForOuterLegends = 25;
forceVectorInput = true;
xTickLimitPadding = 0.5;
pooledStatColor = 'k';
statTextYPosStepSize = 0.1;
statStarYLocOffset = 0.05;

%% Set default values for optional arguments
groupingDefault = [];           % set later
distributionDefault = [];       % set later
dataTypeDefault = 'normal';
errorBarTypeDefault = 'err95';
groupingOrderDefault = 'original';
weightsDefault = [];
meansDefault = [];
errorsDefault = [];
upper95sDefault = [];
lower95sDefault = [];
usePlotSpreadDefault = [];      % set later
jitterWidthDefault = [];        % set later
offsetGroupsDefault = [];       % set later based on statsMode
plotMeanValuesDefault = true;
plotErrorBarsDefault = true;
statsModeDefault = 'withinGroup';
runTTestDefault = true;
runRankTestDefault = true;
xLimitsDefault = [];            % set later
yLimitsDefault = [];            % set later
xTickLocsDefault = [];           % set later
xTickLabelsDefault = {};        % set later
groupingLabelsDefault = '';     % set later
xTickAngleDefault = [];         % set later
xLabelDefault = '';             % no x label by default
yLabelDefault = '';             % no y label by default
titleDefault = '';              % no title by default
markerDefault = '.';
lineWidthDefault = 1;
meanMarkerDefault = 'o';
meanMarkerSizeDefault = 10;
meanLineWidthDefault = 1.5;
errorBarCapSizeDefault = 0;
sigLevelDefault = 0.05;
colorMapDefault = [];           % set later
legendLocationDefault = 'auto'; % set later
axHandleDefault = [];           % axHandle by default

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

% Add optional inputs to the Input Parser
addOptional(iP, 'grouping', groupingDefault, ...
    @(x) validateattributes(x, {'cell', 'string', 'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addOptional(iP, 'condition', distributionDefault, ...
    @(x) validateattributes(x, {'cell', 'string', 'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'DataType', dataTypeDefault, ...
    @(x) any(validatestring(x, validDataTypes)));
addParameter(iP, 'ErrorBarType', errorBarTypeDefault, ...
    @(x) any(validatestring(x, validErrorBarTypes)));
addParameter(iP, 'GroupingOrder', groupingOrderDefault, ...
    @(x) any(validatestring(x, validGroupingOrders)));
addParameter(iP, 'Weights', weightsDefault, ...
    @(x) assert(isempty(x) || isnumeric(x), 'Weights must be empty or a numeric array!'));
addParameter(iP, 'Means', meansDefault, ...
    @(x) isempty(x) || isnumeric(x));
addParameter(iP, 'Errors', errorsDefault, ...
    @(x) isempty(x) || isnumeric(x));
addParameter(iP, 'Upper95s', upper95sDefault, ...
    @(x) isempty(x) || isnumeric(x));
addParameter(iP, 'Lower95s', lower95sDefault, ...
    @(x) isempty(x) || isnumeric(x));
addParameter(iP, 'UsePlotSpread', usePlotSpreadDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'JitterWidth', jitterWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'OffsetGroups', offsetGroupsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotMeanValues', plotMeanValuesDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotErrorBars', plotErrorBarsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'StatsMode', statsModeDefault, ...
    @(x) any(validatestring(x, validStatsModes)));
addParameter(iP, 'RunTTest', runTTestDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RunRankTest', runRankTestDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'XTickLocs', xTickLocsDefault, ...
    @(x) isempty(x) || ischar(x) && (strcmpi(x, 'suppress') || strcmpi(x, 'suppressed')) || ...
        isnumeric(x) && isvector(x));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'XTickLabels', xTickLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'GroupingLabels', groupingLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'XTickAngle', xTickAngleDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'XLabel', xLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Title', titleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Marker', markerDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'LineWidth', lineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'MeanMarker', meanMarkerDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'MeanMarkerSize', meanMarkerSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'MeanLineWidth', meanLineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'ErrorBarCapSize', errorBarCapSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'SigLevel', sigLevelDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>', 0, '<', 1}));
addParameter(iP, 'ColorMap', colorMapDefault);
addParameter(iP, 'LegendLocation', legendLocationDefault, ...
    @(x) all(islegendlocation(x, 'ValidateMode', true)));
addParameter(iP, 'AxesHandle', axHandleDefault);

% Read parsed results from the Input Parser
parse(iP, data, varargin{:});
grouping = iP.Results.grouping;
condition = iP.Results.condition;
dataType = validatestring(iP.Results.DataType, validDataTypes);
errorBarType = validatestring(iP.Results.ErrorBarType, validErrorBarTypes);
groupingOrder = validatestring(iP.Results.GroupingOrder, validGroupingOrders);
weights = iP.Results.Weights;
means = iP.Results.Means;
errors = iP.Results.Errors;
upper95s = iP.Results.Upper95s;
lower95s = iP.Results.Lower95s;
usePlotSpread = iP.Results.UsePlotSpread;
jitterWidth = iP.Results.JitterWidth;
offsetGroups = iP.Results.OffsetGroups;
plotMeanValues = iP.Results.PlotMeanValues;
plotErrorBars = iP.Results.PlotErrorBars;
statsMode = validatestring(iP.Results.StatsMode, validStatsModes);
runTTest = iP.Results.RunTTest;
runRankTest = iP.Results.RunRankTest;
xTickLocs = iP.Results.XTickLocs;
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;
xTickLabels = iP.Results.XTickLabels;
groupingLabels = iP.Results.GroupingLabels;
xTickAngle = iP.Results.XTickAngle;
xLabel = iP.Results.XLabel;
yLabel = iP.Results.YLabel;
plotTitle = iP.Results.Title;
marker = iP.Results.Marker;
lineWidth = iP.Results.LineWidth;
meanMarker = iP.Results.MeanMarker;
meanMarkerSize = iP.Results.MeanMarkerSize;
meanLineWidth = iP.Results.MeanLineWidth;
errorBarCapSize = iP.Results.ErrorBarCapSize;
sigLevel = iP.Results.SigLevel;
colorMap = iP.Results.ColorMap;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);
axHandle = iP.Results.AxesHandle;

% Keep unmatched arguments for the plotSpread() function
otherArguments = struct2arglist(iP.Unmatched);

% Ensure weights match the dimensions of data if provided
if ~isempty(weights)
    assert(isequal(size(weights), size(data)), 'Weights must be the same dimension as data!');
end

%% Preparation
% Decide whether to use the plotSpread tool based on its availability
if isempty(usePlotSpread)
    if exist('plotSpread.m', 'file') == 2
        usePlotSpread = true;
    else
        usePlotSpread = false;
    end
end

% If not compiled, add directories to the search path for required functions
if usePlotSpread && exist('plotSpread.m', 'file') ~= 2 && ~isdeployed
    try
        % Locate the functions directory
        functionsDirectory = locate_functionsdir;

        % Add path for plotSpread()
        addpath_custom(fullfile(functionsDirectory, ...
                                'Downloaded_Functions', 'plotSpread'));
    catch ME
        disp('An error occurred when looking for plotSpread.m:');
        disp(ME.message);
        disp('plotSpread.m will not be used!');
    end
end

% Decide on jitter width if not provided by user
if isempty(jitterWidth)
    if usePlotSpread
        % Let plotSpread.m use its default
        jitterWidth = [];
    else
        % Set default for manual plot
        jitterWidth = defaultJitterWidthNotPlotSpread;
    end
end

% Decide on the grouping vector and labels if data is a matrix or cell array
% Note: If no grouping vector, each column is a group
[grouping, uniqueGroupValues, groupingLabels, data] = ...
    create_default_grouping('Stats', data, 'Grouping', grouping, ...
                            'GroupingLabels', groupingLabels, ...
                            'GroupingLabelPrefix', 'Group');

% Decide on the condition vector and labels if data is a matrix or cell array
% Note: If no condition vector, each column is a condition
[condition, uniqueConditionValues, xTickLabels, data] = ...
    create_default_grouping('Stats', data, 'Grouping', condition, ...
                            'GroupingLabels', xTickLabels, ...
                            'GroupingLabelPrefix', 'Condition');

% Concatenate everything into a single column vector if required
if forceVectorInput
    if ~isempty(weights)
        % Force non-vectors as cell arrays of numeric vectors
        [data, grouping, condition, weights] = ...
            argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', false), ...
                    data, grouping, condition, weights);
        
        % If data and grouping and condition and weights are cell arrays of numeric vectors, pool them
        if iscellnumeric(data) && iscellnumeric(grouping) && iscellnumeric(condition) && iscellnumeric(weights)
            [data, grouping, condition, weights] = ...
                argfun(@(x) vertcat(x{:}), data, grouping, condition, weights);
        end
    else
        % Force non-vectors as cell arrays of numeric vectors
        [data, grouping, condition] = ...
            argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', false), ...
                    data, grouping, condition);
        
        % If data and grouping and condition are cell arrays of numeric vectors, pool them
        if iscellnumeric(data) && iscellnumeric(grouping) && iscellnumeric(condition)
            [data, grouping, condition] = ...
                argfun(@(x) vertcat(x{:}), data, grouping, condition);
        end
    end
end

% Count the number of unique conditions
nConditions = numel(uniqueConditionValues);

% Count the number of unique groups
nGroups = numel(uniqueGroupValues);

% Count the total number of data points
nPoints = numel(data);

% Determine maximum groups per condition to see if shifting is necessary
maxGroupsPerCond = 1;
for iCond = 1:nConditions
    condVal = uniqueConditionValues(iCond);
    numGrp = numel(unique(grouping(condition == condVal)));
    if numGrp > maxGroupsPerCond
        maxGroupsPerCond = numGrp;
    end
end

% Decide whether to offset groups based on statsMode if not explicitly provided
if isempty(offsetGroups)
    if startsWith(statsMode, 'withinGroup', 'IgnoreCase', true)
        offsetGroups = true;
    else
        offsetGroups = false;
    end
end

% Determine group offsets for x-axis to align clusters and prevent overlap based on jitter width
groupOffsets = zeros(nGroups, 1);
if maxGroupsPerCond > 1 && offsetGroups
    % Fallback to default jitter width for spacing if jitterWidth is empty (e.g. plotSpread defaults)
    if isempty(jitterWidth)
        actualJitterWidth = defaultJitterWidthNotPlotSpread;
    else
        actualJitterWidth = jitterWidth;
    end
    groupSpacing = actualJitterWidth * 1.1;
    
    % Calculate total desired spread
    totalSpread = groupSpacing * (nGroups - 1);
    
    % Capping the total spread to 0.8 to strictly prevent overlapping into adjacent condition spaces
    if totalSpread > 0.8
        totalSpread = 0.8;
    end
    
    groupOffsets = linspace(-totalSpread/2, totalSpread/2, nGroups)';
end

% Create a shifted condition array to separate groups horizontally
conditionShifted = zeros(size(condition));
for iGroup = 1:nGroups
    isCurrentGroup = (grouping == uniqueGroupValues(iGroup));
    conditionShifted(isCurrentGroup) = condition(isCurrentGroup) + groupOffsets(iGroup);
end

% Decide whether to update the x tick locations
if ischar(xTickLocs) && (strcmpi(xTickLocs, 'suppress') || strcmpi(xTickLocs, 'suppressed'))
    toUpdateXTicks = false;
else
    toUpdateXTicks = true;
end

% Decide on the x tick locations
if isempty(xTickLocs) || ~toUpdateXTicks
    % Define x tick locations based on the unique condition indices
    xTickLocs = uniqueConditionValues;
end

% Decide on the color map, using the lines map by default
if isempty(colorMap)
    colorMap = @lines;
end
colorMap = decide_on_colormap(colorMap, nGroups, 'ForceCellOutput', true);

% Sort groups internally by condition if requested
if strcmpi(groupingOrder, 'bycondition')
    % Associate each group with its primary condition
    groupCondIdx = zeros(nGroups, 1);
    for iGroup = 1:nGroups
        % Find conditions where this group appears
        condsForGroup = condition(grouping == uniqueGroupValues(iGroup));
        if ~isempty(condsForGroup)
            % Find the index of the first condition natively linked to this group
            groupCondIdx(iGroup) = find(uniqueConditionValues == condsForGroup(1), 1);
        else
            % Assign infinity if no condition is logically associated
            groupCondIdx(iGroup) = Inf;
        end
    end
    
    % Sort by condition index while preserving original order for ties
    [~, sortIdx] = sort(groupCondIdx, 'ascend');
    
    % Reorder the unique group values, labels, and colormap natively
    sortedUniqueGroupValues = uniqueGroupValues(sortIdx);
    if ~isempty(groupingLabels)
        groupingLabels = groupingLabels(sortIdx);
    end
    colorMap = colorMap(sortIdx);
    
    % Remap the main grouping numeric values to strictly reflect the new ordering internally
    newGrouping = zeros(size(grouping));
    for iGroup = 1:nGroups
        oldVal = sortedUniqueGroupValues(iGroup);
        newGrouping(grouping == oldVal) = iGroup;
    end
    grouping = newGrouping;
    uniqueGroupValues = (1:nGroups)';
end

% Set the legend location by utilizing the dedicated subfunction if set to 'auto'
if strcmpi(legendLocation, 'auto')
    legendLocation = set_default_legend_location(nGroups, maxNGroupsForInnerLegend, maxNGroupsForOuterLegends);
end

% Decide on the axes to plot on
axHandle = set_axes_properties('AxesHandle', axHandle);

%% Do the job
% Return immediately if there is no data to plot
if isempty(data)
    handles = struct;
    statsOut = struct;
    return
end

% Retain the current plot on the axes
wasHold = hold_on(axHandle);

% Plot the data points using plotSpread or manually based on preference
if usePlotSpread
    % Prevent plotSpread from plotting its own means if we are plotting custom ones
    if plotMeanValues || plotErrorBars
        % Force otherArguments to be a row cell array to prevent dimension mismatch during concatenation
        otherArguments = [otherArguments(:)', {'showMM', 0}];
    end

    % Call plotSpread to plot the swarm plot using the shifted condition array
    output = plotSpread(axHandle, data, 'distributionIdx', conditionShifted, ...
                            'categoryIdx', grouping, ...
                            'categoryLabels', groupingLabels, ...
                            'categoryColors', colorMap, ...
                            'spreadWidth', jitterWidth, ...
                            otherArguments{:});

    % Reformat and store handles from plotSpread output
    distributions = output{1};
    stats = output{2};
    ax = output{3};
    handles.distributions = distributions;
    handles.stats = stats;
    handles.ax = ax;
    
    % Initialize output graphic handles tracking lists
    handles.eqText = gobjects(0);
    handles.pText = gobjects(0);
    handles.sigMarker = gobjects(0);
else
    % Give each shifted condition index some centered jitter for plotting
    jitterOffset = rand(nPoints, 1) - 0.5;
    if nPoints > 1
        jitterOffset = jitterOffset - mean(jitterOffset);
    else
        jitterOffset = 0;
    end
    conditionWithJitter = conditionShifted + (jitterWidth * jitterOffset);

    % Pre-allocate a graphics object array for plot handles
    distributions = gobjects(nGroups, 1);
    
    % Plot each group manually with the specified marker and color map
    for iGroup = 1:nGroups
        % Get the value for the current group
        currentGroupValue = uniqueGroupValues(iGroup);

        % Find indices for data points belonging to this group
        isCurrentGroup = (grouping == currentGroupValue);

        % Extract the data for this group
        xCoords = conditionWithJitter(isCurrentGroup);
        yCoords = data(isCurrentGroup);

        % Get the color and label for this group
        groupColor = colorMap{iGroup};
        groupLabel = groupingLabels{iGroup};

        % Plot this group's data, ensuring no lines connect markers
        distributions(iGroup) = ...
            plot(axHandle, xCoords, yCoords, 'Marker', marker, ...
                    'Color', groupColor, 'DisplayName', groupLabel, ...
                    'LineStyle', 'none', 'LineWidth', lineWidth, ...
                    otherArguments{:});
    end

    % Export handles to the output structure
    handles.ax = axHandle;
    handles.distributions = distributions;
    handles.stats = [];
    
    % Initialize output graphic handles tracking lists
    handles.eqText = gobjects(0);
    handles.pText = gobjects(0);
    handles.sigMarker = gobjects(0);
end

%% Finalize main plot
% Update x tick locations if desired
if toUpdateXTicks
    xticks(axHandle, xTickLocs);
end

% Update x tick labels if desired
if toUpdateXTicks && ~isempty(xTickLabels)
    xticklabels(axHandle, xTickLabels);
end

% Modify x tick angle
if ~isempty(xTickAngle)
    xtickangle(axHandle, xTickAngle);
end

% Decide on x axis limits based on x tick locations
if isempty(xLimits)
    xLimits = [min(xTickLocs) - xTickLimitPadding, max(xTickLocs) + xTickLimitPadding];
end

% Modify x limits, overriding suppression for single conditions to maintain centering
if ~(ischar(xLimits) && strcmpi(xLimits, 'suppress'))
    xlim(axHandle, xLimits);
elseif ischar(xLimits) && strcmpi(xLimits, 'suppress') && numel(xTickLocs) == 1
    xlim(axHandle, [xTickLocs(1) - xTickLimitPadding, xTickLocs(1) + xTickLimitPadding]);
end

% Modify y limits
if ~(ischar(yLimits) && strcmpi(yLimits, 'suppress')) && ~isempty(yLimits)
    ylim(axHandle, yLimits);
end

% Set x label
if ~isempty(xLabel)
    xlabel(axHandle, xLabel);
end

% Set y label
if ~isempty(yLabel)
    ylabel(axHandle, yLabel);
end

% Set title
if ~isempty(plotTitle)
    title(axHandle, plotTitle);
end

% Generate a legend if there is more than one trace
if ~strcmpi(legendLocation, 'suppress')
    legend(axHandle, 'location', legendLocation);
end

%% Initialize statistics output structure
%% Check for withinGroup variations to allocate accordingly
statsOut = struct();
if startsWith(statsMode, 'withinGroup', 'IgnoreCase', true)
    statsOut.means = NaN(nGroups, nConditions);
    statsOut.errors = NaN(nGroups, nConditions);
    statsOut.lower95s = NaN(nGroups, nConditions);
    statsOut.upper95s = NaN(nGroups, nConditions);
else
    statsOut.means = NaN(1, nConditions);
    statsOut.errors = NaN(1, nConditions);
    statsOut.lower95s = NaN(1, nConditions);
    statsOut.upper95s = NaN(1, nConditions);
end
statsOut.pValueT = [];
statsOut.pValueR = [];

%% Pre-compute Global Mixed Effect Model
% If requested, fit the model once to compute means, asymmetric bounds, and p-values
mixedEffectSuccess = false;
mixedEffectMeans = [];
mixedEffectLower = [];
mixedEffectUpper = [];
mixedEffectPValues = [];
mixedEffectEquation = '';

if strcmpi(statsMode, 'mixedEffect')
    [mixedEffectMeans, mixedEffectLower, mixedEffectUpper, mixedEffectPValues, mixedEffectSuccess, mixedEffectModel, mixedEffectEquation] = ...
        compute_global_mixed_effects(data, condition, grouping, dataType, errorBarType, uniqueConditionValues, (runTTest || runRankTest));
    
    % Store mixed effect results in the output stats structure
    statsOut.mixedEffectMeans = mixedEffectMeans;
    statsOut.mixedEffectLower = mixedEffectLower;
    statsOut.mixedEffectUpper = mixedEffectUpper;
    statsOut.mixedEffectPValues = mixedEffectPValues;
    statsOut.mixedEffectSuccess = mixedEffectSuccess;
    statsOut.mixedEffectModel = mixedEffectModel;
    statsOut.mixedEffectEquation = mixedEffectEquation;

    % If the model failed to fit, gracefully fallback to pooled stats
    if ~mixedEffectSuccess
        warning('Global mixed effect model failed to converge. Falling back to pooled statistics.');
        statsMode = 'pooled';
    end
end

%% Plot statistics
%% Adapt checking logic to encompass the combined within group/cond mode
% Plot means and error bars for each condition and store them regardless of plotting if an output arg is expected
if plotMeanValues || plotErrorBars || nargout > 1
    % Loop through each condition to plot statistics
    for iCond = 1:nConditions
        % Find data and weights for the current condition
        currentCondValue = uniqueConditionValues(iCond);
        isCurrentCond = (condition == currentCondValue);
        dataThisCond = data(isCurrentCond);
        
        % Extract weights for this condition if provided
        if ~isempty(weights)
            weightsThisCond = weights(isCurrentCond);
        else
            weightsThisCond = [];
        end
        
        % Check if statistics should be plotted by group
        if startsWith(statsMode, 'withinGroup', 'IgnoreCase', true)
            % Separate by group for the current condition
            groupingThisCond = grouping(isCurrentCond);
            groupsInThisCond = unique(groupingThisCond);
            nGroupsInThisCond = numel(groupsInThisCond);
            
            % Plot statistics if there are groups present
            if nGroupsInThisCond > 0
                % Split data by individual groups
                dataByGroup = arrayfun(@(x) dataThisCond(groupingThisCond == x), ...
                                        groupsInThisCond, 'UniformOutput', false);
                
                % Split weights by individual groups if provided
                if ~isempty(weightsThisCond)
                    weightsByGroup = arrayfun(@(x) weightsThisCond(groupingThisCond == x), ...
                                            groupsInThisCond, 'UniformOutput', false);
                else
                    weightsByGroup = cell(size(groupsInThisCond));
                end

                % Calculate the x-axis positions for the group means based on computed group offsets
                xMeanPositions = zeros(nGroupsInThisCond, 1);
                for iGroup = 1:nGroupsInThisCond
                    actualGroupIdx = find(uniqueGroupValues == groupsInThisCond(iGroup));
                    xMeanPositions(iGroup) = currentCondValue + groupOffsets(actualGroupIdx);
                end

                % Loop through and plot statistics for each group
                for iGroup = 1:nGroupsInThisCond
                    % Extract the data and weights for the current group
                    groupData = dataByGroup{iGroup};
                    groupWeights = weightsByGroup{iGroup};
                    
                    % Find the true index of the current group
                    actualGroupIdx = find(uniqueGroupValues == groupsInThisCond(iGroup));

                    % Check available custom bounding parameters
                    hasCustomAsymmetric = ~isempty(means) && ~isempty(lower95s) && ~isempty(upper95s);
                    hasCustomSymmetric = ~isempty(means) && ~isempty(errors);

                    % Calculate or extract statistics for the group
                    if hasCustomAsymmetric
                        groupStat = means(actualGroupIdx, iCond);
                        groupLower = lower95s(actualGroupIdx, iCond);
                        groupUpper = upper95s(actualGroupIdx, iCond);
                        groupError = (groupUpper - groupLower) / 2; % Fallback for symmetric record
                    elseif hasCustomSymmetric
                        groupStat = means(actualGroupIdx, iCond);
                        groupError = errors(actualGroupIdx, iCond);
                        groupLower = groupStat - groupError;
                        groupUpper = groupStat + groupError;
                    else
                        % Ensure dimension 1 is explicitly passed when providing optional args to compute_stats
                        if strcmpi(errorBarType, 'err95')
                            [groupStat, groupLower, groupUpper] = ...
                                compute_stats(groupData, {'mean', 'lower95', 'upper95'}, 1, ...
                                        'DataType', dataType, 'Weights', groupWeights, 'IgnoreNan', true);
                            groupError = (groupUpper - groupLower) / 2;
                        else
                            [groupStat, groupError] = ...
                                compute_stats(groupData, {'mean', 'stderr'}, 1, ...
                                        'DataType', dataType, 'Weights', groupWeights, 'IgnoreNan', true);
                            groupLower = groupStat - groupError;
                            groupUpper = groupStat + groupError;
                        end
                    end
                    
                    % Store computed statistics for the group
                    statsOut.means(actualGroupIdx, iCond) = groupStat;
                    statsOut.errors(actualGroupIdx, iCond) = groupError;
                    statsOut.lower95s(actualGroupIdx, iCond) = groupLower;
                    statsOut.upper95s(actualGroupIdx, iCond) = groupUpper;

                    % Retrieve the color for the current group
                    groupColor = colorMap{actualGroupIdx};

                    % Plot error bars
                    if plotErrorBars
                        negErr = groupStat - groupLower;
                        posErr = groupUpper - groupStat;
                        errorbar(xMeanPositions(iGroup), groupStat, negErr, posErr, ...
                                 'Color', groupColor, 'LineWidth', meanLineWidth, ...
                                 'CapSize', errorBarCapSize, 'HandleVisibility', 'off');
                    end
                    
                    % Plot mean values
                    if plotMeanValues
                        plot(xMeanPositions(iGroup), groupStat, meanMarker, ...
                             'MarkerEdgeColor', groupColor, ...
                             'MarkerSize', meanMarkerSize, 'LineWidth', meanLineWidth, ...
                             'HandleVisibility', 'off');
                    end
                end
            end
        else
            % Process and plot statistics pooled across groups or using precomputed mixed effects
            if ~isempty(dataThisCond)

                % Check available custom bounding parameters
                hasCustomAsymmetric = ~isempty(means) && ~isempty(lower95s) && ~isempty(upper95s);
                hasCustomSymmetric = ~isempty(means) && ~isempty(errors);

                % Calculate or extract statistics for the pooled condition
                if hasCustomAsymmetric
                    pooledStat = means(iCond);
                    pooledLower = lower95s(iCond);
                    pooledUpper = upper95s(iCond);
                    pooledError = (pooledUpper - pooledLower) / 2;
                elseif hasCustomSymmetric
                    pooledStat = means(iCond);
                    pooledError = errors(iCond);
                    pooledLower = pooledStat - pooledError;
                    pooledUpper = pooledStat + pooledError;
                elseif strcmpi(statsMode, 'mixedEffect')
                    pooledStat = mixedEffectMeans(iCond);
                    pooledLower = mixedEffectLower(iCond);
                    pooledUpper = mixedEffectUpper(iCond);
                    pooledError = (pooledUpper - pooledLower) / 2;
                else
                    % Ensure dimension 1 is explicitly passed when providing optional args to compute_stats
                    if strcmpi(errorBarType, 'err95')
                        [pooledStat, pooledLower, pooledUpper] = compute_stats(dataThisCond, {'mean', 'lower95', 'upper95'}, 1, ...
                                    'DataType', dataType, 'Weights', weightsThisCond, 'IgnoreNan', true);
                        pooledError = (pooledUpper - pooledLower) / 2;
                    else
                        [pooledStat, pooledError] = compute_stats(dataThisCond, {'mean', 'stderr'}, 1, ...
                                    'DataType', dataType, 'Weights', weightsThisCond, 'IgnoreNan', true);
                        pooledLower = pooledStat - pooledError;
                        pooledUpper = pooledStat + pooledError;
                    end
                end

                % Store computed statistics for the pooled data
                statsOut.means(1, iCond) = pooledStat;
                statsOut.errors(1, iCond) = pooledError;
                statsOut.lower95s(1, iCond) = pooledLower;
                statsOut.upper95s(1, iCond) = pooledUpper;

                % Plot error bars
                if plotErrorBars
                    negErr = pooledStat - pooledLower;
                    posErr = pooledUpper - pooledStat;
                    errorbar(currentCondValue, pooledStat, negErr, posErr, ...
                             'Color', pooledStatColor, 'LineWidth', meanLineWidth, ...
                             'CapSize', errorBarCapSize, 'HandleVisibility', 'off');
                end
                
                % Plot mean values
                if plotMeanValues
                    plot(currentCondValue, pooledStat, meanMarker, ...
                         'MarkerEdgeColor', pooledStatColor, ...
                         'MarkerSize', meanMarkerSize, 'LineWidth', meanLineWidth, ...
                         'HandleVisibility', 'off');
                end
            end
        end
    end
end

%% Evaluate tests and plot equations
%% Calculates proper offset boundaries for condition comparisons and group comparisons seamlessly
% Determine if we have a mixed effect equation to plot
hasMixedEffectEq = strcmpi(statsMode, 'mixedEffect') && mixedEffectSuccess && ~isempty(mixedEffectEquation);
eqOffset = 0;
if hasMixedEffectEq
    eqOffset = 1;
end

% Set up common check for within-group modes
isWithinGroupMode = startsWith(statsMode, 'withinGroup', 'IgnoreCase', true);

% Run statistical tests across conditions or against null hypotheses
if (runTTest || runRankTest) || hasMixedEffectEq
    % Check if each group only has a single condition
    isSingleConditionPerGroup = false;
    if isWithinGroupMode
        nConditionsPerGroup = arrayfun(@(x) numel(unique(condition(grouping == x))), uniqueGroupValues);
        if all(nConditionsPerGroup == 1)
            isSingleConditionPerGroup = true;
        end
    end

    % Expand y-axis limits on the top to show test results and equations without interfering with data
    yLims = ylim(handles.ax);
    yRange = diff(yLims);
    
    % Calculate the total number of lines to reserve at the top of the plot
    numTestsToRun = runTTest + runRankTest;
    
    % Reserve lines for between-condition tests
    linesForCondTests = 0;
    if nConditions >= 2 && ~isSingleConditionPerGroup 
        if isWithinGroupMode
            linesForCondTests = numTestsToRun * nGroups;
        else
            linesForCondTests = numTestsToRun;
        end
    elseif nConditions == 1 || isSingleConditionPerGroup
        % For 1-sample tests
        linesForCondTests = numTestsToRun;
    end
    
    % Reserve lines for between-group tests (when comparing adjacent groups within conditions)
    linesForGroupTests = 0;
    if strcmpi(statsMode, 'withinGroupAndCond') && nGroups >= 2
        linesForGroupTests = numTestsToRun * (nGroups - 1);
    end
    
    % Tally the grand total of lines to shift plot limits
    numLinesAtTop = linesForCondTests + linesForGroupTests + eqOffset;

    % Scale the y-axis if space is needed
    if numLinesAtTop > 0
        % Calculate the fraction of the plot height needed for the statistical text and stars
        fractionNeeded = numLinesAtTop * statTextYPosStepSize + statStarYLocOffset;

        % Ensure the fraction needed is physically possible (less than 100% of the plot)
        if fractionNeeded < 1
            % Calculate the new total range needed to accommodate both the data and the text
            newYRange = yRange / (1 - fractionNeeded);

            % Expand the upper y-limit by adding the new range to the existing lower limit
            yLims(2) = yLims(1) + newYRange;

            % Apply the newly calculated limits to the current axes
            ylim(handles.ax, yLims);

            % Update the stored range variable to reflect the newly expanded axes
            yRange = newYRange;
        end
    end

    % Plot the fitted equation on the plot after y limits are finalized to prevent overlap
    % and store the returned text handle
    if hasMixedEffectEq
        handles.eqText = plot_text(mixedEffectEquation, 'TextLocation', 'topleft', 'AxesHandle', handles.ax);
    end

    % Only proceed to evaluating stats if tests were requested
    if (runTTest || runRankTest)
        if nConditions >= 2 && ~isSingleConditionPerGroup
            % Get the values for the adjacent conditions to compare
            cond1Values = uniqueConditionValues(1:end-1);
            cond2Values = uniqueConditionValues(2:end);
            
            % Pack parameters for the arrayfun iteration
            testParams.data = data;
            testParams.condition = condition;
            testParams.grouping = grouping;
            testParams.uniqueGroupValues = uniqueGroupValues;
            testParams.nGroups = nGroups;
            testParams.statsMode = statsMode;
            testParams.runTTest = runTTest;
            testParams.runRankTest = runRankTest;
            testParams.dataType = dataType;
            testParams.sigLevel = sigLevel;
            testParams.yLims = yLims;
            testParams.eqOffset = eqOffset;
            testParams.statTextYPosStepSize = statTextYPosStepSize;
            testParams.statStarYLocOffset = statStarYLocOffset;
            testParams.colorMap = colorMap;
            testParams.mixedEffectPValues = mixedEffectPValues;
            testParams.axHandle = handles.ax;
            testParams.groupOffsets = groupOffsets;
            testParams.isWithinGroupMode = isWithinGroupMode;
            
            % Iterate through each pair of conditions and plot statistics
            [~, pTsOut, pRsOut, hTestsOut] = arrayfun(@(idx) evaluate_and_plot_2sample_tests(cond1Values(idx), cond2Values(idx), idx, testParams), (1:length(cond1Values))', 'UniformOutput', false);

            % Store the concatenated p-values into the stats output structure
            statsOut.pValueT = cat(2, pTsOut{:});
            statsOut.pValueR = cat(2, pRsOut{:});

            % Accumulate graphic handles from 2-sample tests
            for iTest = 1:numel(hTestsOut)
                handles.pText = [handles.pText; hTestsOut{iTest}.pText(:)];
                handles.sigMarker = [handles.sigMarker; hTestsOut{iTest}.sigMarker(:)];
            end

        elseif nConditions == 1 || isSingleConditionPerGroup
            % Single condition testing against a null hypothesis
            
            % Pre-calculate Y relative positions so they horizontally align across all groups/conditions
            yRelT = 1 - statTextYPosStepSize * (1 + eqOffset);
            yRelR = 1 - statTextYPosStepSize * (1 + runTTest + eqOffset);
            
            % Pre-allocate outputs for 1-sample tests
            if isWithinGroupMode
                statsOut.pValueT = NaN(nGroups, 1);
                statsOut.pValueR = NaN(nGroups, 1);
            else
                statsOut.pValueT = NaN(1, 1);
                statsOut.pValueR = NaN(1, 1);
            end
            
            % Execute tests separated by group
            if isWithinGroupMode
                % Loop through each group and test separately
                for iGroup = 1:nGroups
                    % Extract the value for the current group
                    currentGroupValue = uniqueGroupValues(iGroup);
                    
                    % Find the condition for this group
                    groupConditionMask = grouping == currentGroupValue;
                    condValue = unique(condition(groupConditionMask));
                    
                    % Skip if no condition
                    if isempty(condValue)
                        continue;
                    end
                    condValue = condValue(1); 
                    
                    % Get data for this group in the specific condition
                    groupData = data(groupConditionMask & (condition == condValue));

                    % Skip testing if data is missing for this group
                    if isempty(groupData)
                        continue; 
                    end

                    % Evaluate metrics against null using subfunction 
                    [pT, pR, isAppropriateForTTest] = compute_1sample_stats(groupData, dataType);
                    
                    % Store computed group p-values
                    statsOut.pValueT(iGroup, 1) = pT;
                    statsOut.pValueR(iGroup, 1) = pR;

                    % Get the color for this group
                    groupColor = colorMap{iGroup};

                    % Calculate the x-axis position for the group stat to align exactly with the jitter spread
                    xPosText = condValue + groupOffsets(iGroup);

                    % Run and plot parametric test if requested
                    if runTTest && ~isnan(pT)
                        hT = plot_test_result(pT, 'PString', 'p_t', ...
                                    'XLocText', xPosText, 'XLocStar', xPosText, ...
                                    'YLocTextRel', yRelT, 'YLocStarRel', yRelT + statStarYLocOffset, ...
                                    'SigLevel', sigLevel, 'IsAppropriate', isAppropriateForTTest, ...
                                    'ColorSignificant', groupColor, ...
                                    'AxesHandle', handles.ax);

                        % Accumulate returned handles into the main handles structure
                        handles.pText = [handles.pText; hT.pText(:)];
                        handles.sigMarker = [handles.sigMarker; hT.sigMarker(:)];
                    end

                    % Run and plot non-parametric test if requested
                    if runRankTest && ~isnan(pR)
                        hR = plot_test_result(pR, 'PString', 'p_r', ...
                                    'XLocText', xPosText, 'XLocStar', xPosText, ...
                                    'YLocTextRel', yRelR, 'YLocStarRel', yRelR + statStarYLocOffset, ...
                                    'SigLevel', sigLevel, 'IsAppropriate', ~isAppropriateForTTest, ...
                                    'ColorSignificant', groupColor, ...
                                    'AxesHandle', handles.ax);

                        % Accumulate returned handles into the main handles structure
                        handles.pText = [handles.pText; hR.pText(:)];
                        handles.sigMarker = [handles.sigMarker; hR.sigMarker(:)];
                    end
                end
            else
                % Pool groups and test single condition against null or use mixed effects
                condValue = uniqueConditionValues(1);
                isCond = (condition == condValue);
                condData = data(isCond);
                
                % Define x-position for the text directly above the condition
                xPosText = condValue;

                % Execute tests if data is present
                if ~isempty(condData)
                    % Evaluate metrics against null using subfunction 
                    if strcmpi(statsMode, 'mixedEffect')
                        pT = mixedEffectPValues(1);
                        pR = NaN;
                        isAppropriateForTTest = true;
                    else
                        [pT, pR, isAppropriateForTTest] = compute_1sample_stats(condData, dataType);
                    end
                    
                    % Store computed pooled p-values
                    statsOut.pValueT(1, 1) = pT;
                    statsOut.pValueR(1, 1) = pR;

                    % Run and plot parametric test if requested
                    if runTTest && ~isnan(pT)
                        hT = plot_test_result(pT, 'PString', 'p_t', ...
                                    'XLocText', xPosText, 'XLocStar', xPosText, ...
                                    'YLocTextRel', yRelT, 'YLocStarRel', yRelT + statStarYLocOffset, ...
                                    'SigLevel', sigLevel, 'IsAppropriate', isAppropriateForTTest, ...
                                    'AxesHandle', handles.ax);

                        % Accumulate returned handles into the main handles structure
                        handles.pText = [handles.pText; hT.pText(:)];
                        handles.sigMarker = [handles.sigMarker; hT.sigMarker(:)];
                    end

                    % Run and plot non-parametric test if requested
                    if runRankTest && ~isnan(pR)
                        hR = plot_test_result(pR, 'PString', 'p_r', ...
                                    'XLocText', xPosText, 'XLocStar', xPosText, ...
                                    'YLocTextRel', yRelR, 'YLocStarRel', yRelR + statStarYLocOffset, ...
                                    'SigLevel', sigLevel, 'IsAppropriate', ~isAppropriateForTTest, ...
                                    'AxesHandle', handles.ax);

                        % Accumulate returned handles into the main handles structure
                        handles.pText = [handles.pText; hR.pText(:)];
                        handles.sigMarker = [handles.sigMarker; hR.sigMarker(:)];
                    end
                end
            end
        end

        % Evaluate tests across adjacent groups within the same condition when requested
        if strcmpi(statsMode, 'withinGroupAndCond') && nGroups >= 2
            % Pass parameters down
            testParams.data = data;
            testParams.condition = condition;
            testParams.grouping = grouping;
            testParams.uniqueGroupValues = uniqueGroupValues;
            testParams.nGroups = nGroups;
            testParams.runTTest = runTTest;
            testParams.runRankTest = runRankTest;
            testParams.dataType = dataType;
            testParams.sigLevel = sigLevel;
            testParams.yLims = yLims;
            testParams.eqOffset = eqOffset;
            testParams.statTextYPosStepSize = statTextYPosStepSize;
            testParams.statStarYLocOffset = statStarYLocOffset;
            testParams.axHandle = handles.ax;
            testParams.groupOffsets = groupOffsets;
            testParams.linesForCondTests = linesForCondTests;
            
            % Iterate through conditions and test adjacent groups
            [~, pTsBgOut, pRsBgOut, hTestsBgOut] = arrayfun(@(idx) evaluate_and_plot_2sample_groups_within_cond(uniqueConditionValues(idx), idx, testParams), (1:length(uniqueConditionValues))', 'UniformOutput', false);

            % Export between-group comparisons natively out to stats structure
            statsOut.pValueT_betweenGroups = cat(2, pTsBgOut{:});
            statsOut.pValueR_betweenGroups = cat(2, pRsBgOut{:});

            % Accumulate graphic handles from between-group tests
            for iTest = 1:numel(hTestsBgOut)
                handles.pText = [handles.pText; hTestsBgOut{iTest}.pText(:)];
                handles.sigMarker = [handles.sigMarker; hTestsBgOut{iTest}.sigMarker(:)];
            end
        end
    end
end

% Hold off
hold_off(wasHold, axHandle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [status, pTs, pRs, hTestOut] = evaluate_and_plot_2sample_tests(cond1Value, cond2Value, pairIdx, p)
%% Evaluates and plots 2-sample tests for a given pair of conditions

% Initialize dummy status variable and pre-allocate outputs
status = true;
hTestOut.pText = gobjects(0, 1);
hTestOut.sigMarker = gobjects(0, 1);

if p.isWithinGroupMode
    pTs = NaN(p.nGroups, 1);
    pRs = NaN(p.nGroups, 1);
else
    pTs = NaN(1, 1);
    pRs = NaN(1, 1);
end

% Set the initial baseline for y positions depending on the previously established limits and offsets
yRange = diff(p.yLims);
yPos = p.yLims(2) - p.eqOffset * p.statTextYPosStepSize * yRange; 

% Execute tests separated by group
if p.isWithinGroupMode
    % Loop through each group and test separately
    for iGroup = 1:p.nGroups
        % Extract the value for the current group
        currentGroupValue = p.uniqueGroupValues(iGroup);
        
        % Get data for this group in condition 1 & 2
        group1Data = p.data((p.grouping == currentGroupValue) & (p.condition == cond1Value));
        group2Data = p.data((p.grouping == currentGroupValue) & (p.condition == cond2Value));

        % Skip testing if data is missing for this group
        if isempty(group1Data) || isempty(group2Data)
            continue; 
        end

        % Evaluate metrics using subfunction 
        [pT, pR, isAppropriateForTTest] = compute_2sample_stats(group1Data, group2Data, p.dataType);
        
        % Store group p-values
        pTs(iGroup) = pT;
        pRs(iGroup) = pR;

        % Get the color for this group
        groupColor = p.colorMap{iGroup};

        % Calculate x locations for the two conditions being compared to connect lines
        xLoc1 = cond1Value + p.groupOffsets(iGroup);
        xLoc2 = cond2Value + p.groupOffsets(iGroup);
        xPosText = mean([xLoc1, xLoc2]);

        % Run and plot parametric test if requested
        if p.runTTest && ~isnan(pT)
            % Calculate relative Y position for the parametric test text
            yPos = yPos - p.statTextYPosStepSize * yRange; 
            yRel = (yPos - p.yLims(1)) / yRange;
            
            % Plot the parametric test result
            hT = plot_test_result(pT, 'PString', 'p_t', ...
                        'XLocText', xPosText, 'XLocStar', xPosText, ...
                        'YLocTextRel', yRel, 'YLocStarRel', yRel + p.statStarYLocOffset, ...
                        'SigLevel', p.sigLevel, 'IsAppropriate', isAppropriateForTTest, ...
                        'ColorSignificant', groupColor, ...
                        'AxesHandle', p.axHandle);

            % Accumulate returned handles into the local handles structure
            hTestOut.pText = [hTestOut.pText; hT.pText(:)];
            hTestOut.sigMarker = [hTestOut.sigMarker; hT.sigMarker(:)];

            % Draw a horizontal line below the p-value text to show which conditions are compared
            yLineData = yPos - 0.015 * yRange;
            plot(p.axHandle, [xLoc1, xLoc2], [yLineData, yLineData], 'Color', groupColor, 'LineWidth', 1, 'HandleVisibility', 'off');
        end

        % Run and plot non-parametric test if requested
        if p.runRankTest && ~isnan(pR)
            % Calculate relative Y position for the non-parametric test text
            yPos = yPos - p.statTextYPosStepSize * yRange; 
            yRel = (yPos - p.yLims(1)) / yRange;

            % Plot the non-parametric test result
            hR = plot_test_result(pR, 'PString', 'p_r', ...
                        'XLocText', xPosText, 'XLocStar', xPosText, ...
                        'YLocTextRel', yRel, 'YLocStarRel', yRel + p.statStarYLocOffset, ...
                        'SigLevel', p.sigLevel, 'IsAppropriate', ~isAppropriateForTTest, ...
                        'ColorSignificant', groupColor, ...
                        'AxesHandle', p.axHandle);

            % Accumulate returned handles into the local handles structure
            hTestOut.pText = [hTestOut.pText; hR.pText(:)];
            hTestOut.sigMarker = [hTestOut.sigMarker; hR.sigMarker(:)];

            % Draw a horizontal line below the p-value text to show which conditions are compared
            yLineData = yPos - 0.015 * yRange;
            plot(p.axHandle, [xLoc1, xLoc2], [yLineData, yLineData], 'Color', groupColor, 'LineWidth', 1, 'HandleVisibility', 'off');
        end
    end
else
    % Pool groups and test between conditions or use precomputed mixed effects
    cond1Mask = (p.condition == cond1Value);
    cond2Mask = (p.condition == cond2Value);
    cond1Data = p.data(cond1Mask);
    cond2Data = p.data(cond2Mask);
    
    % Calculate x locations for the two conditions being compared for drawing the connecting lines
    xLoc1 = cond1Value;
    xLoc2 = cond2Value;
    xPosText = mean([xLoc1, xLoc2]);

    % Execute tests if both conditions have data
    if ~isempty(cond1Data) && ~isempty(cond2Data)
        % Evaluate metrics using subfunction 
        if strcmpi(p.statsMode, 'mixedEffect')
            pT = p.mixedEffectPValues(pairIdx);
            pR = NaN;
            isAppropriateForTTest = true;
        else
            [pT, pR, isAppropriateForTTest] = compute_2sample_stats(cond1Data, cond2Data, p.dataType);
        end
        
        % Store pooled p-values
        pTs(1) = pT;
        pRs(1) = pR;

        % Run and plot parametric test if requested
        if p.runTTest && ~isnan(pT)
            % Calculate relative Y position for the parametric test text
            yPos = yPos - p.statTextYPosStepSize * yRange; 
            yRel = (yPos - p.yLims(1)) / yRange;
            
            % Plot the parametric test result
            hT = plot_test_result(pT, 'PString', 'p_t', ...
                        'XLocText', xPosText, 'XLocStar', xPosText, ...
                        'YLocTextRel', yRel, 'YLocStarRel', yRel + p.statStarYLocOffset, ...
                        'SigLevel', p.sigLevel, 'IsAppropriate', isAppropriateForTTest, ...
                        'AxesHandle', p.axHandle);

            % Accumulate returned handles into the local handles structure
            hTestOut.pText = [hTestOut.pText; hT.pText(:)];
            hTestOut.sigMarker = [hTestOut.sigMarker; hT.sigMarker(:)];

            % Draw a horizontal line below the p-value text to show which conditions are compared
            yLineData = yPos - 0.015 * yRange;
            plot(p.axHandle, [xLoc1, xLoc2], [yLineData, yLineData], 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
        end

        % Run and plot non-parametric test if requested
        if p.runRankTest && ~isnan(pR)
            % Calculate relative Y position for the non-parametric test text
            yPos = yPos - p.statTextYPosStepSize * yRange; 
            yRel = (yPos - p.yLims(1)) / yRange;

            % Plot the non-parametric test result
            hR = plot_test_result(pR, 'PString', 'p_r', ...
                        'XLocText', xPosText, 'XLocStar', xPosText, ...
                        'YLocTextRel', yRel, 'YLocStarRel', yRel + p.statStarYLocOffset, ...
                        'SigLevel', p.sigLevel, 'IsAppropriate', ~isAppropriateForTTest, ...
                        'AxesHandle', p.axHandle);

            % Accumulate returned handles into the local handles structure
            hTestOut.pText = [hTestOut.pText; hR.pText(:)];
            hTestOut.sigMarker = [hTestOut.sigMarker; hR.sigMarker(:)];

            % Draw a horizontal line below the p-value text to show which conditions are compared
            yLineData = yPos - 0.015 * yRange;
            plot(p.axHandle, [xLoc1, xLoc2], [yLineData, yLineData], 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [status, pTs, pRs, hTestOut] = evaluate_and_plot_2sample_groups_within_cond(condValue, condIdx, p)
%% Evaluates and plots 2-sample tests between adjacent groups within a given condition

% Initialize dummy status variable and pre-allocate outputs
status = true;
hTestOut.pText = gobjects(0, 1);
hTestOut.sigMarker = gobjects(0, 1);

pTs = NaN(p.nGroups - 1, 1);
pRs = NaN(p.nGroups - 1, 1);

% Set the initial baseline for y positions depending on the previously established limits and offsets 
yRange = diff(p.yLims);
yPos = p.yLims(2) - (p.eqOffset + p.linesForCondTests) * p.statTextYPosStepSize * yRange; 

% Loop through adjacent groups and test separately
for iGroup = 1:p.nGroups-1
    % Extract the values for the two groups
    g1Val = p.uniqueGroupValues(iGroup);
    g2Val = p.uniqueGroupValues(iGroup+1);
    
    % Get data for the two groups in the current condition
    group1Data = p.data((p.grouping == g1Val) & (p.condition == condValue));
    group2Data = p.data((p.grouping == g2Val) & (p.condition == condValue));

    % Skip testing if data is missing for either group
    if isempty(group1Data) || isempty(group2Data)
        continue; 
    end

    % Evaluate metrics using subfunction 
    [pT, pR, isAppropriateForTTest] = compute_2sample_stats(group1Data, group2Data, p.dataType);
    
    % Store group p-values
    pTs(iGroup) = pT;
    pRs(iGroup) = pR;

    % Calculate x locations for the two groups being compared to connect lines
    xLoc1 = condValue + p.groupOffsets(iGroup);
    xLoc2 = condValue + p.groupOffsets(iGroup+1);
    xPosText = mean([xLoc1, xLoc2]);

    % Run and plot parametric test if requested
    if p.runTTest && ~isnan(pT)
        % Calculate relative Y position for the parametric test text
        yPos = yPos - p.statTextYPosStepSize * yRange; 
        yRel = (yPos - p.yLims(1)) / yRange;
        
        % Plot the parametric test result
        hT = plot_test_result(pT, 'PString', 'p_t', ...
                    'XLocText', xPosText, 'XLocStar', xPosText, ...
                    'YLocTextRel', yRel, 'YLocStarRel', yRel + p.statStarYLocOffset, ...
                    'SigLevel', p.sigLevel, 'IsAppropriate', isAppropriateForTTest, ...
                    'AxesHandle', p.axHandle);

        % Accumulate returned handles into the local handles structure
        hTestOut.pText = [hTestOut.pText; hT.pText(:)];
        hTestOut.sigMarker = [hTestOut.sigMarker; hT.sigMarker(:)];

        % Draw a horizontal line below the p-value text to show which groups are compared
        yLineData = yPos - 0.015 * yRange;
        plot(p.axHandle, [xLoc1, xLoc2], [yLineData, yLineData], 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
    end

    % Run and plot non-parametric test if requested
    if p.runRankTest && ~isnan(pR)
        % Calculate relative Y position for the non-parametric test text
        yPos = yPos - p.statTextYPosStepSize * yRange; 
        yRel = (yPos - p.yLims(1)) / yRange;

        % Plot the non-parametric test result
        hR = plot_test_result(pR, 'PString', 'p_r', ...
                    'XLocText', xPosText, 'XLocStar', xPosText, ...
                    'YLocTextRel', yRel, 'YLocStarRel', yRel + p.statStarYLocOffset, ...
                    'SigLevel', p.sigLevel, 'IsAppropriate', ~isAppropriateForTTest, ...
                    'AxesHandle', p.axHandle);

        % Accumulate returned handles into the local handles structure
        hTestOut.pText = [hTestOut.pText; hR.pText(:)];
        hTestOut.sigMarker = [hTestOut.sigMarker; hR.sigMarker(:)];

        % Draw a horizontal line below the p-value text to show which groups are compared
        yLineData = yPos - 0.015 * yRange;
        plot(p.axHandle, [xLoc1, xLoc2], [yLineData, yLineData], 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [condMeans, condLower, condUpper, pairPTs, success, fittedModel, modelEquation] = compute_global_mixed_effects(data, condition, group, dataType, errorBarType, uniqueConds, doTests)
%% Fits a global mixed effect model to estimate marginal means and pairwise p-values

% Initialize default return variables
pairPTs = [];
success = false;
fittedModel = [];
modelEquation = '';

% Create logical index of valid data points (no NaNs)
validIdx = ~isnan(data) & ~isnan(condition) & ~isnan(group);

% Filter data, condition, and group to remove NaNs
data = data(validIdx);
condition = condition(validIdx);
group = group(validIdx);

% Calculate total number of unique conditions
numConds = length(uniqueConds);

% Initialize condition limits arrays
condMeans = NaN(numConds, 1);
condLower = NaN(numConds, 1);
condUpper = NaN(numConds, 1);

% Initialize pairwise p-values if tests are requested
if doTests
    if numConds >= 2
        pairPTs = NaN(numConds - 1, 1);
    else
        pairPTs = NaN;
    end
end

% Return early if data is insufficient for mixed effect analysis
if numel(data) < 2 || numel(unique(group)) < 2
    return;
end

% Attempt to construct and fit the mixed effect model
try
    % Build a standardized categorical table forcing uniform condition order
    conditionCategorical = categorical(condition, uniqueConds);
    groupCategorical = categorical(group);
    dataTable = table(data, conditionCategorical, groupCategorical, 'VariableNames', {'y', 'condition', 'g'});

    % Choose the fitting procedure based on the data type
    switch lower(dataType)
        case 'normal'
            % Fit linear mixed-effects model for normal data
            fittedModel = fitlme(dataTable, 'y ~ condition + (1|g)');
        case 'proportion'
            % Fit generalized linear mixed-effects model for binomial data
            fittedModel = fitglme(dataTable, 'y ~ condition + (1|g)', 'Distribution', 'Binomial');
        case 'ratio'
            % Log-transform the response variable for ratio data
            dataTable.y = log(dataTable.y);

            % Fit linear mixed-effects model on log-transformed data
            fittedModel = fitlme(dataTable, 'y ~ condition + (1|g)');
    end

    % Extract degrees of freedom from the fitted model
    degFree = fittedModel.DFE;

    % Estimate Marginal Means and Errors by querying model directly
    newTable = table(categorical(uniqueConds, uniqueConds), 'VariableNames', {'condition'});

    % Assign a dummy group for prediction (marginalized later)
    newTable.g = repmat(groupCategorical(1), numConds, 1);

    % Handle predictions based on the underlying data distribution
    if strcmpi(dataType, 'proportion')
        % Predict conditional on fixed effects only for proportion data
        [yPred, yCi] = predict(fittedModel, newTable, 'Conditional', false);
        condMeans = yPred;
        
        % Calculate final error boundaries based on the requested error bar type
        if strcmpi(errorBarType, 'err95')
            condLower = yCi(:,1);
            condUpper = yCi(:,2);
        else
            % Approximate symmetric error for standard error bar plotting
            stdErr = (yCi(:,2) - yCi(:,1)) / (2 * 1.95996); % 95% CI Z-score approximation
            condLower = condMeans - stdErr;
            condUpper = condMeans + stdErr;
        end
    elseif strcmpi(dataType, 'ratio')
        % Predict conditional on fixed effects only for log-scale data
        [yPred, yCi] = predict(fittedModel, newTable, 'Conditional', false);
 
        % Exponentiate predictions to return to the original ratio scale
        condMeans = exp(yPred);
        
        % Calculate final error boundaries based on the requested error bar type
        if strcmpi(errorBarType, 'err95')
            condLower = exp(yCi(:,1));
            condUpper = exp(yCi(:,2));
        else
            % Get critical t-value for 95% confidence interval
            tValue = tinv(0.975, degFree);

            % Estimate standard error on the log scale
            stdErrLog = (yCi(:,2) - yCi(:,1)) / (2 * tValue);

            % Transform error back using the Delta method
            stdErr = condMeans .* stdErrLog;
            condLower = condMeans - stdErr;
            condUpper = condMeans + stdErr;
        end
    else
        % Predict conditional on fixed effects only for standard normal data
        [yPred, yCi] = predict(fittedModel, newTable, 'Conditional', false);

        % Set the direct predictions as the marginal means
        condMeans = yPred;
        
        % Calculate final error boundaries based on the requested error bar type
        if strcmpi(errorBarType, 'err95')
            condLower = yCi(:,1);
            condUpper = yCi(:,2);
        else
            % Get critical t-value for 95% confidence interval
            tValue = tinv(0.975, degFree);

            % Estimate standard error using the t-value
            stdErr = (yCi(:,2) - yCi(:,1)) / (2 * tValue);
            condLower = condMeans - stdErr;
            condUpper = condMeans + stdErr;
        end
    end

    % Extract Pairwise P-Values mathematically evaluating contrast coefficient subsets 
    if doTests
        % Perform tests based on the number of conditions
        if numConds >= 2
            for iCond = 1:(numConds - 1)
                % Initialize contrast vector for the hypothesis test
                contrastVector = zeros(1, fittedModel.NumCoefficients);

                % Configure contrast coefficients based on pair sequence
                if iCond == 1
                    % Pair 1 matches condition 2 relative to intercept
                    contrastVector(2) = 1;
                else
                    % Adjacent Pair evaluation isolates difference delta
                    contrastVector(iCond) = -1;
                    contrastVector(iCond+1) = 1;
                end

                % Compute p-value using the contrast coefficient test
                pairPTs(iCond) = coefTest(fittedModel, contrastVector);
            end
        else
            % Perform single-sample test to see if the intercept coefficient significantly differs from zero
            pairPTs(1) = coefTest(fittedModel, 1, 0);
        end
    end

    % Extract fixed effects coefficients to construct the equation string
    coefNames = fittedModel.CoefficientNames;
    coefVals = fittedModel.fixedEffects;

    % Initialize cell array for formatted equation terms
    equationTerms = cell(length(coefNames), 1);

    % Format each term in the equation string
    for iCoef = 1:length(coefNames)
        if strcmp(coefNames{iCoef}, '(Intercept)')
            % Insert intercept value as a base number
            equationTerms{iCoef} = sprintf('%.3f', coefVals(iCoef));
        else
            % Remove the categorical prefix for a cleaner variable string
            cleanName = strrep(coefNames{iCoef}, 'condition_', 'cond');
            equationTerms{iCoef} = sprintf('%.3f*%s', coefVals(iCoef), cleanName);
        end
    end

    % Combine terms and append the random effects notation
    modelEquation = ['y = ', strjoin(equationTerms, ' + '), ' + (1|group)'];

    % Mark the overall fitting and extraction process as successful
    success = true;
catch ME
    % Silently fail output to trigger immediate fallback procedure upstream
    disp(['Mixed effects model failed to converge globally: ', ME.message]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pT, pR, isAppropriateForTTest] = compute_1sample_stats(groupData, dataType)
%% Computes one-sample statistical tests against a null hypothesis based on data type

% Pre-allocate statistical outputs 
pT = NaN; 
pR = NaN; 
isAppropriateForTTest = false;

% Remove NaNs for proper length checking and calculations
validData = groupData(~isnan(groupData));

% Execute the appropriate one-sample test based on the specified data type
switch lower(dataType)
    case 'normal'
        % Null hypothesis is 0 for standard metrics
        nullVal = 0;
        
        % Validate normality assumption 
        isAppropriateForTTest = test_normality(validData);
        
        % Check if sufficient data length exists to perform calculations
        if numel(validData) >= 2
            [~, pT] = ttest(validData, nullVal);
            pR = signrank(validData, nullVal);
        end
        
    case 'ratio'
        % Null hypothesis is 1 for ratio comparisons
        nullVal = 1;
        
        % Validate normality on the logarithmic scale 
        logData = log(validData);
        isAppropriateForTTest = test_normality(logData);
        
        % Check if sufficient data length exists to perform calculations
        if numel(validData) >= 2
            % Perform parametric test against log(1) == 0
            [~, pT] = ttest(logData, log(nullVal));
            
            % Perform non-parametric test against raw ratio == 1
            pR = signrank(validData, nullVal);
        end
        
    case 'proportion'
        % Null hypothesis is 0.5 (random chance)
        nullVal = 0.5;
        numTrials = numel(validData);
        numSuccesses = sum(validData == 1); 
        
        % Use heuristic for normal distribution approximation appropriateness
        isAppropriateForTTest = (numTrials * nullVal >= 5) && (numTrials * (1 - nullVal) >= 5);
        
        % Check if sufficient data length exists to perform calculations
        if numTrials >= 1
            % Parametric equivalent (Z-test)
            zStat = (numSuccesses / numTrials - nullVal) / sqrt(nullVal * (1 - nullVal) / numTrials);
            pT = 2 * normcdf(-abs(zStat));
            
            % Non-parametric equivalent (Exact Binomial test)
            pR = 2 * min(binocdf(numSuccesses, numTrials, nullVal), ...
                         1 - binocdf(numSuccesses - 1, numTrials, nullVal));
            pR = min(pR, 1); % Ensure probabilities don't surpass 100%
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pT, pR, isAppropriateForTTest] = compute_2sample_stats(group1Data, group2Data, dataType)
%% Computes two-sample statistical tests between distributions based on data type

% Pre-allocate statistical outputs 
pT = NaN; 
pR = NaN; 
isAppropriateForTTest = false;

% Remove NaNs for proper length checking and calculations
validData1 = group1Data(~isnan(group1Data));
validData2 = group2Data(~isnan(group2Data));

switch lower(dataType)
    case 'normal'
        % Validate normality assumption in both independent groups
        isAppropriateForTTest = test_normality(validData1) && test_normality(validData2);
        
        % Check if sufficient data length exists to perform calculations
        if numel(validData1) >= 2 && numel(validData2) >= 2
            [~, pT] = ttest2(validData1, validData2);
            pR = ranksum(validData1, validData2);
        end
        
    case 'ratio'
        % Validate normality assumption in both independent groups on logarithmic scale
        logData1 = log(validData1);
        logData2 = log(validData2);
        isAppropriateForTTest = test_normality(logData1) && test_normality(logData2);
        
        % Check if sufficient data length exists to perform calculations
        if numel(validData1) >= 2 && numel(validData2) >= 2
            [~, pT] = ttest2(logData1, logData2);
            pR = ranksum(validData1, validData2);
        end
        
    case 'proportion'
        % Count trials and successes for each group
        numTrials1 = numel(validData1);
        numTrials2 = numel(validData2);
        numSuccesses1 = sum(validData1 == 1);
        numSuccesses2 = sum(validData2 == 1);
        
        % Check if sufficient data length exists to perform calculations
        if numTrials1 >= 1 && numTrials2 >= 1
            % Calculate pooled probability for the Z-test
            pooledProb = (numSuccesses1 + numSuccesses2) / (numTrials1 + numTrials2);
            
            % Use heuristic for normal distribution approximation appropriateness
            isAppropriateForTTest = (numTrials1 * pooledProb >= 5) && (numTrials1 * (1 - pooledProb) >= 5) && ...
                                    (numTrials2 * pooledProb >= 5) && (numTrials2 * (1 - pooledProb) >= 5);
            
            % Parametric equivalent (Two-proportion Z-test)
            if pooledProb == 0 || pooledProb == 1
                pT = 1; % Data pools perfectly without deviation
            else
                standardError = sqrt(pooledProb * (1 - pooledProb) * (1 / numTrials1 + 1 / numTrials2));
                zStat = ((numSuccesses1 / numTrials1) - (numSuccesses2 / numTrials2)) / standardError;
                pT = 2 * normcdf(-abs(zStat));
            end
            
            % Non-parametric equivalent (Fisher's exact test)
            contingencyTable = [numSuccesses1, numTrials1 - numSuccesses1; ...
                                numSuccesses2, numTrials2 - numSuccesses2];
            try
                [~, pR] = fishertest(contingencyTable);
            catch
                pR = NaN; % Fallback if exact testing fails computation matrix criteria
            end
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function legendLocation = set_default_legend_location(nGroups, maxInner, maxOuter) 
%% Sets the default legend location based on the number of groups

if nGroups > 1 && nGroups <= maxInner
    legendLocation = 'best';
elseif nGroups > maxInner && nGroups <= maxOuter
    legendLocation = 'eastoutside';
else
    legendLocation = 'suppress';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%