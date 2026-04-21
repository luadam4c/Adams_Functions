function handles = plot_grouped_jitter (data, varargin)
%% Plots a jitter plot colored by group from data (uses plotSpread)
% Usage: handles = plot_grouped_jitter (data, grouping (opt), condition (opt), varargin)
% Explanation:
%       This function creates a grouped jitter plot (swarm plot) where data 
%       points are spread along the x-axis to avoid overlap. It supports two 
%       levels of grouping: 
%           1. 'Condition' (mapped to x-axis positions)
%           2. 'Grouping' (mapped to color)
%       It optionally overlays descriptive statistics (mean and error bars) 
%       and performs statistical tests (unpaired t-test or rank-sum test) 
%       between conditions, either pooled across groups or separated by group.
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
%           'XTickLabels', {'Baseline', 'Post-Stim'});
%
%       % Example 3: Customizing limits and statistics
%       plot_grouped_jitter(data, grouping, condition, ...
%           'YLimits', [-5, 10], ...
%           'RunTTest', true, ...
%           'StatisticsByGroup', true, ...
%           'PlotMeanValues', true);
%
%       randVec1 = randi(10, 10, 1);
%       randVec2 = randi(10, 10, 1) + 10;
%       data = [randVec1, randVec2];
%       plot_grouped_jitter(data)
%       plot_grouped_jitter(data, 'UsePlotSpread', true)
%       plot_grouped_jitter(data, 'UsePlotSpread', false)
%       plot_grouped_jitter(data, 'UsePlotSpread', false, 'JitterWidth', 0.5)
%
%       data = [randn(50,1);randn(50,1)+3.5]*[1 1];
%       grouping = [[ones(50,1);zeros(50,1)],[randi([0,1],[100,1])]];
%       plot_grouped_jitter(data, grouping, 'RunTTest', true)
%       plot_grouped_jitter(data, grouping, 'UsePlotSpread', false, 'RunRankTest', true)
%
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
%                                           perform TODO
%                       'ratio'         - expect ratios
%                                           compute TODO
%                                           perform TODO
%                   default == 'normal'
%                   - 'Weights': optional weight vector for calculating 
%                       weighted standard errors
%                   must be empty or a numeric array of the same size as data
%                   default == []
%                   - 'Means': optional means to plot instead of calculating
%                   must be a numeric array 
%                   default == []
%                   - 'Errors': optional errors to plot instead of calculating
%                   must be a numeric array 
%                   default == []
%                   - 'UsePlotSpread': whether to use plotSpread.m
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == [] (true if plotSpread.m is found)
%                   - 'JitterWidth': width of the jitter
%                   must be a non-negative scalar
%                   default == [] (plotSpread default) or 0.3 (manual plot)
%                   - 'PlotMeanValues': whether to plot the mean values
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotErrorBars': whether to plot error bars 
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'StatisticsByGroup': whether to compute statistics 
%                                          by group (true) or pooled (false)
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
%                               'northeast' if nGroups is 2~9
%                               'eastoutside' if nGroups is 10+
%                   - 'AxesHandle': axes handle to plot on
%                   must be a empty or an axes object handle
%                   default == set in set_axes_properties.m
%                   - Any other parameter-value pair for plotSpread() or plot()
%
% Requires:
%       cd/addpath_custom.m
%       cd/create_default_grouping.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/decide_on_colormap.m
%       cd/hold_off.m
%       cd/hold_on.m
%       cd/locate_functionsdir.m
%       cd/nanstderr.m
%       cd/plot_test_result.m
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

%% Define hard-coded parameters for the function
validDataTypes = {'normal', 'proportion'};
maxNGroupsForInnerLegend = 8;
maxNGroupsForOuterLegends = 25;
forceVectorInput = true;       % Consider making this into an optional argument
defaultJitterWidthNotPlotSpread = 0.3;
xTickLimitPadding = 0.5;
groupMeanSpreadWidth = 0.4;
pooledStatColor = 'k';
statTextYPosStepSize = 0.1;
statStarYLocOffset = 0.05;

%% Set default values for optional arguments
groupingDefault = [];           % set later
distributionDefault = [];       % set later
dataTypeDefault = 'normal';
weightsDefault = [];
meansDefault = [];
errorsDefault = [];
usePlotSpreadDefault = [];      % set later
jitterWidthDefault = [];        % set later
plotMeanValuesDefault = true;
plotErrorBarsDefault = true;
statisticsByGroupDefault = true;
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
addParameter(iP, 'Weights', weightsDefault, ...
    @(x) assert(isempty(x) || isnumeric(x), 'Weights must be empty or a numeric array!'));
addParameter(iP, 'Means', meansDefault, ...
    @(x) isempty(x) || isnumeric(x));
addParameter(iP, 'Errors', errorsDefault, ...
    @(x) isempty(x) || isnumeric(x));
addParameter(iP, 'UsePlotSpread', usePlotSpreadDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'JitterWidth', jitterWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'PlotMeanValues', plotMeanValuesDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotErrorBars', plotErrorBarsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'StatisticsByGroup', statisticsByGroupDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
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
weights = iP.Results.Weights;
means = iP.Results.Means;
errors = iP.Results.Errors;
usePlotSpread = iP.Results.UsePlotSpread;
jitterWidth = iP.Results.JitterWidth;
plotMeanValues = iP.Results.PlotMeanValues;
plotErrorBars = iP.Results.PlotErrorBars;
statisticsByGroup = iP.Results.StatisticsByGroup;
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

% Set legend location based on number of groups
% TODO: Use set_default_legend_location.m
if strcmpi(legendLocation, 'auto')
    if nGroups > 1 && nGroups <= maxNGroupsForInnerLegend
        legendLocation = 'northeast';
    elseif nGroups > maxNGroupsForInnerLegend && ...
            nGroups <= maxNGroupsForOuterLegends
        legendLocation = 'eastoutside';
    else
        legendLocation = 'suppress';
    end
end

% Decide on the axes to plot on
axHandle = set_axes_properties('AxesHandle', axHandle);

%% Do the job
% Return immediately if there is no data to plot
if isempty(data)
    handles = struct;
    return
end

% Retain the current plot on the axes
wasHold = hold_on(axHandle);

% Plot the data points using plotSpread or manually based on preference
if usePlotSpread
    % Prevent plotSpread from plotting its own means if we are plotting custom ones
    if plotMeanValues || plotErrorBars
        otherArguments = [otherArguments, {'showMM', 0}];
    end

    % Call plotSpread to plot the swarm plot
    output = plotSpread(axHandle, data, 'distributionIdx', condition, ...
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
else
    % Give each condition index some jitter for plotting
    conditionWithJitter = condition + (jitterWidth * (rand(nPoints, 1) - 0.5));

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

% Modify x limits
if ~(ischar(xLimits) && strcmpi(xLimits, 'suppress'))
    xlim(axHandle, xLimits);
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

%% Plot statistics
% Plot means and error bars for each condition
if plotMeanValues || plotErrorBars
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
        if statisticsByGroup
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

                % Calculate the x-axis positions for the group means
                if nGroupsInThisCond == 1
                    xMeanPositions = currentCondValue;
                else
                    offsets = linspace(-groupMeanSpreadWidth/2, groupMeanSpreadWidth/2, nGroupsInThisCond);
                    xMeanPositions = currentCondValue + offsets;
                end

                % Loop through and plot statistics for each group
                for iGroup = 1:nGroupsInThisCond
                    % Extract the data and weights for the current group
                    groupData = dataByGroup{iGroup};
                    groupWeights = weightsByGroup{iGroup};
                    
                    % Find the true index of the current group
                    actualGroupIdx = find(uniqueGroupValues == groupsInThisCond(iGroup));

                    % Calculate or extract standard error and mean
                    if ~isempty(means) && ~isempty(errors)
                        groupStat = means(actualGroupIdx, iCond);
                        groupSem = errors(actualGroupIdx, iCond);
                    else
                        [groupSem, groupStat] = nanstderr(groupData, 'DataType', dataType, 'Weights', groupWeights);
                    end
                    
                    % Retrieve the color for the current group
                    groupColor = colorMap{actualGroupIdx};

                    % Plot error bars
                    if plotErrorBars
                        errorbar(xMeanPositions(iGroup), groupStat, groupSem, ...
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
            % Process and plot statistics pooled across groups
            if ~isempty(dataThisCond)
                % Calculate or extract standard error and mean
                if ~isempty(means) && ~isempty(errors)
                    pooledStat = means(iCond);
                    pooledSem = errors(iCond);
                else
                    [pooledSem, pooledStat] = nanstderr(dataThisCond, 'DataType', dataType, 'Weights', weightsThisCond);
                end

                % Plot error bars
                if plotErrorBars
                    errorbar(currentCondValue, pooledStat, pooledSem, ...
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

% Run statistical tests across conditions
if (runTTest || runRankTest) && nConditions >= 2
    % Get the values for the first two conditions to compare
    cond1Value = uniqueConditionValues(1:end-1);
    cond2Value = uniqueConditionValues(2:end);

    % Get current y-axis limits to position text
    yLims = ylim(handles.ax);
    yRange = diff(yLims);
    yPos = yLims(2); 

    % Define x-position for the text
    xPosText = mean([cond1Value, cond2Value]);
    
    % Execute tests separated by group
    if statisticsByGroup
        % Loop through each group and test separately
        for iGroup = 1:nGroups
            % Extract the value for the current group
            currentGroupValue = uniqueGroupValues(iGroup);
            
            % Get data for this group in condition 1 & 2
            group1Data = data((grouping == currentGroupValue) & (condition == cond1Value));
            group2Data = data((grouping == currentGroupValue) & (condition == cond2Value));

            % Skip testing if data is missing for this group
            if isempty(group1Data) || isempty(group2Data)
                continue; 
            end

            % Check if data is normal to select the appropriate test
            isNormal1 = test_normality(group1Data);
            isNormal2 = test_normality(group2Data);
            isAppropriateForTTest = isNormal1 && isNormal2;

            % Get the color for this group
            groupColor = colorMap{iGroup};

            % Run and plot t-test if requested
            if runTTest
                [~, p_t] = ttest2(group1Data, group2Data);
                yPos = yPos - statTextYPosStepSize * yRange; 
                yRel = (yPos - yLims(1)) / yRange;
                
                hT = plot_test_result(p_t, 'PString', 'p_t', ...
                            'XLocText', xPosText, 'XLocStar', xPosText, ...
                            'YLocTextRel', yRel, 'YLocStarRel', yRel + statStarYLocOffset, ...
                            'SigLevel', sigLevel, 'IsAppropriate', isAppropriateForTTest);

                % Override color for texts
                set(hT.pText, 'Color', groupColor); 
            end

            % Run and plot rank-sum test if requested
            if runRankTest
                p_r = ranksum(group1Data, group2Data);
                yPos = yPos - statTextYPosStepSize * yRange; 
                yRel = (yPos - yLims(1)) / yRange;

                hR = plot_test_result(p_r, 'PString', 'p_r', ...
                            'XLocText', xPosText, 'XLocStar', xPosText, ...
                            'YLocTextRel', yRel, 'YLocStarRel', yRel + statStarYLocOffset, ...
                            'SigLevel', sigLevel, 'IsAppropriate', ~isAppropriateForTTest);

                % Override color for texts
                set(hR.pText, 'Color', groupColor); 
            end
        end
    else
        % Pool groups and test between conditions
        cond1Data = data(condition == cond1Value);
        cond2Data = data(condition == cond2Value);
        
        % Execute tests if both conditions have data
        if ~isempty(cond1Data) && ~isempty(cond2Data)
            % Check if pooled data is normal
            isNormal1 = test_normality(cond1Data);
            isNormal2 = test_normality(cond2Data);
            isAppropriateForTTest = isNormal1 && isNormal2;
            
            % Set color for pooled statistics text
            statsColor = pooledStatColor; 

            % Run and plot t-test if requested
            if runTTest
                [~, p_t] = ttest2(cond1Data, cond2Data);
                yPos = yPos - statTextYPosStepSize * yRange; 
                yRel = (yPos - yLims(1)) / yRange;
                
                hT = plot_test_result(p_t, 'PString', 'p_t', ...
                            'XLocText', xPosText, 'XLocStar', xPosText, ...
                            'YLocTextRel', yRel, 'YLocStarRel', yRel + statStarYLocOffset, ...
                            'SigLevel', sigLevel, 'IsAppropriate', isAppropriateForTTest);
                set(hT.pText, 'Color', statsColor); 
            end

            % Run and plot rank-sum test if requested
            if runRankTest
                p_r = ranksum(cond1Data, cond2Data);
                yPos = yPos - statTextYPosStepSize * yRange; 
                yRel = (yPos - yLims(1)) / yRange;

                hR = plot_test_result(p_r, 'PString', 'p_r', ...
                            'XLocText', xPosText, 'XLocStar', xPosText, ...
                            'YLocTextRel', yRel, 'YLocStarRel', yRel + statStarYLocOffset, ...
                            'SigLevel', sigLevel, 'IsAppropriate', ~isAppropriateForTTest);
                set(hR.pText, 'Color', statsColor); 
            end
        end
    end
end

% Hold off
hold_off(wasHold, axHandle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%