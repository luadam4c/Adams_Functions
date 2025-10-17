function handles = plot_grouped_jitter (data, varargin)
%% Plots a jitter plot colored by group from data (uses plotSpread)
% Usage: handles = plot_grouped_jitter (data, grouping (opt), condition (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
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
%       varargin    - 'UsePlotSpread': whether to use plotSpread.m
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
%                   - 'YLabel': label for the y axis, 
%                   must be a string scalar or a character vector 
%                   default == none
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
% TODO: Implement cell array input

%% Hard-coded parameters
maxNGroupsForInnerLegend = 8;
maxNGroupsForOuterLegends = 25;
forceVectorInput = true;       % Consider making this into an optional argument
defaultJitterWidthNotPlotSpread = 0.3;
markerDefault = '.';           % TODO: Make this an optional argument
meanMarker = 'o';              % TODO: Make this an optional argument
meanMarkerSize = 10;           % TODO: Make this an optional argument
meanLineWidth = 1.5;           % TODO: Make this an optional argument
errorBarCapSize = 0;           % TODO: Make this an optional argument
sigLevel = 0.05;               % Significance level for tests

%% Default values for optional arguments
groupingDefault = [];           % set later
distributionDefault = [];       % set later
usePlotSpreadDefault = [];      % set later
jitterWidthDefault = [];        % set later
plotMeanValuesDefault = true;
plotErrorBarsDefault = true;
runTTestDefault = true;
runRankTestDefault = true;
xLimitsDefault = [];            % set later
xTickLocsDefault = [];           % set later
xTickLabelsDefault = {};        % set later
groupingLabelsDefault = '';     % set later
xTickAngleDefault = [];         % set later
yLabelDefault = '';             % no y label by default
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
addParameter(iP, 'UsePlotSpread', usePlotSpreadDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'JitterWidth', jitterWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'PlotMeanValues', plotMeanValuesDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotErrorBars', plotErrorBarsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RunTTest', runTTestDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RunRankTest', runRankTestDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'XTickLocs', xTickLocsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'XTickLabels', xTickLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'GroupingLabels', groupingLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'XTickAngle', xTickAngleDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ColorMap', colorMapDefault);
addParameter(iP, 'LegendLocation', legendLocationDefault, ...
    @(x) all(islegendlocation(x, 'ValidateMode', true)));
addParameter(iP, 'AxesHandle', axHandleDefault);

% Read from the Input Parser
parse(iP, data, varargin{:});
grouping = iP.Results.grouping;
condition = iP.Results.condition;
usePlotSpread = iP.Results.UsePlotSpread;
jitterWidth = iP.Results.JitterWidth;
plotMeanValues = iP.Results.PlotMeanValues;
plotErrorBars = iP.Results.PlotErrorBars;
runTTest = iP.Results.RunTTest;
runRankTest = iP.Results.RunRankTest;
xTickLocs = iP.Results.XTickLocs;
xLimits = iP.Results.XLimits;
xTickLabels = iP.Results.XTickLabels;
groupingLabels = iP.Results.GroupingLabels;
xTickAngle = iP.Results.XTickAngle;
yLabel = iP.Results.YLabel;
colorMap = iP.Results.ColorMap;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);
axHandle = iP.Results.AxesHandle;

% Keep unmatched arguments for the plotSpread() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Decide on whether to use plotSpread.m
if isempty(usePlotSpread)
    if exist('plotSpread.m', 'file') == 2
        usePlotSpread = true;
    else
        usePlotSpread = false;
    end
end

%% If not compiled, add directories to search path for required functions
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


% Decide on the grouping vector and possibly labels if the data is a 
%   matrix or cell array
% If no grouping vector, each column is a group
[grouping, uniqueGroupValues, groupingLabels, data] = ...
    create_default_grouping('Stats', data, 'Grouping', grouping, ...
                            'GroupingLabels', groupingLabels, ...
                            'GroupingLabelPrefix', 'Group');

% Decide on the condition vector and possibly labels if the data is a
%   matrix or cell array
% If no condition vector, each column is a condition
[condition, uniqueConditionValues, xTickLabels, data] = ...
    create_default_grouping('Stats', data, 'Grouping', condition, ...
                            'GroupingLabels', xTickLabels, ...
                            'GroupingLabelPrefix', 'Condition');

% Concatenate everything into a single column vector
if forceVectorInput
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

% Count the number of unique conditions
nConditions = numel(uniqueConditionValues);

% Count the number of groups
nGroups = numel(uniqueGroupValues);

% Count the number of data points
nPoints = numel(data);

% Decide on the x tick locations
if isempty(xTickLocs) || ischar(xTickLocs) && strcmpi(xTickLocs, 'suppress')
    % Define x tick locations based on the unique condition indices
    xTickLocs = uniqueConditionValues;
end

% Decide whether to update the x tick locations
if ischar(xTickLocs) && strcmpi(xTickLocs, 'suppress')
    toUpdateXTicks = false;
else
    toUpdateXTicks = true;
end

% Decide on the color map, using the lines map by default
if isempty(colorMap)
    colorMap = @lines;
end
colorMap = decide_on_colormap(colorMap, nGroups, 'ForceCellOutput', true);

% Set legend location based on number of subplots
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
% Return if there is no data
if isempty(data)
    handles = struct;
    return
end

% Hold on
wasHold = hold_on(axHandle);

% Plot the data points
if usePlotSpread
    % Don't show plotSpread's means if we are plotting our own
    if plotMeanValues || plotErrorBars
        otherArguments = [otherArguments, {'showMM', 0}];
    end

    output = plotSpread(axHandle, data, 'distributionIdx', condition, ...
                            'categoryIdx', grouping, ...
                            'categoryLabels', groupingLabels, ...
                            'categoryColors', colorMap, ...
                            'spreadWidth', jitterWidth, ...
                            otherArguments{:});

    % Reformat handles outputs
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
    
    % Use the built-in plot() function to plot each group with 
    %   marker set by markerDefault and the desired color map
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
            plot(axHandle, xCoords, yCoords, markerDefault, ...
                    'Color', groupColor, 'DisplayName', groupLabel, ...
                    'LineStyle', 'none', otherArguments{:});
    end

    % Export handles
    handles.ax = axHandle;
    handles.distributions = distributions;
    handles.stats = []; % No stats are calculated in this manual version
end

%% Plot statistics
% Plot means and error bars for each condition
if plotMeanValues || plotErrorBars
    for iCond = 1:nConditions
        currentCondValue = uniqueConditionValues(iCond);
        isCurrentCond = (condition == currentCondValue);
        dataThisCond = data(isCurrentCond);
        groupingThisCond = grouping(isCurrentCond);
        groupsInThisCond = unique(groupingThisCond);
        nGroupsInThisCond = numel(groupsInThisCond);
        
        if nGroupsInThisCond > 0
            dataByGroup = arrayfun(@(x) dataThisCond(groupingThisCond == x), ...
                                    groupsInThisCond, 'UniformOutput', false);
            spreadWidth = 0.4;
            if nGroupsInThisCond == 1
                xMeanPositions = currentCondValue;
            else
                offsets = linspace(-spreadWidth/2, spreadWidth/2, nGroupsInThisCond);
                xMeanPositions = currentCondValue + offsets;
            end

            for iGroup = 1:nGroupsInThisCond
                groupData = dataByGroup{iGroup};
                groupMean = mean(groupData, 'omitnan');
                groupSem = std(groupData, 'omitnan') / sqrt(numel(groupData));
                groupColor = colorMap{uniqueGroupValues == groupsInThisCond(iGroup)};

                if plotErrorBars
                    errorbar(xMeanPositions(iGroup), groupMean, groupSem, ...
                             'Color', groupColor, 'LineWidth', meanLineWidth, ...
                             'CapSize', errorBarCapSize, 'HandleVisibility', 'off');
                end
                if plotMeanValues
                    plot(xMeanPositions(iGroup), groupMean, meanMarker, ...
                         'MarkerEdgeColor', groupColor, ...
                         'MarkerSize', meanMarkerSize, 'LineWidth', meanLineWidth, ...
                         'HandleVisibility', 'off');
                end
            end
        end
    end
end

% Run statistical tests across conditions, separated by group
if (runTTest || runRankTest) && nConditions >= 2
    % Get the values for the first two conditions to compare
    cond1Value = uniqueConditionValues(1:end-1);
    cond2Value = uniqueConditionValues(2:end);

    % Get current y-axis limits to position text
    yLims = ylim(handles.ax);
    yRange = diff(yLims);
    yPos = yLims(2); % Start at the top

    % Define x-position for the text (midway between the two conditions)
    xPosText = mean([cond1Value, cond2Value]);

    % Loop through each group (color)
    for iGroup = 1:nGroups
        currentGroupValue = uniqueGroupValues(iGroup);
        
        % Get data for this group in condition 1 & 2
        group1Data = data((grouping == currentGroupValue) & (condition == cond1Value));
        group2Data = data((grouping == currentGroupValue) & (condition == cond2Value));

        if isempty(group1Data) || isempty(group2Data)
            continue; % Skip if data is missing for this group in either condition
        end

        % Check if data is normal (required for t-test)
        isNormal1 = test_normality(group1Data);
        isNormal2 = test_normality(group2Data);
        isAppropriateForTTest = isNormal1 && isNormal2;

        % Get the color for this group
        groupColor = colorMap{iGroup};

        if runTTest
            [~, p_t] = ttest2(group1Data, group2Data);
            yPos = yPos - 0.1 * yRange; % Move text down
            yRel = (yPos - yLims(1)) / yRange;
            
            % Plot t-test result
            hT = plot_test_result(p_t, 'PString', 'p_t', ...
                        'XLocText', xPosText, 'XLocStar', xPosText, ...
                        'YLocTextRel', yRel, 'YLocStarRel', yRel + 0.05, ...
                        'SigLevel', sigLevel, 'IsAppropriate', isAppropriateForTTest);
            set(hT.pText, 'Color', groupColor); % Override color for texts
        end

        if runRankTest
            p_r = ranksum(group1Data, group2Data);
            yPos = yPos - 0.1 * yRange; % Move text down
            yRel = (yPos - yLims(1)) / yRange;

            % Plot rank-sum test result
            hR = plot_test_result(p_r, 'PString', 'p_r', ...
                        'XLocText', xPosText, 'XLocStar', xPosText, ...
                        'YLocTextRel', yRel, 'YLocStarRel', yRel + 0.05, ...
                        'SigLevel', sigLevel, 'IsAppropriate', ~isAppropriateForTTest);
            set(hR.pText, 'Color', groupColor); % Override color for texts
        end
    end
end


%% Finalize plot
% Update x tick locations if desired
if toUpdateXTicks
    xticks(xTickLocs);
end

% Update x tick labels
if ~isempty(xTickLabels)
    xticklabels(xTickLabels);
end

% Modify x tick angle
if ~isempty(xTickAngle)
    xtickangle(xTickAngle);
end

% Decide on x axis limits based on x tick locations
if isempty(xLimits)
    xLimits = [min(xTickLocs) - 0.5, max(xTickLocs) + 0.5];
end

% Modify x limits
if ~(ischar(xLimits) && strcmpi(xLimits, 'suppress'))
    xlim(xLimits);
end

% Set y label
if ~isempty(yLabel)
    ylabel(yLabel);
end

% Generate a legend if there is more than one trace
if ~strcmpi(legendLocation, 'suppress')
    legend('location', legendLocation);
end

% Hold off
hold_off(wasHold, axHandle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Force data values as a matrix
valuesEachCondition = force_column_cell(data);

% Count the number of conditions
nConditions = numel(valuesEachCondition);

% Count the number of groups
nGroups = numel(valuesEachCondition{1});

% Create group numbers
groupNumbers = transpose(1:nGroups);

% Decide on the grouping labels
if isempty(groupingLabels)
    groupingLabels = create_labels_from_numbers(groupNumbers, 'Prefix', 'Group');
end

% Reorganize cell array
valuesEachGroup = ...
    extract_columns(valuesEachCondition, 'TreatCnvAsColumn', true, ...
                    'OutputMode', 'single');

handles = ...
    cellfun(@(a, b) plotSpread(a, 'distributionColors', colorMap(b, :), ...
                                otherArguments{:}), ...
            valuesEachGroup, num2cell(groupNumbers), 'UniformOutput', false);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
