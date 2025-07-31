function handles = plot_grouped_jitter (data, varargin)
%% Plots a jitter plot colored by group from data (uses plotSpread)
% Usage: handles = plot_grouped_jitter (data, grouping (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%       randVec1 = randi(10, 10, 1);
%       randVec2 = randi(10, 10, 1) + 10;
%       data = [randVec1, randVec2];
%       plot_grouped_jitter(data)
%
%       data = [randn(50,1);randn(50,1)+3.5]*[1 1];
%       grouping = [[ones(50,1);zeros(50,1)],[randi([0,1],[100,1])]];
%       plot_grouped_jitter(data, grouping)
%
% Outputs:
%       handles     - handles to TODO
%                   specified as a structure
%
% Arguments:
%       data        - data table or data vectors
%                   Note: The dimension with fewer elements is taken as 
%                           the parameter
%                   must be a table or a numeric array
%                       or a cell array of numeric vectors
%       grouping    - (opt) group assignment for each data point
%                   must be an array of one the following types:
%                       'cell', 'string', numeric', 'logical', 
%                           'datetime', 'duration'
%                   default == the column number for a 2D array
%       varargin    - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == [0.5, nConditions + 0.5]
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
%                   - Any other parameter-value pair for violinplot()
%
% Requires:
%       cd/addpath_custom.m
%       cd/create_default_grouping.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/decide_on_colormap.m
%       cd/hold_on.m
%       cd/hold_off.m
%       cd/locate_functionsdir.m
%       cd/struct2arglist.m
%       ~/Downloaded_Functions/plotSpread/plotSpread.m
%
% Used by:
%       cd/m3ha_plot_grouped_jitter.m
%       cd/m3ha_plot_figure08.m
%       cd/m3ha_simulate_population.m

% File History:
% 2020-04-13 Modified from plot_violin.m
% 2020-04-18 Now uses the categoryidx option in plotSpread

%% Hard-coded parameters
maxNGroupsForInnerLegend = 8;
maxNGroupsForOuterLegends = 25;

%% Default values for optional arguments
groupingDefault = [];           % set later
xLimitsDefault = [];            % set later
xTickLabelsDefault = {};        % set later
groupingLabelsDefault = '';     % set later
xTickAngleDefault = [];         % set later
yLabelDefault = '';             % no y label by default
colorMapDefault = [];           % set later
legendLocationDefault = 'auto'; % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% If not compiled, add directories to search path for required functions
if exist('violinplot.m', 'file') ~= 2 && ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for violinplot()
    addpath_custom(fullfile(functionsDirectory, ...
                            'Downloaded_Functions', 'plotSpread'));
end

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

% Add parameter-value pairs to the Input Parser
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

% Read from the Input Parser
parse(iP, data, varargin{:});
grouping = iP.Results.grouping;
xLimits = iP.Results.XLimits;
xTickLabels = iP.Results.XTickLabels;
groupingLabels = iP.Results.GroupingLabels;
xTickAngle = iP.Results.XTickAngle;
yLabel = iP.Results.YLabel;
colorMap = iP.Results.ColorMap;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);

% Keep unmatched arguments for the violinplot() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Decide on the grouping vector and possibly labels
[grouping, uniqueGroupValues, groupingLabels, data] = ...
    create_default_grouping('Stats', data, 'Grouping', grouping, ...
                            'GroupingLabels', groupingLabels);

% Count the number of conditions
nConditions = size(data, 2);

% Count the number of groups
nGroups = numel(uniqueGroupValues);

% Decide on the x tick labels
if isempty(xTickLabels)
    xTickLabels = create_labels_from_numbers(1:nConditions, ...
                                                'Prefix', 'Condition');
end

% Decide on x axis limits
if isempty(xLimits)
    xLimits = [0.5, nConditions + 0.5];
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

%% Do the job
% Return if there is no data
if isempty(data)
    handles = struct;
    return
end

% Hold on
wasHold = hold_on;

% Plot the data points for each cell
output = plotSpread(data, 'categoryIdx', grouping, ...
                        'categoryLabels', groupingLabels, ...
                        'categoryColors', colorMap, otherArguments{:});

% Reformat handles outputs
distributions = output{1};
stats = output{2};
ax = output{3};
handles.distributions = distributions;
handles.stats = stats;
handles.ax = ax;

% Modify x limits
if ~(ischar(xLimits) && strcmpi(xLimits, 'suppress'))
    xlim(xLimits);
end

% Update x tick labels
if ~isempty(xTickLabels)
    xticklabels(xTickLabels);
end

% Modify x tick angle
if ~isempty(xTickAngle)
    xtickangle(xTickAngle);
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
hold_off(wasHold);

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
