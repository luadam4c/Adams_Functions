function handles = plot_grouped_jitter (data, varargin)
%% Plots a violin (unpaired comparison) plot from data
% Usage: handles = plot_grouped_jitter (data, varargin)
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
%       varargin    - 'RelativeBandwidth': band width relative to the range
%                   must be a numeric scalar
%                   default == 0.1
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == [0.5, nConditions + 0.5]
%                   - 'XTickLabels': x tick labels in place of parameter values
%                   must be a cell array of character vectors/strings
%                   default == {}
%                   - 'XTickAngle': angle for parameter tick labels
%                   must be a numeric scalar
%                   default == 0
%                   - 'YLabel': label for the y axis, 
%                   must be a string scalar or a character vector 
%                   default == none
%                   - 'ColorMap': a color map for each group
%                   must be a numeric array with 3 columns
%                   default == set in decide_on_colormap.m
%                   - Any other parameter-value pair for violinplot()
%
% Requires:
%       cd/addpath_custom.m
%       cd/apply_iteratively.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/decide_on_colormap.m
%       cd/extract_columns.m
%       cd/hold_on.m
%       cd/hold_off.m
%       cd/locate_functionsdir.m
%       cd/struct2arglist.m
%       ~/Downloaded_Functions/Violinplot/plotSpread.m
%
% Used by:
%       cd/m3ha_plot_jitter.m

% File History:
% 2020-04-13 Modified from plot_violin.m
% 

%% Hard-coded parameters
% TODO: Make optional argument
groupingLabels = {};

%% Default values for optional arguments
relativeBandWidthDefault = 0.1;
xLimitsDefault = [];            % set later
xTickLabelsDefault = {};        % set later
xTickAngleDefault = [];         % set later
yLabelDefault = '';             % no y label by default
colorMapDefault = [];           % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% If not compiled, add directories to search path for required functions
if ~isdeployed
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

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'RelativeBandWidth', relativeBandWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'XTickLabels', xTickLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'XTickAngle', xTickAngleDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ColorMap', colorMapDefault);

% Read from the Input Parser
parse(iP, data, varargin{:});
relativeBandWidth = iP.Results.RelativeBandWidth;
xLimits = iP.Results.XLimits;
xTickLabels = iP.Results.XTickLabels;
xTickAngle = iP.Results.XTickAngle;
yLabel = iP.Results.YLabel;
colorMap = iP.Results.ColorMap;

% Keep unmatched arguments for the violinplot() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Force data values as a cell array
valuesEachCondition = force_column_cell(data);

% Count the number of conditions
nConditions = numel(valuesEachCondition);

% Count the number of groups
nGroups = numel(valuesEachCondition{1});

% Create group numbers
groupNumbers = transpose(1:nGroups);

% Decide on the x tick labels
if isempty(xTickLabels)
    xTickLabels = create_labels_from_numbers(1:nConditions, ...
                                                'Prefix', 'Condition');
end

% Decide on the grouping labels
if isempty(groupingLabels)
    groupingLabels = create_labels_from_numbers(1:nGroups, 'Prefix', 'Group');
end

% Decide on x axis limits
if isempty(xLimits)
    xLimits = [0.5, nConditions + 0.5];
end

% Decide on the color map
colorMap = decide_on_colormap(colorMap, nGroups);

% Compute range of all values
rangeValues = apply_iteratively(@max, valuesEachCondition) - ...
                apply_iteratively(@min, valuesEachCondition);

% Reorganize cell array
valuesEachGroup = extract_columns(valuesEachCondition, 'TreatCnvAsColumn', true, ...
                                    'OutputMode', 'single');

%% Do the job
% Return if there is no data
if isempty(valuesEachCondition)
    handles = struct;
    return
end

% Hold on
wasHold = hold_on;

% Plot the data points for each cell
handles = cellfun(@(a, b) plotSpread(a, ...
                                'distributionColors', colorMap(b, :), ...
                                otherArguments{:}), ...
                valuesEachGroup, num2cell(groupNumbers), 'UniformOutput', false);

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

% Hold off
hold_off(wasHold);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%