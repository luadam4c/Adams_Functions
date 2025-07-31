function violins = plot_violin (data, varargin)
%% Plots a violin (unpaired comparison) plot from data
% Usage: violins = plot_violin (data, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       randVec1 = randi(10, 10, 1);
%       randVec2 = randi(10, 10, 1) + 10;
%       data = [randVec1, randVec2];
%       plot_violin(data)
%       plot_violin(data, 'MedianColorMap', 'GreenYellow')
%
% Outputs:
%       violins     - handles to plotted violins objects
%                   specified as a Violin object array
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
%                   default == [0.5, nGroups + 0.5]
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == uses compute_axis_limits.m
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
%                   - 'MedianColorMap': a color map for the median(s)
%                   must be a numeric array with 3 columns
%                   default == set in decide_on_colormap.m
%                   - Any other parameter-value pair for violinplot()
%
% Requires:
%       cd/addpath_custom.m
%       cd/apply_iteratively.m
%       cd/create_error_for_nargin.m
%       cd/decide_on_colormap.m
%       cd/locate_functionsdir.m
%       cd/struct2arglist.m
%       cd/force_data_as_matrix.m
%       ~/Downloaded_Functions/Violinplot/violinplot.m
%
% Used by:
%       cd/compare_events_pre_post_stim.m
%       cd/m3ha_plot_figure08.m
%       cd/m3ha_plot_violin.m
%       cd/m3ha_simulate_population.m

% File History:
% 2019-12-30 Moved from m3ha_plot_violin.m
% 2020-06-26 Added 'YLimits' as an optional argument

%% Hard-coded parameters
% Note: The following must be consistent with violinplot.m
% violinplotMedianEdgeColor = [0.5, 0.5, 0.5];

%% Default values for optional arguments
relativeBandWidthDefault = 0.1;
xLimitsDefault = [];            % set later
yLimitsDefault = [];            % set later
xTickLabelsDefault = {};        % set later
xTickAngleDefault = [];         % set later
yLabelDefault = '';             % no y label by default
colorMapDefault = [];           % set later
medianColorMapDefault = [];     % st later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% If not compiled, add directories to search path for required functions
if exist('violinplot.m', 'file') ~= 2 && ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for violinplot()
    addpath_custom(fullfile(functionsDirectory, ...
                            'Downloaded_Functions', 'Violinplot'));
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
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'XTickLabels', xTickLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'XTickAngle', xTickAngleDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ColorMap', colorMapDefault);
addParameter(iP, 'MedianColorMap', medianColorMapDefault);

% Read from the Input Parser
parse(iP, data, varargin{:});
relativeBandWidth = iP.Results.RelativeBandWidth;
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;
xTickLabels = iP.Results.XTickLabels;
xTickAngle = iP.Results.XTickAngle;
yLabel = iP.Results.YLabel;
colorMap = iP.Results.ColorMap;
medianColorMap = iP.Results.MedianColorMap;

% Keep unmatched arguments for the violinplot() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Force data values as a numeric matrix 
%   where each group is a column and each row is a sample
[dataValues, groupLabelsAuto] = force_data_as_matrix(data);

% Count the number of groups
nGroups = size(dataValues, 2);

% Decide on the x tick labels
if isempty(xTickLabels)
    xTickLabels = groupLabelsAuto;
end

% If there is any column with all NaN values, remove it from the violin plot
if any(all(isnan(dataValues), 1))
    % Test whether each column has no values
    noValues = all(isnan(dataValues), 1);

    % Remove those columns
    disp('Groups with all NaN values are removed!');
    dataValues = dataValues(:, ~noValues);
    xTickLabels = xTickLabels(~noValues);
    nGroups = size(dataValues, 2);
end

% Decide on x axis limits
if isempty(xLimits)
    xLimits = [0.5, nGroups + 0.5];
end

% Decide on y axis limits
if isempty(yLimits)
    yLimits = compute_axis_limits(dataValues, 'y');
end

% Decide on the color map
colorMap = decide_on_colormap(colorMap, nGroups);

% Compute range of all values
rangeValues = apply_iteratively(@max, dataValues) - ...
                apply_iteratively(@min, dataValues);

%% Do the job
% Return if there is no data
if isempty(dataValues)
    violins = Violin.empty;
    return
end

% Compute the bandwidth for the kernel density estimates
bandWidth = relativeBandWidth * rangeValues;

% Make sure the bandwidth is nonzero
if bandWidth == 0
    bandWidth = 1e-9;
end

% Plot a violin plot
% violinplot(dataValues, xTickLabels);
violins = violinplot(dataValues, xTickLabels, 'BandWidth', bandWidth, ...
                        otherArguments{:});
%{
while bandWidth > bandWidth * 10 ^ -2
    try
        violins = violinplot(dataValues, xTickLabels, 'BandWidth', bandWidth, ...
                                otherArguments{:});
        break
    catch
        relativeBandWidth = relativeBandWidth / 10;
        fprintf('Warning: Relative bandwidth changed to %g!\n', relativeBandWidth);
        bandWidth = bandWidth / 10;
        clf;
    end
end
%}
                        
% Apply the color map
for iGroup = 1:nGroups
    violins(iGroup).ViolinColor = colorMap(iGroup, :);
end

% Modify x limits
if ~(ischar(xLimits) && strcmpi(xLimits, 'suppress'))
    xlim(xLimits);
end

% Modify y limits
if ~(ischar(yLimits) && strcmpi(yLimits, 'suppress') || isempty(yLimits))
    ylim(yLimits);
end

% Modify x tick angle
if ~isempty(xTickAngle)
    xtickangle(xTickAngle);
end

% Set y label
if ~isempty(yLabel)
    ylabel(yLabel);
end

% Set median color if requested
if ~isempty(medianColorMap)
    medianColorMap = decide_on_colormap(medianColorMap, nGroups);

    % Find all median scatters using the edge color
%    medianScatters = findobj(gca, 'Type', 'Scatter', ...
%                            'MarkerEdgeColor', violinplotMedianEdgeColor);

    % Updated median scatter face color
    for iGroup = 1:nGroups
        if ishandle(violins(iGroup).MedianPlot)
            violins(iGroup).MedianColor = medianColorMap(iGroup, :);
        end
%        thisColor = medianColorMap(iGroup, :);
%        set(medianScatters(iGroup), 'MarkerFaceColor', medianColorMap(iGroup, :));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
