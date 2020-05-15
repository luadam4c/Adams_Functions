function handles = plot_grouped_scatter (xValues, yValues, varargin)
%% Plot a grouped scatter plot with 95% confidence ellipses
% Usage: handles = plot_grouped_scatter (xValues, yValues, grouping (opt), varargin)
% Explanation:
%       TODO
%
%
% Example(s):
%       xVec = randi(10, 10, 1);
%       yVec = randi(10, 10, 1) + 10;
%       grouping = [ones(5, 1); zeros(5, 1)];
%       plot_grouped_scatter(xVec, yVec)
%       plot_grouped_scatter(xVec, yVec, grouping)
%
% Outputs:
%       handles     - handles to TODO
%                   specified as a structure
%
% Arguments:
%       xValues     - xValues
%                   must be an array
%       yValues     - yValues
%                   must be an array
%       grouping    - (opt) group assignment for each data point
%                   must be an array
%                   default == the column number for a 2D array
%       varargin    - 'PlotOnly': whether to plot the scatters only
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotEllipse': whether to plot 95% condifence ellipses
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%       TODO
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == 'suppress'
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == 'suppress'
%                   - 'XUnits': x-axis units
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == 'unit'
%                   - 'XLabel': label for the x axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == strcat('XValue (', xUnits, ')')
%                   - 'YUnits': y-axis units
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == 'unit'
%                   - 'YLabel': label for the y axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == strcat('YValue (', yUnits, ')')
%                   - 'GroupingLabels': labels for the groupings, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == {'Group #1', 'Group #2', ...}
%                   - 'ColorMap': a color map that provides a different
%                                   color for each group
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
%                   - 'FigTitle': title for the figure
%                   must be a string scalar or a character vector
%                   default == strcat(yLabel, ' vs. ', xLabel)
%                   - 'FigHandle': figure handle for created figure
%                   must be a empty or a figure object handle
%                   default == []
%                   - 'FigNumber': figure number for creating figure
%                   must be a positive integer scalar
%                   default == []

%                   - 'OutFolder': directory to save figure, 
%                                   e.dots. 'output'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.dots., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - Any other parameter-value pair for the gscatter() function
%
% Requires:
%       cd/compute_axis_limits.m
%       cd/compute_confidence_ellipse.m
%       cd/construct_fullpath.m
%       cd/create_default_grouping.m
%       cd/hold_off.m
%       cd/hold_on.m
%       cd/isemptycell.m
%       cd/islegendlocation.m
%       cd/ispositiveintegerscalar.m
%       cd/isfigtype.m
%       cd/save_all_figtypes.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/m3ha_simulate_population.m
%

% File History:
% 2017-12-13 Modified from plot_grouped_histogram.m
% 2018-05-27 Fixed the case when nGroups is NaN or 0
% 2020-05-14 Added many optional arguments
% TODO: Merge with plot_correlation

%% Hard-coded parameters
maxInFigure = 8;                % maximum number of groups to keep the legend
                                %   inside the figure

%% TODO: make the following optional arguments with given default
ellipseNPoints = 1000;              % 1000 points
ellipseLineStyle = '-';             % dashed line
ellipseLineWidth = 1;

%% Default values for optional arguments
groupingDefault = [];           % set later
plotOnlyDefault = false;            % setup default labels by default
plotEllipseDefault = true;          % whether to plot ellipses by default
confidenceLevelDefault = 95;        % default confidence level (%)
xScaleDefault = 'linear';
yScaleDefault = 'linear';
xLimitsDefault = [];
yLimitsDefault = [];
xUnitsDefault = 'unit';             % the default x-axis units
xLabelDefault = '';                 % set later
yUnitsDefault = 'unit';             % the default y-axis units
yLabelDefault = '';                 % set later
groupingLabelsDefault = '';         % set later
colorMapDefault = [];               % set later
legendLocationDefault = 'auto';     % set later
figTitleDefault = '';               % set later
figHandleDefault = [];              % no existing figure by default
figNumberDefault = [];              % no figure number by default
outFolderDefault = '';              % default directory to save figure
figNameDefault = '';                % don't save figure by default
figTypesDefault = 'png';            % save as png file by default
markerSizeDefault = [];
markerTypeDefault = 'o';
markerLineWidthDefault = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'xValues');
addRequired(iP, 'yValues');

% Add optional inputs to the Input Parser
addOptional(iP, 'grouping', groupingDefault);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PlotOnly', plotOnlyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotEllipse', plotEllipseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ConfidenceLevel', confidenceLevelDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'XScale', xScaleDefault, ...
    @(x) any(validatestring(x, {'linear', 'log'})));
addParameter(iP, 'YScale', yScaleDefault, ...
    @(x) any(validatestring(x, {'linear', 'log'})));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'XUnits', xUnitsDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'XLabel', xLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'YUnits', yUnitsDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'GroupingLabels', groupingLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ColorMap', colorMapDefault);
addParameter(iP, 'LegendLocation', legendLocationDefault, ...
    @(x) all(islegendlocation(x, 'ValidateMode', true)));
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigHandle', figHandleDefault);
addParameter(iP, 'FigNumber', figNumberDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'FigNumber must be a empty or a positive integer scalar!'));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));
addParameter(iP, 'MarkerSize', markerSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'MarkerType', markerTypeDefault, ...
    @(x) any(validatestring(x, {'o', '+', '*', '.', 'x', 's', 'd', ...
                                '^', 'v', '>', '<', 'p', 'fig'})));
addParameter(iP, 'MarkerLineWidth', markerLineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));

% Read from the Input Parser
parse(iP, xValues, yValues, varargin{:});
grouping = iP.Results.grouping;
plotOnly = iP.Results.PlotOnly;
plotEllipse = iP.Results.PlotEllipse;
confidenceLevel = iP.Results.ConfidenceLevel;
xScale = iP.Results.XScale;
yScale = iP.Results.YScale;
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;
xUnits = iP.Results.XUnits;
xLabel = iP.Results.XLabel;
yUnits = iP.Results.YUnits;
yLabel = iP.Results.YLabel;
groupingLabels = iP.Results.GroupingLabels;
colorMap = iP.Results.ColorMap;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);
figTitle = iP.Results.FigTitle;
figHandle = iP.Results.FigHandle;
figNumber = iP.Results.FigNumber;
outFolder = iP.Results.OutFolder;
figName = iP.Results.FigName;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
markerSize = iP.Results.MarkerSize;
markerType = iP.Results.MarkerType;
markerLineWidth = iP.Results.MarkerLineWidth;

% Keep unmatched arguments for the gscatter() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Change default if plotting histogram only
if plotOnly
    xLimits = 'suppress';
    yLimits = 'suppress';
    xLabel = 'suppress';
    yLabel = 'suppress';
    figTitle = 'suppress';
    legendLocation = 'suppress';
    plotEllipse = false;
end

% Set the default x-axis label
if isempty(xLabel)
    xLabel = strcat('XValue (', xUnits, ')');
end

% Set the default y-axis label
if isempty(yLabel)
    yLabel = strcat('YValue (', yUnits, ')');
end

% Set the default figure title
if isempty(figTitle)
    figTitle = replace([yLabel, ' vs. ', xLabel], '_', '\_');
end

% If the figure name is not a full path, create full path
if ~isempty(figName)
    figName = construct_fullpath(figName, 'Directory', outFolder);
end

% Decide on the grouping vector from the xValues
[groupingFromX, ~, ~, xValues] = ...
    create_default_grouping('Stats', xValues, 'Grouping', grouping, ...
                            'GroupingLabels', groupingLabels, ...
                            'ToLinearize', true);

% Decide on the grouping vector and possibly labels from the yValues
[groupingFromY, groupValues, groupingLabels, yValues] = ...
    create_default_grouping('Stats', yValues, 'Grouping', grouping, ...
                            'GroupingLabels', groupingLabels, ...
                            'ToLinearize', true);

% Make sure the grouping vectors match
if ~isequal(groupingFromX, groupingFromY)
    error('xValues and yValues don''t match!');
else
    grouping = groupingFromY;
end

% Create ellipse labels
if plotEllipse
    ellipseLabels = strcat(groupingLabels, 'Ellipse');
end

% Linearize y values as well
yValues = yValues(:);

% Count the number of groups
nGroups = numel(groupValues);

% Don't plot any if there are no groups
if nGroups == 0
    handles.fig = gobjects(1);
    handles.dots = gobjects(1);
    return
end

% Decide on the color map
colorMap = decide_on_colormap(colorMap, nGroups);
colorMapCell = decide_on_colormap(colorMap, nGroups, 'ForceCellOutput', true);

% Set legend location based on number of groups
if strcmpi(legendLocation, 'auto')
    if nGroups > 1 && nGroups <= maxInFigure
        legendLocation = 'northeast';
    elseif nGroups > maxInFigure
        legendLocation = 'eastoutside';
    else
        legendLocation = 'suppress';
    end
end

%% Compute confidence ellipses
if plotEllipse
    % Compute all confidence ellipses
    ellipseTable = compute_confidence_ellipse(xValues, yValues, grouping, ...
                            xScale, yScale, confidenceLevel, ellipseNPoints);

    % Extract from table
    xEllipses = ellipseTable.xEllipses;
    yEllipses = ellipseTable.yEllipses;

    % Compute default axis limits
    if isempty(xLimits)
        xLimits = compute_axis_limits({xValues, xEllipses}, ...
                                        'x', 'Coverage', 90);
    end
    if isempty(yLimits)
        yLimits = compute_axis_limits({yValues, yEllipses}, ...
                                        'y', 'Coverage', 90);
    end
end

%% Plot and save scatter plot
% Decide on the figure to plot on
fig = set_figure_properties('FigHandle', figHandle, 'FigNumber', figNumber);

% Store hold status and hold on
wasHold = hold_on;

% Plot grouped scatter plot
dots = gscatter(xValues, yValues, grouping, ...
                colorMap, markerType, otherArguments{:});

% Set the legend labels for each Line object
arrayfun(@(x) set(dots(x), 'DisplayName', groupingLabels{x}), 1:nGroups);

% Update marker size
if ~isempty(markerSize)
    set(dots, 'MarkerSize', markerSize);
end

% Update marker line width
if ~isempty(markerLineWidth)
    set(dots, 'LineWidth', markerLineWidth);
end

% Update x and y axis scales
set(gca, 'XScale', xScale, 'YScale', yScale);

% Set x axis limits
if ~isempty(xLimits) && xLimits(1) ~= xLimits(2) && ...
        ~strcmpi(xLimits, 'suppress')
    xlim(xLimits);
end

% Set y axis limits
if ~isempty(yLimits) && yLimits(1) ~= yLimits(2) && ...
        ~strcmpi(yLimits, 'suppress')
    ylim(yLimits);
end

% Generate a legend if there is more than one group
if ~strcmpi(legendLocation, 'suppress')
    lgd = legend('Location', legendLocation);
    set(lgd, 'AutoUpdate', 'off', 'Interpreter', 'none');
end

% Plot a 95% confidence ellipse for each group
if plotEllipse
    % Decide on the ellipses to plot
    toPlot = ~isemptycell(xEllipses) & ~isemptycell(yEllipses);

    % Print the groups without ellipses
    arrayfun(@(x) fprintf('No ellipse will be plotted for Group #%d!\n', x), ...
            find(~toPlot));

    % Plot ellipses
    ellipses = ...
        cellfun(@(x, y, c, d) plot(x, y, 'Color', c, ...
                                'LineStyle', ellipseLineStyle, ...
                                'LineWidth', ellipseLineWidth, ...
                                'DisplayName', d), ...
                xEllipses(toPlot), yEllipses(toPlot), ...
                colorMapCell(toPlot), ellipseLabels(toPlot));
end

% Generate an x-axis label
if ~strcmpi(xLabel, 'suppress')
    xlabel(xLabel);
end

% Generate a y-axis label
if ~strcmpi(yLabel, 'suppress')
    ylabel(yLabel);
end

% Generate a title
if ~strcmpi(figTitle, 'suppress')
    title(figTitle);
end

% Hold off if it was originally so
hold_off(wasHold);

% Save the figure if requested
if ~isempty(figName)
    % Save the figure in all file types requested
    save_all_figtypes(fig, figName, figTypes);

    % Close figure
    close(fig);
end

%% Save in output
handles.fig = fig;
handles.dots = dots;
handles.ellipses = ellipses;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ellipseTable = compute_confidence_ellipse (xValues, yValues, ...
                    grouping, xScale, yScale, confidenceLevel, ellipseNPoints)
%% Computes confidence ellipses from data
%% TODO: Pull out as its own function
% TODO: Make xScale, yScale, confidenceLevel, ellipseNPoints optional arguments
%   Requires:
%       cd/array_fun.m
%       cd/plot_ellipse.m
%       cd/unique_custom.m

%% Preparation
% Compute the critical value for the chi-squared distribution
%   to have a cumulative probability of the confidence level
criticalValue = chi2inv(confidenceLevel/100, 2);

% Get all unique grouping values
uniqueGroupValues = unique_custom(grouping, 'IgnoreNaN', true);

%% Fit each group with a bivariate Gaussian
% Use the sample means and sample covariance matrices, 
%   the maximum likelihood estimates of the corresponding parameters 
%   in a bivariate Gaussian fit, to compute the center, 
%   half lengths, and angles of the confidenceLevel % confidence ellipse
[sampleMeans, sampleCovariances, eigenValues, eigenVectors, ...
        halfLengths, angles, xEllipses, yEllipses] = ...
    array_fun(@(groupValue) compute_confidence_ellipse_helper(...
                            xValues(grouping == groupValue, :), ...
                            yValues(grouping == groupValue, :), ...
                            xScale, yScale, criticalValue, ellipseNPoints), ...
            uniqueGroupValues, 'UniformOutput', false);

%% Return outputs in a table
ellipseTable = table(sampleMeans, sampleCovariances, ...
                    eigenValues, eigenVectors, ...
                    halfLengths, angles, xEllipses, yEllipses);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sampleMeans, sampleCovariances, eigenValues, eigenVectors, ...
            halfLengths, angles, xEllipses, yEllipses] = ...
                compute_confidence_ellipse_helper (xValues, yValues, ...
                                xScale, yScale, criticalValue, ellipseNPoints)
% Computes a confidence ellipse from data

% Transform data to the scale used for plotting
if strcmp(xScale, 'log')
    xScaled = log10(xValues);
else
    xScaled = xValues;
end
if strcmp(yScale, 'log')
    yScaled = log10(yValues);
else
    yScaled = yValues;
end

% Transform data to column vectors
xScaled = xScaled(:);
yScaled = yScaled(:);

% Place data in two columns
twoColumns = [xScaled, yScaled];

% Compute the sample mean for non-NaN data
sampleMeans = nanmean(twoColumns);

% Compute the sample covariance for non-NaN data
sampleCovariances = nancov(twoColumns);

% Compute the eigenvalues of the covariance matrix
eigenValues = eig(sampleCovariances);

if any(eigenValues<= 0)
    fprintf('There seems to be a direct correlation!\n');
    fprintf('Covariances: \n');
    disp(sampleCovariances);
    fprintf('Eigenvalues: \n');
    disp(eigenValues);
    fprintf('\n');
end

% Find the half lengths and angle of rotation for the ellipse
%   Note: if an eigenvalue is zero, there is direct correlation and
%           no ellipse will be plotted
if length(eigenValues) >= 2 && all(eigenValues> 0)
    % Compute half lengths for the 95% confidence ellipse
    halfLengths = sqrt(eigenValues.* criticalValue);

    % Compute the eigenvectors of the covariance
    [eigenVectors, ~] = eig(sampleCovariances);

    % Compute the angle of rotation in radians
    %   This is the angle between the x axis and the first eigenvector
    angles = atan(eigenVectors(2, 1) / eigenVectors(1, 1));
end

% Find the x and y values for the ellipse
if length(sampleMeans) >= 2 && ~isempty(halfLengths) && ...
        ~isempty(angles)
    % Obtain the x and y values of the ellipse on the scaled plot
    [~, xPlot, yPlot] = ...
        plot_ellipse(sampleMeans, halfLengths, ...
                    angles, 'NPoints', ellipseNPoints, ...
                    'ToPlot', false);

    % Convert back to original scale
    if strcmp(xScale, 'log')
        xEllipses = 10 .^ (xPlot);
    else
        xEllipses = xPlot;
    end
    if strcmp(yScale, 'log')
        yEllipses = 10 .^ (yPlot);
    else
        yEllipses = yPlot;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
