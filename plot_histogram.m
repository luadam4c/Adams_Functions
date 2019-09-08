function [bars, fig] = plot_histogram (varargin)
%% Plots a histogram labelling out of range values differently
% Usage: [bars, fig] = plot_histogram (X (opt), varargin)
% Explanation:
%       Automatically combines the counts of X outside of the finite range 
%           of edges on the left or on the right to a bin on the left or 
%           on the right, respectively.
%       Note: The bar() function is used for the main histogram
%               unless 'UseBuiltIn' is set to true
%
% Example(s):
%       plot_histogram([-100, randn(1, 100) + 4, 100])
%       plot_histogram([randn(100, 3); [100, -200, 300]])
%       plot_histogram('Counts', (1:5)', 'Edges', (1:6)')
%
% Outputs:
%       bars        - handles to Bar objects
%                       bars(1:nGroups) - main histogram
%                       bars(nGroups+1) - left out of range bar if any
%                       bars(nGroups+2) - right out of range bar if any
%                   specified as a Bar object array
%       fig         - figure handle for the created figure
%                   specified as a figure object handle
% Side Effects:
%       Plots a histogram
% Arguments:
%       X           - (opt) data to distribute among bins
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%       varargin    - 'PlotOutliers': whether to plot outliers separately
%                   must be logical 1 (true) or 0 (false)
%                   default == true
%                   - 'UseBuiltIn': whether to use built in histogram() function
%                                   Note: this will not work if data is grouped
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Counts': bin counts, with each group 
%                                   being a different column
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%                   default == returned by compute_grouped_histcounts(stats)
%                   - 'Edges': bin edges
%                   must be a vector of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%                   default == returned by compute_grouped_histcounts(stats)
%                   - 'Grouping': group assignment for each data point
%                   must be an array of one the following types:
%                       'cell', 'string', numeric', 'logical', 
%                           'datetime', 'duration'
%                   default == the column number for a 2D array
%                   - 'GroupedStyle': histogram groupedStyle
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'side-by-side'  - same as 'grouped' bar graph
%                       'stacked'       - same as 'stacked' bar graph
%                       'overlapped'    - overlapped histograms
%                   default == 'stacked'
%                   - 'SpecialColor': color of expanded bins
%                   must be a 3-element numeric vector:
%                   default == [0 0.8 0.8] (light blue)
%                   - 'OutlierMethod': method for determining outliers
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'boxplot'   - same as the Matlab function boxplot()
%                       'isoutlier' - Use the built-in isoutlier function
%                       'fiveStds'  - Take out data points 
%                                       more than 5 standard deviations away
%                       'threeStds' - Take out data points 
%                                       more than 3 standard deviations away
%                       'twoStds'   - Take out data points 
%                                       more than 2 standard deviations away
%                   default == 'isoutlier'
%                   - 'XLimits': x-axis limits
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == minimum and maximum edges of bins
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == 'suppress'
%                   - 'XUnits': x-axis units
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == 'unit'
%                   - 'XLabel': label for the time axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == strcat('Value (', xUnits, ')')
%                   - 'YLabel': label(s) for the y axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector
%                   default == 'Count'
%                   - 'XTickLocs': locations of X ticks
%                   must be 'suppress' or a numeric vector
%                   default == TODO
%                   - 'XTickLabels': labels for each raster
%                   must be 'suppress' or a cell array of character/string arrays
%                   default == TODO
%                   - 'GroupingLabels': labels for the groupings, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == {'Group #1', 'Group #2', ...}
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
%                   default == strcat('Distribution of ', xLabel)
%                   - 'FigHandle': figure handle for created figure
%                   must be a empty or a figure object handle
%                   default == []
%                   - 'FigNumber': figure number for creating figure
%                   must be a positive integer scalar
%                   default == []
%                   - 'OutFolder': directory to save figure, 
%                                   e.g. 'output'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - Any other parameter-value pair for the bar() function
%
% Requires:
%       cd/argfun.m
%       cd/compute_centers_from_edges.m
%       cd/compute_grouped_histcounts.m
%       cd/construct_fullpath.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/extract_subvectors.m
%       cd/force_column_vector.m
%       cd/force_row_vector.m
%       cd/plot_grouped_histogram.m
%       cd/isfigtype.m
%       cd/isnumericvector.m
%       cd/remove_outliers.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/test_var_difference.m
%
% Used by:    
%       cd/plot_psth.m
%       /home/Matlab/Marks_Functions/paula/Oct2017/zgRasterFigureMaker.m
%       /media/adamX/m3ha/data_dclamp/initial_slopes.m

% File History:
% 2017-12-12 Created by Adam Lu
% 2018-06-05 Made edges an optional parameter and make the default dependent
%               on the isoutlier() and histcounts() functions
% 2018-06-11 Now uses the remove_outliers.m function
% 2019-01-15 Now uses plot_grouped_histogram.m by default
%               and added 'UseBuiltIn' as an optional parameter (default false)
% 2019-01-15 Made 'PlotOutliers' an optional parameter with default true
%       and rename as just plot_histogram.m
% 2019-02-24 Fixed bug when 'Counts' is passed in 
% 2019-03-14 Now returns empty plot if there is no data
% 2019-03-14 Fixed bug when xTickLabelNums is empty
% 2019-08-13 Fixed bug when histogram() is used
% 2019-09-08 Made X an optional argument (allow passing in of just 
%               'Counts' and 'Edges')
% 2019-09-08 Now uses set_figure_properties.m

%% Hard-coded parameters
validOutlierMethods = {'boxplot', 'isoutlier', ...
                        'fiveStds', 'threeStds', 'twoStds'};
validGroupedStyles = {'side-by-side', 'stacked', 'overlapped'};
minNXTicks = 5;

%% Default values for optional arguments
XDefault = [];                          % set later
plotOutliersDefault = true;             % plot outliers by default
useBuiltInDefault = false;              % use plot_grouped_histogram by default
countsDefault = [];                     % set later
edgesDefault = [];                      % set later
groupingDefault = [];                   % set later
groupedStyleDefault = 'stacked';        % grouped bars are stacked by default
specialColorDefault = [0, 0.8, 0.8];    % light blue
outlierMethodDefault = 'isoutlier';     % use built-in isoutlier function
xLimitsDefault = [];                    % set later
yLimitsDefault = 'suppress';            % set later
xUnitsDefault = 'unit';                 % the default x-axis units
xLabelDefault = '';                     % set later
yLabelDefault = 'Count';                % set later
xTickLocsDefault = [];                  % set later
xTickLabelsDefault = {};                % set later
groupingLabelsDefault = '';             % set later
legendLocationDefault = 'auto';         % set later
figTitleDefault = '';                   % set later
figHandleDefault = [];                  % no existing figure by default
figNumberDefault = [];                  % no figure number by default
outFolderDefault = '';                  % default directory to save figure
figNameDefault = '';                    % don't save figure by default
figTypesDefault = 'png';                % save as png file by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add optional inputs to the Input Parser
addOptional(iP, 'X', XDefault, ...      % data to distribute among bins
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Counts', countsDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'Edges', edgesDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'Grouping', groupingDefault, ...
    @(x) validateattributes(x, {'cell', 'string', 'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'GroupedStyle', groupedStyleDefault, ...
    @(x) any(validatestring(x, validGroupedStyles)));
addParameter(iP, 'PlotOutliers', plotOutliersDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'UseBuiltIn', useBuiltInDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SpecialColor', specialColorDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 3}));
addParameter(iP, 'OutlierMethod', outlierMethodDefault, ...
    @(x) any(validatestring(x, validOutlierMethods)));
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
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'XTickLocs', xTickLocsDefault, ...
    @(x) assert(ischar(x) && strcmpi(x, 'suppress') || isnumericvector(x), ...
        'XTickLocs must be ''suppress'' or a numeric vector!'));
addParameter(iP, 'XTickLabels', xTickLabelsDefault, ...
    @(x) assert(ischar(x) && strcmpi(x, 'suppress') || ...
                iscell(x) && all(cellfun(@(x) ischar(x) || isstring(x), x)), ...
        'XTickLabels must be ''suppress'' or a cell array of character/string arrays!'));
addParameter(iP, 'GroupingLabels', groupingLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
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

% Read from the Input Parser
parse(iP, varargin{:});
X = iP.Results.X;
plotOutliers = iP.Results.PlotOutliers;
useBuiltIn = iP.Results.UseBuiltIn;
counts = iP.Results.Counts;
edgesUser = iP.Results.Edges;
grouping = iP.Results.Grouping;
groupedStyle = validatestring(iP.Results.GroupedStyle, validGroupedStyles);
specialColor = iP.Results.SpecialColor;
outlierMethod = validatestring(iP.Results.OutlierMethod, validOutlierMethods);
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;
xUnits = iP.Results.XUnits;
xLabel = iP.Results.XLabel;
yLabel = iP.Results.YLabel;
xTickLocs = iP.Results.XTickLocs;
xTickLabels = iP.Results.XTickLabels;
groupingLabels = iP.Results.GroupingLabels;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);
figTitle = iP.Results.FigTitle;
figHandle = iP.Results.FigHandle;
figNumber = iP.Results.FigNumber;
outFolder = iP.Results.OutFolder;
figName = iP.Results.FigName;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the bar() or histogram() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Get the current MATLAB release
matlabRelease = version('-release');

% Get the current MATLAB release year
matlabYear = str2double(matlabRelease(1:4));

% Make sure there is data
if isempty(X) && (isempty(counts) || isempty(edgesUser))
    bars = gobjects;
    fig = gobjects;
    disp('There is no data to plot!');
    return
end

% Force rows as columns
X = force_column_vector(X, 'IgnoreNonVectors', true);
grouping = force_column_vector(grouping, 'TreatCellAsArray', true, ...
                                'IgnoreNonVectors', true);
edgesUser = force_column_vector(edgesUser, 'IgnoreNonVectors', true);

% If the figure name is passed in, make sure a new figure is created
if ~isempty(figName)
    % If the figure name is not a full path, create full path
    figName = construct_fullpath(figName, 'Directory', outFolder);

    % Create a new figure
    alwaysNew = true;
else
    alwaysNew = false;
end

% Decide on the figure to plot on
fig = set_figure_properties('FigHandle', figHandle, 'FigNumber', figNumber, ...
                            'AlwaysNew', alwaysNew);

%% Identify edges if not provided
if isempty(edgesUser)
    if plotOutliers
        % Find rows that contain outliers
        [~, rowsToKeep] = remove_outliers(X, 'OutlierMethod', outlierMethod);

        % Remove outliers if any
        XTrimmed = X(rowsToKeep, :);
        if ~isempty(grouping)
            groupingTrimmed = grouping(rowsToKeep);
        else
            groupingTrimmed = grouping;
        end

        % Find the proper bin edges for trimmed data
        [~, edgesUser] = compute_grouped_histcounts(XTrimmed, groupingTrimmed);
    else
        % Find the bin edges for all data
        [~, edgesUser] = compute_grouped_histcounts(X, grouping);
    end
else
    % Force as a column vector
    edgesUser = force_column_vector(edgesUser, 'IgnoreNonVectors', true);
end

%% Create histogram
% Extract finite part of edges
edgesFinite = edgesUser(isfinite(edgesUser));

% Expand edges for full histogram bin counts
edgesExpanded = [-Inf; edgesFinite; Inf];

% Compute histogram bincounts with expanded edges
%   Note: the first bin is always < the first non-Inf number in edges
%   Note: the last bin is always >= the last non-Inf number in edges
if isempty(counts)
    counts = compute_grouped_histcounts(X, grouping, 'Edges', edgesExpanded);
else
    % Force as a column vector
    counts = force_column_vector(counts, 'IgnoreNonVectors', true);

    % Create a row of zeros
    zeroRow = zeros(1, size(counts, 2));
    
    % Add zeros for the expanded bins
    counts = [zeroRow; counts; zeroRow];
end

% Decide whether the histogram will be expanded on the left and right
%   Note: Either there is data out of range or the user decides to expand
toExpandOnTheLeft = any(counts(1, :) > 0) || edgesUser(1) == -Inf;
toExpandOnTheRight = any(counts(end, :) > 0) || edgesUser(end) == Inf;

% Check for out of range data and adjust the bincounts and edges
if toExpandOnTheLeft && toExpandOnTheRight
    % The new edges are simply the expanded edges
    edgesNew = edgesExpanded;
elseif toExpandOnTheLeft && ~toExpandOnTheRight
    % Remove the last bin
    counts(end, :) = [];

    % The new edges excludes Inf
    edgesNew = edgesExpanded(1:end-1);
elseif ~toExpandOnTheLeft && toExpandOnTheRight
    % Remove the first bin
    counts(1, :) = [];

    % The new edges excludes -Inf
    edgesNew = edgesExpanded(2:end);
else
    % Remove the first and last bins
    counts(1, :) = [];
    counts(end, :) = [];

    % The new edges excludes -Inf and Inf
    edgesNew = edgesExpanded(2:end-1);
end

% The left edges of the histogram exclude the last edge
leftEdges = edgesNew(1:end-1);

% Determine the left edges for plotting 
%   If expanded to the left, replace -Inf with a finite left edge
leftEdgesPlot = leftEdges;          % initialize to be the same as left edges
if edgesNew(1) == -Inf              % data out of range on the left
    % Use the left most finite bin width to set the width for the first bin
    leftMostBinWidth = edgesExpanded(3) - edgesExpanded(2);

    % Update the left edge of the first bin
    leftEdgesPlot(1) = leftEdges(2) - leftMostBinWidth;
end

% Determine the right edge of the histogram
if edgesNew(end) == Inf             % data out of range on the right
    % Note: The right most finite bin width is used 
    %   by bar() to set the width for the last bin
    rightMostBinWidth = edgesNew(end-1) - edgesNew(end-2);
    
    % Update the right edge of the last bin
    rightMostEdgePlot = edgesNew(end-1) + rightMostBinWidth;
else                                % nothing out of range on the right
    % The right most edge is the right most finite edge
    rightMostEdgePlot = edgesNew(end);
end

% Update edges for plotting
edgesPlot = [leftEdgesPlot; rightMostEdgePlot];

% Set xLimits if not specified
if isempty(xLimits)
    xLimits = [edgesPlot(1), edgesPlot(end)];
end

% Plot histogram
if useBuiltIn
    if matlabYear >= 2017
        % Force as row vectors to be consistent with histogram()
        [edgesPlot, counts] = argfun(@force_row_vector, edgesPlot, counts);

        % Plot histogram with histogram()
        bars = histogram('BinEdges', edgesPlot, 'BinCounts', counts, ...
                      'DisplayName', 'data', otherArguments{:});
                                        % available for R2017a and beyond
    else
        % Plot histogram by using the bar() function in the 'histc' groupedStyle
        bars = bar(leftEdgesPlot, counts, 'histc', ...
                'DisplayName', 'data', otherArguments{:});
    end

    % Set y axis limits
    if ~isempty(yLimits) && ~strcmpi(yLimits, 'suppress')
        ylim(yLimits);
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
else
    % Allow the option to plot a grouped histogram
    [bars, fig] = ...
        plot_grouped_histogram('Counts', counts, 'Edges', edgesPlot, ...
                    'Style', groupedStyle, 'YLimits', yLimits, ...
                    'XUnits', xUnits, 'XLabel', xLabel, 'YLabel', yLabel, ...
                    'GroupingLabels', groupingLabels, ...
                    'LegendLocation', legendLocation, ...
                    'FigTitle', figTitle, 'FigHandle', fig, ...
                    'FigNumber', figNumber, ...
                    otherArguments{:});
end

% Count the number of bars plotted
nGroups = numel(bars);

% Set default x tick locations
if isempty(xTickLocs)
    % Initialize x tick locations with old locations
    xTicks = get(gca, 'XTick'); 

    % If too few x ticks within x limits, find new tick locations
    if ~isempty(xLimits) && ...
            sum(xTicks >= xLimits(1) & xTicks <= xLimits(2)) < minNXTicks
        % Use the minimum number of x ticks
        xTicks = extract_subvectors(edgesFinite, 'Windows', xLimits, ...
                                    'MaxNum', minNXTicks);
    end

    % Remove x ticks that are beyond finite range of edges
    %   and initialize x tick labels with these numbers
    xTicks = xTicks(xTicks >= edgesFinite(1));
    xTicks = xTicks(xTicks <= edgesFinite(end));
    xTickLabelNums = xTicks;

    % Update x ticks to include where -Inf and Inf would be placed
    if edgesNew(1) == -Inf              % data out of range on the left
        % Add -Inf as first XTick at edgesPlot(1)
        xTicks = [edgesPlot(1), xTicks];
        xTickLabelNums = [-Inf, xTickLabelNums];
    end
    if edgesNew(end) == Inf              % data out of range on the right
        % Add Inf as last XTick at edgesPlot(end)
        xTicks = [xTicks, edgesPlot(end)];
        xTickLabelNums = [xTickLabelNums, Inf];
    end
else
    xTicks = xTickLocs;
    xTickLabelNums = xTicks;
end

% Update the x tick locations
set(gca, 'XTick', xTicks);

% Create default x tick labels
if isempty(xTickLabels)
    % Create x tick labels using xTickLabelNums
    xTickLabels = create_labels_from_numbers(xTickLabelNums, 'Precision', 2);
end

% Update x tick labels
set(gca, 'XTickLabel', xTickLabels);

% Store hold status
wasHold = ishold;

% Hold on
hold on;

% Plot expanded bins
if ~isempty(xTickLabelNums) && xTickLabelNums(1) == -Inf
    if useBuiltIn
        bars(nGroups + 1) = ...
            histogram(edgesPlot(1) * ones(1, counts(1)), ...
                        edgesPlot(1:2), ...
                        'FaceAlpha', 1, 'FaceColor', specialColor, ...
                        'DisplayName', 'data too small', otherArguments{:});
    else
        bars(nGroups + 1) = ...
            bar(mean(edgesPlot(1:2)), sum(counts(1, :)), ...
                        edgesPlot(2) - edgesPlot(1), ...
                        'FaceAlpha', 1, 'FaceColor', specialColor, ...
                        'DisplayName', 'data too small', otherArguments{:});
    end
else
    bars(nGroups + 1) = gobjects(1);
end
if ~isempty(xTickLabelNums) && xTickLabelNums(end) == Inf
    if useBuiltIn
        bars(nGroups + 2) = ...
            histogram(edgesPlot(end-1) * ones(1, counts(end)), ...
                        edgesPlot(end-1:end), ...
                        'FaceAlpha', 1, 'FaceColor', specialColor, ...
                        'DisplayName', 'data too large', otherArguments{:});
    else
        bars(nGroups + 2) = ...
            bar(mean(edgesPlot(end-1:end)), sum(counts(end, :)), ...
                        edgesPlot(end) - edgesPlot(end-1), ...
                        'FaceAlpha', 1, 'FaceColor', specialColor, ...
                        'DisplayName', 'data too large', otherArguments{:});
    end
else
    bars(nGroups + 2) = gobjects(1);
end

% Hold off if it was that way before
if ~wasHold
    hold off;
end

% Set x axis limits
if ~strcmpi(xLimits, 'suppress')
    xlim(xLimits);
end

% Save the figure
if ~isempty(figName)
    % Save the figure in all file types requested
    save_all_figtypes(fig, figName, figTypes);

    % Close figure
    close(fig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if isempty(xTicks)
    xTicks = edgesPlot;
end

% Expand edges
edgesExpanded = edges;
if edgesExpanded(1) ~= -Inf
    if iscolumn(edgesExpanded)
        edgesExpanded = [-Inf; edgesExpanded];
    else
        edgesExpanded = [-Inf, edgesExpanded];
    end
end
if edgesExpanded(end) ~= Inf
    if iscolumn(edgesExpanded)
        edgesExpanded = [edgesExpanded; Inf];
    else
        edgesExpanded = [edgesExpanded, Inf];
    end
end

% Initialize flags

if edgesNew(1) == -Inf              % data out of range on the left
    if ismember(edgesPlot(1), xTicks)   % xTicks already include the left end
        % First XTickLabel should be -Inf
        xTickLabelNums(1) = -Inf;
    else
        % Need to add -Inf as first XTick at edgesPlot(1)
        xTicks = [edgesPlot(1), xTicks];
        xTickLabelNums = [-Inf, xTickLabelNums];
    end
end
if edgesNew(end) == Inf              % data out of range on the right
    if ismember(edgesPlot(end), xTicks) % xTicks already include the right end
        % Last XTickLabel should be Inf
        xTickLabelNums(end) = Inf;
    else
        % Need to add Inf as last XTick at edgesPlot(end)
        xTicks = [xTicks, edgesPlot(end)];
        xTickLabelNums = [xTickLabelNums, Inf];
    end
end

if matlabYear >= 2017
    b.CData(1, :) = specialColor;      % valid for at least R2017a
end
if matlabYear >= 2017
    b.CData(1, :) = specialColor;      % valid for at least R2017a
end

nStds = str2double(outlierMethod(1));

%       bars           - the histogram returned as a Bar object
%                   specified as a Patch (R2015a) or Bar (R2017a) object
%       h1          - the histogram for the isolated expanded left bar if exists
%                   specified as a Histogram object
%       h2          - the histogram for the isolated expanded right bar if exists
%                   specified as a Histogram object

% Use the built-in histcounts function to find the proper bin edges
[~, edges] = histcounts(XTrimmed);
counts = histcounts(X, edgesExpanded);

if iscolumn(xTickLabelNums)
    xTickLabels = cellfun(@num2str, mat2cell(xTickLabelNums, ...
                    ones(1, length(xTickLabelNums), 1)), ...
                    'UniformOutput', false);
else
    xTickLabels = cellfun(@num2str, mat2cell(xTickLabelNums, ...
                    1, ones(1, length(xTickLabelNums))), ...
                    'UniformOutput', false);
end

% Initialize x tick locations with bin centers of edgesPlot
xTicks = compute_centers_from_edges(edgesPlot); 

% Create xTickLabelNums that reflect true bin values
xTickLabelNums = xTicks;
if edgesNew(1) == -Inf
    xTickLabelNums(1) = -Inf;
end
if edgesNew(end) == Inf
    xTickLabelNums(end) = Inf;
end

barWidth = 1;                   % set this to 1 to make bars touch each other

if iscolumn(edgesFinite)
    edgesExpanded = [-Inf; edgesFinite; Inf];
else
    edgesExpanded = [-Inf, edgesFinite, Inf];
end
if iscolumn(edgesNew)
    edgesPlot = [leftEdgesPlot; rightMostEdgePlot];
else
    edgesPlot = [leftEdgesPlot, rightMostEdgePlot];
end

if numel(edgesFinite) >= 2
    binWidth = edgesFinite(2) - edgesFinite(1);
    xTicksLeft = xLimits(1) + binWidth*10;
else
    xTicksLeft = xLimits(1);
end

if ~isempty(figHandle)
    fig = figure(figHandle);
elseif ~isempty(figNumber)
    fig = figure(figNumber);
elseif ~isempty(figName)
    fig = figure;
else
    fig = gcf;
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
