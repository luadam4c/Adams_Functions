function [b, h1, h2] = histogram_include_outofrange(X, edges, varargin)
%% Plot a histogram including out of range values
% Usage: [b, h1, h2] = histogram_include_outofrange(X, edges, varargin)
% Explanation:
%       Automatically combines the counts of X outside of the finite range 
%           of edges on the left or on the right to a bin on the left or 
%           on the right, respectively.
%       Note: The bar() function is actually used for the main histogram
% Outputs:
%       b           - the histogram returned as a Bar object
%                   specified as a Patch (R2015a) or Bar (R2017a) object
%       h1          - the histogram for the isolated expanded left bar if exists
%                   specified as a Histogram object
%       h2          - the histogram for the isolated expanded right bar if exists
%                   specified as a Histogram object
% Side Effects:
%       Plots a histogram
% Arguments:    
%       X           - data to distribute among bins
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%       edges       - bin edges
%                   must be a vector of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%       varargin    - 'SpecialColor': color of expanded bins
%                   must be a 3-element numeric vector:
%                   default == [0 0.8 0.8] (light blue)
%                   - 'XLimits': x-axis limits
%                   must be a two element vector of one the following types:
%                       'numeric', 'logical''datetime', 'duration'
%                   default == minimum and maximum edges of bins
%
% Used by:    
%       TODO: place any custom functions/scripts that uses this function here
%       /TODO:dir/TODO:file
%
% File History:
% 2017-12-12 Created by Adam Lu
% 

%% Default values for optional arguments
xLimitsDefault = [];
specialColorDefault = [0 0.8 0.8];          % light blue

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help histogram_include_outofrange'' for usage']);
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'histogram_include_outofrange';
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'X', ...                        % data to distribute among bins
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'nonempty'}));
addRequired(iP, 'edges', ...                    % bin edges
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SpecialColor', specialColorDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 3}));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) validateattributes(x, {'numeric', 'categorical', ...
        'datetime', 'duration'}, {'vector', 'numel', 2}));

% Read from the Input Parser
parse(iP, X, edges, varargin{:});
xLimits = iP.Results.XLimits;
specialColor = iP.Results.SpecialColor;

% Display warning message if some inputs are unmatched
if ~isempty(fieldnames(iP.Unmatched))
    fprintf('WARNING: The following name-value pairs could not be parsed: \n');
    disp(iP.Unmatched);
end

%% Perform job
% Extract finite part of edges
edgesFinite = edges(isfinite(edges));

% Expand edges
if iscolumn(edgesFinite)
    edgesExpanded = [-Inf; edgesFinite; Inf];
else
    edgesExpanded = [-Inf, edgesFinite, Inf];
end

% Compute histogram bincounts with expanded edges
%   Note: the first bin is always < the first non-Inf number in edges
%   Note: the last bin is always >= the last non-Inf number in edges
counts = histcounts(X, edgesExpanded);

% Check for out of range data and adjust the bincounts and edges
expandedOnTheLeft = false;  % whether histogram will be expanded on the left
expandedOnTheRight = false; % whether histogram will be expanded on the right
if (counts(1) > 0 || edges(1) == -Inf) && (counts(end) > 0 || edges(end) == Inf)
        % if data out of range or user specifies so on both sides
    % Histogram will be expanded on both sides
    expandedOnTheLeft = true;
    expandedOnTheRight = true;

    % The new edges are simply the expanded edges
    edgesNew = edgesExpanded;
elseif counts(1) > 0 || edges(1) == -Inf
        % if data out of range or user specifies so on the left side only
    % Histogram will be expanded on the left
    expandedOnTheLeft = true;

    % Remove the last bin
    counts(end) = [];

    % The new edges excludes Inf
    edgesNew = edgesExpanded(1:end-1);
elseif counts(end) > 0 || edges(end) == Inf
        % if data out of range or user specifies so on the right side only
    % Histogram will be expanded on the right
    expandedOnTheRight = true;

    % Remove the first bin
    counts(1) = [];

    % The new edges excludes -Inf
    edgesNew = edgesExpanded(2:end);
else    % if no data out of range
    % Remove the first and last bins
    counts(1) = [];
    counts(end) = [];

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
    % OBSERVATION: The right most finite bin width is used 
    %   by bar() to set the width for the last bin
    rightMostBinWidth = edgesNew(end-1) - edgesNew(end-2);
    
    % Update the right edge of the last bin
    rightMostEdgePlot = edgesNew(end-1) + rightMostBinWidth;
else                                % nothing out of range on the right
    % The right most edge is the right most finite edge
    rightMostEdgePlot = edgesNew(end);
end

% Update edges for plotting
if iscolumn(edgesNew)
    edgesPlot = [leftEdgesPlot; rightMostEdgePlot];
else
    edgesPlot = [leftEdgesPlot, rightMostEdgePlot];
end

% Set xLimits if not specified
if isempty(xLimits)
    xLimits = [edgesPlot(1), edgesPlot(end)];
end

% Plot histogram by using the bar() function in the 'histc' style
b = bar(leftEdgesPlot, counts, 'histc');

% Initialize XTick locations with current locations
xTicks = get(gca, 'XTick'); 

% Remove XTicks that are beyond finite range of edges
%   and initialize XTickLabels with these numbers
xTicks = xTicks(xTicks >= edgesFinite(1));
xTicks = xTicks(xTicks <= edgesFinite(end));
xTickLabelNums = xTicks;

% Update xTicks to include where -Inf and Inf would be placed
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
set(gca, 'XTick', xTicks);
% xticks(xTicks);               % valid for R2016a and beyond

% Update xTickLabels using xTickLabelNums
if iscolumn(xTickLabelNums)
    xTickLabels = cellfun(@num2str, mat2cell(xTickLabelNums, ...
                    ones(1, length(xTickLabelNums), 1)), ...
                    'UniformOutput', false);
else
    xTickLabels = cellfun(@num2str, mat2cell(xTickLabelNums, ...
                    1, ones(1, length(xTickLabelNums))), ...
                    'UniformOutput', false);
end
set(gca, 'XTickLabel', xTickLabels);
% xticklabels(xTickLabels);     % valid for R2016a and beyond

% Change bar color of expanded bins to special color
if ~ishold
    hold on;
    wasHold = false;
else
    wasHold = true;
end
if xTickLabelNums(1) == -Inf
    h1 = histogram(edgesPlot(1) * ones(1, counts(1)), ...
                    edgesPlot(1:2), ...
                    'FaceAlpha', 1, 'FaceColor', specialColor);
    % b.CData(1, :) = specialColor;      % valid for at least R2017a
else
    h1 = [];
end
if xTickLabelNums(end) == Inf
    h2 = histogram(edgesPlot(end-1) * ones(1, counts(end)), ...
                    edgesPlot(end-1:end), ...
                    'FaceAlpha', 1, 'FaceColor', specialColor);
    % b.CData(1, :) = specialColor;      % valid for at least R2017a
else
    h2 = [];
end
if ~wasHold
    hold off;
end

% Update x axis limits
xlim(xLimits);

% Newer versions of MATLAB (at least 2017a) can use histogram()
% TODO: Check version?
% histogram('BinEdges', edgesPlot, 'BinCounts', counts);

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

%}

