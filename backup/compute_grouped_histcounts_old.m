function varargout = compute_grouped_histcounts (stats, varargin)
%% Computes bin counts and edges from grouped data
% Usage: varargout = compute_grouped_histcounts (stats, grouping (opt), varargin)
% Explanation:
%       This is similar to histcounts() but returns a 2-D array
%           if a grouping vector is provided
%
% Example(s):
%       randVec = randn(100, 1);
%       stats1 = [randVec, randVec + 1, randVec - 1];
%       [counts, edges, centers] = compute_grouped_histcounts(stats1)
%       stats2 = [randVec; randVec + 1; randVec - 1];
%       grouping2 = [repmat({'Mark'}, 100, 1); repmat({'Peter'}, 100, 1); repmat({'Katie'}, 100, 1)];
%       [counts, edges, centers] = compute_grouped_histcounts(stats2, grouping2)
%
% Outputs:
%       counts      - bin counts, with each group being a different column
%                   specified as a an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%       edges       - bin edges used
%                   specified as a vector of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%       binCenters  - bin centers
%                   specified as a vector of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%
% Arguments:
%       stats       - data to distribute among bins
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%       grouping    - (opt) group assignment for each data point
%                   must be an array of one the following types:
%                       'cell', 'string', numeric', 'logical', 
%                           'datetime', 'duration'
%                   default == the column number for a 2D array
%       varargin    - 'Edges': bin edges
%                   must be a vector of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%                   default == automatic detection of 
%                   - 'FixedEdges': numbers that must exist in bin edges
%                   must be a numeric, logical, datetime or duration vector
%                   default == []
%                   - Any other parameter-value pair for histcounts()
%
% Requires:
%       cd/adjust_edges.m
%       cd/argfun.m
%       cd/compute_bins.m
%       cd/compute_centers_from_edges.m
%       cd/create_default_grouping.m
%       cd/create_error_for_nargin.m
%       cd/force_column_vector.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/compute_activation_profile.m
%       cd/compute_psth.m
%       cd/plot_grouped_histogram.m
%       cd/plot_histogram.m
%       cd/plot_psth.m

% File History:
% 2019-01-15 Moved from plot_grouped_histogram.m
% 2019-10-05 Now does not adjust edges to fixed edges if edges is provided
% 2019-10-12 Fixed bug in computing bin edges

%% Hard-coded parameters

%% Default values for optional arguments
groupingDefault = [];           % set later
edgesDefault = [];              % set later
fixedEdgesDefault = [];

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
addRequired(iP, 'stats', ...
    @(x) validateattributes(x, {'cell', 'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'grouping', groupingDefault, ...
    @(x) validateattributes(x, {'cell', 'string', 'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Edges', edgesDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'FixedEdges', fixedEdgesDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));

% Read from the Input Parser
parse(iP, stats, varargin{:});
grouping = iP.Results.grouping;
edges = iP.Results.Edges;
fixedEdges = iP.Results.FixedEdges;

% Keep unmatched arguments for the histcounts() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Decide on the grouping vector
[grouping, uniqueGroupValues] = ...
    create_default_grouping('Stats', stats, 'Grouping', grouping);

% Count the number of groups
nGroups = numel(uniqueGroupValues);

% Force non-vectors as cell arrays of numeric vectors
[stats, grouping] = ...
    argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', false), ...
            stats, grouping);

% If stats and grouping are cell arrays of numeric vectors, pool them
if iscellnumeric(stats) && iscellnumeric(grouping)
    [stats, grouping] = argfun(@(x) vertcat(x{:}), stats, grouping);
end

%% Break up stats into a cell array of vectors
statsCell = arrayfun(@(x) stats(grouping == uniqueGroupValues(x)), ...
                    transpose(1:nGroups), 'UniformOutput', false);

%% Compute default bin edges for all data if not provided
if isempty(edges)
    % Compute bin edges for each group
    [~, edgesAll] = ...
        cellfun(@(x) compute_bins(x, 'Edges', edges, ...
                            'FixedEdges', fixedEdges, otherArguments{:}), ...
                statsCell, 'UniformOutput', false);

    % Compute the minimum bin width across groups
    minBinWidth = min(extract_elements(edgesAll, 'firstdiff'));

    % Compute the minimum and maximum edges across groups
    minEdges = min(extract_elements(edgesAll, 'first'));
    maxEdges = max(extract_elements(edgesAll, 'last'));

    % Compute the range of edges
    rangeEdges = maxEdges - minEdges;

    % Compute the number of bins
    nBins = ceil(rangeEdges / minBinWidth);

    % Create bin edges that works for all data
    edges = transpose(linspace(minEdges, maxEdges, nBins + 1));
    
    %% Modify edges to include fixed edges if necessary
    edges = adjust_edges(edges, 'FixedEdges', fixedEdges);
end

%% Compute the bin counts for each group based on these edges
%   Note: edges are not updated here
counts = cellfun(@(x) compute_bins(x, 'Edges', edges, otherArguments{:}), ...
                statsCell, 'UniformOutput', false);
counts = horzcat(counts{:});

%% Compute the bin centers
if nargout >= 3
    binCenters = compute_centers_from_edges(edges);
end

%% Outputs
varargout{1} = counts;
varargout{2} = edges;
if nargout >= 3
    varargout{3} = binCenters;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Get all unique group values
uniqueGroupValues = unique_custom(grouping, 'IgnoreNaN', true);

% Count the number of groups
nGroups = numel(uniqueGroupValues);

% Find all unique grouping values
if iscellnumeric(grouping)
    % Get all unique grouping values from each vector
    uniqueGroupValues = cellfun(@(x) unique_custom(x, 'IgnoreNaN', true), ...
                                grouping, 'UniformOutput', false);

    % Replace any empty values with a new value
else
    % Get unique grouping values
    uniqueGroupValues = unique_custom(x, 'IgnoreNaN', true);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
