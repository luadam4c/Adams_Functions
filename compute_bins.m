function [counts, edges] = compute_bins (stats, varargin)
%% Computes bin counts and edges from a vector
% Usage: [counts, edges] = compute_bins (stats, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [counts, edges] = compute_bins(rand(100, 1))
%       [counts, edges] = compute_bins(rand(100, 1), 'FixedEdges', 0.5)
%
% Outputs:
%       counts      - bin counts
%                   specified as a numeric array
%       edges       - bin edges
%                   specified as a numeric, logical, datetime or duration vector
% Arguments:
%       stats       - a statistics vector
%                   must be a numeric, logical, datetime or duration vector
%       varargin    - 'Edges': bin edges
%                   must be a numeric, logical, datetime or duration vector
%                   default == []
%                   - 'FixedEdges': numbers that must exist in bin edges
%                   must be a numeric, logical, datetime or duration vector
%                   default == []
%                   - Any other parameter-value pair for histcounts()
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/force_column_vector.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/compute_grouped_histcounts.m
%       cd/compute_psth.m
%       cd/compute_spike_histogram.m
%       cd/plot_psth.m

% File History:
% 2018-12-28 Moved from plot_grouped_histogram.m
% 2019-09-08 Added 'FixedEdges' as an optional argument
% 

%% Hard-coded parameters

%% Default values for optional arguments
edgesDefault = [];
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
    @(x) validateattributes(x, {'numeric', 'logical', ...
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
edges = iP.Results.Edges;
fixedEdges = iP.Results.FixedEdges;

% Keep unmatched arguments for the histcounts() function
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
% Compute bin counts and edges
if ~isempty(edges)
    % Use provided bin edges
    [counts, edges] = histcounts(stats, edges, otherArguments{:});
else
    % Use default bin edges
    [counts, edges] = histcounts(stats, otherArguments{:});
end

% If the edges do not contain a fixed edge, shift so that it does
if ~isempty(fixedEdges)
    if ~all(ismember(fixedEdges, edges))
        % Sort the fixed edges
        fixedEdges = sort(fixedEdges);

        % Get the center fixed edge
        centerEdge = extract_elements(fixedEdges, 'center');

        % Count the number of edges greater than the center fixed edge
        nEdgesRight = length(find(edges > centerEdge));

        % Count the number of edges less than the center fixed edge
        nEdgesLeft = length(find(edges < centerEdge));

        % Extract the average bin width
        binWidth = nanmean(diff(edges));

        % Compute a new bin width if necessary
        diffs = diff(fixedEdges);
        if ~isempty(diffs) && ~all(mod(diffs, binWidth) == 0);
            % TODO
        end

        % Compute the new bin limits
        minEdge = centerEdge - nEdgesLeft * binWidth;
        maxEdge = centerEdge + nEdgesRight * binWidth;

        % Compute the new bin edges
        edgesNew = minEdge:binWidth:maxEdge;

        % Compute bins again
        [counts, edges] = histcounts(stats, edgesNew, otherArguments{:});
    end
end

% Force output as column vectors
[counts, edges] = argfun(@force_column_vector, counts, edges);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
