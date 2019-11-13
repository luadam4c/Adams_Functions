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
%       cd/adjust_edges.m
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/force_column_vector.m
%       cd/rmfield_custom.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/compute_grouped_histcounts.m
%       cd/compute_spike_histogram.m

% File History:
% 2018-12-28 Moved from plot_grouped_histogram.m
% 2019-09-08 Added 'FixedEdges' as an optional argument
% 

%% Hard-coded parameters
fieldsInConflictWithBinEdges = {'NumBins', 'BinWidth', 'BinMethod', 'BinLimits'};

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
otherArgumentsStruct = iP.Unmatched;

%% Preparation
if ~isempty(edges)
    rmfield_custom(otherArgumentsStruct, fieldsInConflictWithBinEdges);    
end

% Convert to an arguments list
otherArgumentsCell = struct2arglist(otherArgumentsStruct);

%% Do the job
% Compute bin counts and edges
if ~isempty(edges)
    % Use provided bin edges
    [counts, edges] = histcounts(stats, edges, otherArgumentsCell{:});
else
    % Use default bin edges
    [counts, edges] = histcounts(stats, otherArgumentsCell{:});
end

% If the edges do not contain a fixed edge, shift so that it does
if ~isempty(fixedEdges)
    % Update edges if necessary
    [edgesNew, isUpdated] = adjust_edges(edges, 'FixedEdges', fixedEdges);

    % Compute bins again if edges are updated
    if isUpdated
        % Remove any 
        rmfield_custom(otherArgumentsStruct, fieldsInConflictWithBinEdges);    

        % Update arguments list
        otherArgumentsCell = struct2arglist(otherArgumentsStruct);

        % Compute bins again
        [counts, edges] = histcounts(stats, edgesNew, otherArgumentsCell{:});
    end
end

% Force output as column vectors
[counts, edges] = argfun(@force_column_vector, counts, edges);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
