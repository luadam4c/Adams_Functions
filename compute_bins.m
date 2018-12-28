function [counts, edges] = compute_bins (stats, varargin)
%% Computes bin counts and edges from a vector
% Usage: [counts, edges] = compute_bins (stats, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
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
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/force_column_numeric.m
%
% Used by:
%       cd/plot_grouped_histogram.m

% File History:
% 2018-12-28 Moved from plot_grouped_histogram.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
edgesDefault = [];

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

% Read from the Input Parser
parse(iP, stats, varargin{:});
edges = iP.Results.Edges;

% Keep unmatched arguments for the TODO function
otherArguments = iP.Unmatched;

%% Do the job
% Compute bin counts and edges
if ~isempty(edges)
    % Use provided bin edges
    [counts, edges] = histcounts(stats, edges);
else
    % Use default bin edges
    [counts, edges] = histcounts(stats);
end

% Force output as column vectors
[counts, edges] = argfun(@force_column_numeric, counts, edges);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%