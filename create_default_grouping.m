function [grouping, groupingLabels] = create_default_grouping (stats, varargin)
%% Creates default grouping vectors and grouping labels from data
% Usage: [grouping, groupingLabels] = create_default_grouping (stats, grouping, groupingLabels, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       grouping        - final group assignment for each data point
%       groupingLabels  - final group labels
% Arguments:
%       stats       - data to distribute among bins
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%       grouping        - (opt) group assignment for each data point
%                       must be an array of one the following types:
%                           'cell', 'string', numeric', 'logical', 
%                               'datetime', 'duration'
%                       default == the column number for a 2D array
%       groupingLabels  - labels for the groupings if not to return default
%                       must be a string scalar or a character vector 
%                           or a cell array of strings or character vectors
%                       default == {'Group #1', 'Group #2', ...}
%
% Requires:
%       cd/convert_to_rank.m
%       cd/create_error_for_nargin.m
%       cd/create_grouping_by_columns.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/compute_grouped_histcounts.m
%       cd/plot_grouped_histogram.m

% File History:
% 2019-01-15 Moved from plot_grouped_histogram.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
groupingDefault = [];           % set later
groupingLabelsDefault = '';     % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'stats', ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'grouping', groupingDefault, ...
    @(x) validateattributes(x, {'cell', 'string', 'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addOptional(iP, 'groupingLabels', groupingLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));

% Read from the Input Parser
parse(iP, stats, varargin{:});
grouping = iP.Results.grouping;
groupingLabels = iP.Results.groupingLabels;

%% Do the job
if isempty(grouping)
    % Create a grouping if not provided
    grouping = create_grouping_by_columns(stats);
elseif iscellstr(grouping) || isstring(grouping) 
    % Use these for grouping labels
    if isempty(groupingLabels)
        % Make unique strings the grouping labels
        groupingLabels = unique(grouping);
    end

    % Create a numeric grouping vector based on the order in the grouping labels
    grouping = convert_to_rank(grouping, 'RankedElements', groupingLabels, ...
                                'SearchMode', 'substrings');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%