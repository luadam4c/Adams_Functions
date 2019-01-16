function varargout = create_default_grouping (varargin)
%% Creates numeric grouping vectors and grouping labels from data, counts or original non-numeric grouping vectors
% Usage: varargout = create_default_grouping (varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       grouping        - final numeric group assignment for each data entry
%       groupingLabels  - final group labels
% Arguments:
%       varargin    - 'Grouping': group assignment for each data point
%                   must be an array of one the following types:
%                       'cell', 'string', numeric', 'logical', 
%                           'datetime', 'duration'
%                   default == the column number for a 2D array
%                   - 'GroupingLabels' - labels for the groupings 
%                                           if not to return default
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == {'Group #1', 'Group #2', ...}
%                   - 'Stats': (opt) data to distribute among bins
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%                   - 'Counts': (opt) bin counts
%                   must be a numeric array
%                   
% Requires:
%       cd/convert_to_rank.m
%       cd/create_error_for_nargin.m
%       cd/create_grouping_by_columns.m
%       cd/create_labels_from_numbers.m
%       cd/force_column_vector.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/compute_grouped_histcounts.m
%       cd/plot_grouped_histogram.m

% File History:
% 2019-01-15 Moved from plot_grouped_histogram.m
% 

%% Hard-coded parameters
% TODO: Make these optional arguments
groupingLabelPrefix = '';

%% Default values for optional arguments
groupingDefault = [];           % set later
groupingLabelsDefault = '';     % set later
statsDefault = [];              % set later
countsDefault = [];             % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Grouping', groupingDefault, ...
    @(x) validateattributes(x, {'cell', 'string', 'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'GroupingLabels', groupingLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'Stats', statsDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'Counts', countsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, varargin{:});
grouping = iP.Results.Grouping;
groupingLabels = iP.Results.GroupingLabels;
stats = iP.Results.Stats;
counts = iP.Results.Counts;

%% Do the job
if isempty(grouping)
    if ~isempty(stats)
        % Force rows as a columns
        stats = force_column_vector(stats, 'TreatCellAsArray', true, ...
                                    'IgnoreNonVectors', true);

        % Create a grouping vector from the columns
        grouping = create_grouping_by_columns(stats);
    elseif ~isempty(counts)
        % Force rows as a columns
        stats = force_column_vector(counts, 'TreatCellAsArray', true, ...
                                    'IgnoreNonVectors', true);

        % Create a grouping vector from the columns
        grouping = create_grouping_by_columns(counts);
    else
        grouping = NaN;
    end
elseif iscellstr(grouping) || isstring(grouping) 
    % Use these for grouping labels
    if isempty(groupingLabels)
        % Make unique strings the grouping labels
        groupingLabels = unique(grouping);
    end

    % Create a numeric grouping vector based on the order in the grouping labels
    grouping = convert_to_rank(grouping, 'RankedElements', groupingLabels, ...
                                'SearchMode', 'substrings');

else
    % Do nothing
end

% Set the default grouping labels
if nargout >= 2 && isempty(groupingLabels)
    % Get all unique group values
    groupValues = unique(grouping);

    % Create grouping labels from unique values
    groupingLabels = create_labels_from_numbers(groupValues, ...
                                        'Prefix', groupingLabelPrefix);
end

%% Outputs
varargout{1} = grouping;
if nargout >= 2
    varargout{2} = groupingLabels;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%