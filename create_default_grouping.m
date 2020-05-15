function varargout = create_default_grouping (varargin)
%% Creates numeric grouping vectors and grouping labels from data, counts or original non-numeric grouping vectors
% Usage: [grouping, uniqueGroupValues, groupingLabels, stats] = ...
%                   create_default_grouping (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [v, u, l, s] = create_default_grouping('Stats', magic(3))
%       [v, u, l, s] = create_default_grouping('Stats', {1:5, 2:3, 6:10})
%       [v, u, l, s] = create_default_grouping('Stats', {1:5, 1:2; 1:3, 1:4})
%       [v, u, l, s] = create_default_grouping('Stats', {{1:5}, {1:3, 1:4}})
%       [v, u, l, s] = create_default_grouping('Counts', magic(3))
%       [v, u, l, s] = create_default_grouping('Grouping', {'cat', 'dog', 'rabbit'})
%
% Outputs:
%       grouping            - final numeric group assignment for each data entry
%       uniqueGroupValues   - unique grouping values
%       groupingLabels      - final group labels
%       stats               - reorganized stats array 
%                               that is the same dimesions as grouping
%
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
%                   - 'Stats': data to distribute among bins
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%                       or a cell array of such
%                   default == []
%                   - 'Counts': bin counts
%                   must be a numeric array
%                   default == []
%                   - 'TreatCellAsArray': whether to treat a cell array
%                                           as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatCellNumAsArray': whether to treat a cell array
%                                       of numeric arrays as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatCellStrAsArray': whether to treat a cell array
%                                       of character arrays as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'ToLinearize': whether to linearize a non-vector array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   
% Requires:
%       cd/convert_to_rank.m
%       cd/count_vectors.m
%       cd/create_error_for_nargin.m
%       cd/create_grouping_by_vectors.m
%       cd/create_labels_from_numbers.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
%       cd/iscellnumeric.m
%       cd/iscellnumericvector.m
%       cd/isnum.m
%       cd/struct2arglist.m
%       cd/union_over_cells.m
%       cd/unique_custom.m
%
% Used by:
%       cd/compute_grouped_histcounts.m
%       cd/compute_psth.m
%       cd/plot_grouped_histogram.m
%       cd/plot_grouped_jitter.m
%       cd/plot_grouped_scatter.m
%       cd/plot_swd_histogram.m

% File History:
% 2019-01-15 Moved from plot_grouped_histogram.m
% 2019-10-03 Now only treats cellstrs (but not cell arrays in general) as arrays
% 2019-10-04 Now uses vector numbers as grouping labels 
%               if grouping is a set of vectors
% 2019-10-04 Made uniqueGroupValues the second argument
% 2020-04-18 Now accepts cell arrays of cell arrays or numeric arrays
% 2020-04-18 Added 'ToLinearize' as an optional argument
% 2020-05-14 Fixed infinite loop

%% Hard-coded parameters
% TODO: Make these optional arguments
groupingLabelPrefix = 'Group';
ignoreEmpty = true;
useVectorCounts = false;

%% Default values for optional arguments
groupingDefault = [];           % set later
groupingLabelsDefault = '';     % set later
statsDefault = [];              % set later
countsDefault = [];             % set later
treatCellAsArrayDefault = false;% treat cell arrays as many arrays by default
treatCellNumAsArrayDefault = false; 
treatCellStrAsArrayDefault = true;  % treat cell arrays of character arrays
                                    %   as an array by default
toLinearizeDefault = false;     % whether to linearize a nonvector array

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
    @(x) isnum(x) || iscell(x));
addParameter(iP, 'Counts', countsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'TreatCellAsArray', treatCellAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellNumAsArray', treatCellNumAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellStrAsArray', treatCellStrAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ToLinearize', toLinearizeDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
grouping = iP.Results.Grouping;
groupingLabels = iP.Results.GroupingLabels;
stats = iP.Results.Stats;
counts = iP.Results.Counts;
treatCellAsArray = iP.Results.TreatCellAsArray;
treatCellNumAsArray = iP.Results.TreatCellNumAsArray;
treatCellStrAsArray = iP.Results.TreatCellStrAsArray;
toLinearize = iP.Results.ToLinearize;

%% Preparation
% If data and grouping is a cell array of cell arrays, reformat
if iscell(stats) && iscell(stats{1})
    if isempty(grouping)
        [groupingReorg, ~, ~, statsReorg] = ...
            cellfun(@(x) create_default_grouping('Stats', x, ...
                                'Counts', counts, ...
                                'TreatCellAsArray', treatCellAsArray, ...
                                'TreatCellNumAsArray', treatCellNumAsArray, ...
                                'TreatCellStrAsArray', treatCellStrAsArray, ...
                                'ToLinearize', true), ...
                    stats, 'UniformOutput', false);

        [varargout{1:nargout}] = ...
            create_default_grouping('Grouping', groupingReorg, ...
                                'Stats', statsReorg, ...
                                'GroupingLabels', groupingLabels, ...
                                'TreatCellAsArray', treatCellAsArray, ...
                                'TreatCellNumAsArray', treatCellNumAsArray, ...
                                'TreatCellStrAsArray', treatCellStrAsArray, ...
                                'ToLinearize', toLinearize);
        return
    else
        error('Not implemented yet!');
    end
end

% Force as a matrix if a cell array of numeric vectors
if iscellnumericvector(stats) || iscellnumericvector(grouping)
    [stats, grouping] = ...
        argfun(@(x) force_matrix(x, 'TreatCellAsArray', treatCellAsArray, ...
                                'TreatCellNumAsArray', treatCellNumAsArray, ...
                                'TreatCellStrAsArray', treatCellStrAsArray), ...
                                stats, grouping);
end

%% Do the job
if isempty(grouping)
    if ~isempty(stats)
        % Create a grouping vector from the vector numbers
        grouping = create_grouping_by_vectors(stats, ...
                                'TreatCellAsArray', treatCellAsArray, ...
                                'TreatCellNumAsArray', treatCellNumAsArray, ...
                                'TreatCellStrAsArray', treatCellStrAsArray);
    elseif ~isempty(counts)
        % Create a grouping vector from the column numbers
        grouping = create_grouping_by_vectors(counts);
    else
        grouping = NaN;
    end
elseif iscellstr(grouping) || isstring(grouping) 
    % Use these for grouping labels
    if isempty(groupingLabels)
        % Make unique strings the grouping labels
        groupingLabels = unique_custom(grouping);
        
        if ignoreEmpty
            groupingLabels(isemptycell(groupingLabels)) = [];
        end
    end

    % Create a numeric grouping vector based on the order in the grouping labels
    grouping = convert_to_rank(grouping, 'RankedElements', groupingLabels, ...
                                'SearchMode', 'substrings');

else
    % Do nothing
end

% Linearize as column vector if requested
%   Note: Must do this after default grouping vector creation
if toLinearize
    [stats, grouping] = ...
        argfun(@(x) force_column_vector(x, 'TreatCellAsArray', treatCellAsArray, ...
                                'TreatCellNumAsArray', treatCellNumAsArray, ...
                                'TreatCellStrAsArray', treatCellStrAsArray, ...
                                'ToLinearize', true), ...
                                stats, grouping);
end

% Determine unique group values
if nargout >= 2
    if useVectorCounts
        % Count the number of vectors
        nVectors = count_vectors(grouping);

        % The vectors numbers should be the grouping values assigned
        uniqueGroupValues = 1:nVectors;
    else
        % Get all group values
        allGroupValues = union_over_cells(grouping);

        % Get all unique grouping values
        uniqueGroupValues = unique_custom(allGroupValues, 'IgnoreNaN', true);
    end
end

% Set the default grouping labels
if nargout >= 3 && isempty(groupingLabels)
    % Create grouping labels from unique values
    groupingLabels = create_labels_from_numbers(uniqueGroupValues, ...
                                        'Prefix', groupingLabelPrefix);
end

%% Outputs
varargout{1} = grouping;
if nargout >= 2
    varargout{2} = uniqueGroupValues;
end
if nargout >= 3
    varargout{3} = groupingLabels;
end
if nargout >= 4
    varargout{4} = stats;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Force rows as a columns
stats = force_column_vector(stats, 'IgnoreNonVectors', true, ...
                        'TreatCellAsArray', treatCellAsArray, ...
                        'TreatCellNumAsArray', treatCellNumAsArray, ...
                        'TreatCellStrAsArray', treatCellStrAsArray);
% Force rows as a columns
counts = force_column_vector(counts, 'IgnoreNonVectors', true, ...
                        'TreatCellAsArray', treatCellAsArray, ...
                        'TreatCellNumAsArray', treatCellNumAsArray, ...
                        'TreatCellStrAsArray', treatCellStrAsArray);

if iscellnumeric(grouping) || isnumeric(grouping) && ~isvector(grouping)

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%