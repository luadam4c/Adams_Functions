function [combTrace, paramsUsed] = ...
                compute_combined_trace (traces, combineMethod, varargin)
%% Computes a combined trace from a set of traces
% Usage: [combTrace, paramsUsed] = ...
%               compute_combined_trace (traces, combineMethod, varargin)
% Explanation:
%       TODO
% Example(s):
%       vecs = {[1;3;4], [2;2;5], [6;6;6]};
%       vecs2 = [[1;3;4], [2;2;5], [6;6;6]];
%       compute_combined_trace(vecs, 'mean', 'Group', {'a', 'a', 'b'})
%       compute_combined_trace(vecs, 'max', 'Group', {'a', 'a', 'b'})
%       compute_combined_trace(vecs, 'min', 'Group', {'a', 'a', 'b'})
%       compute_combined_trace(vecs, 'bootmean', 'Group', {'a', 'a', 'b'})
%       compute_combined_trace(vecs2, 'bootmean', 'Group', {'a', 'a', 'b'})
% Outputs:
%       combTrace       - the combined trace(s)
%                           If grouped, a cell array is returned
%                               with the result from each group in each cell
%                       specified as a numeric column vector
% Arguments:    
%       traces          - traces to average
%                       Note: If a non-vector array, each column is a vector
%                       must be a numeric array or a cell array
%       combineMethod   - method for combining traces
%                       must be an unambiguous, case-insensitive match to one of: 
%                           'average' or 'mean' - take the average
%                           'maximum'   - take the maximum
%                           'minimum'   - take the minimum
%                           'all'       - take the logical AND
%                           'any'       - take the logical OR
%                           'first'     - take the first trace
%                           'last'      - take the last trace
%                           'bootmeans' - bootstrapped averages
%                           'bootmax'   - bootstrapped maximums
%                           'bootmin'   - bootstrapped minimums
%       varargin    - 'NSamples': number of samples in the average trace
%                   must be a nonnegative integer scalar
%                   default == minimum of the lengths of all traces
%                   - 'AlignMethod': method for alignment
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'leftAdjust'  - align to the left
%                       'rightAdjust' - align to the right
%                   default == 'leftAdjust'
%                   - 'TreatRowAsMatrix': whether to treat a row vector
%                                           as many one-element vectors
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Grouping': a grouping vector used to group traces
%                   must be a vector
%                   default == []
%                   
% Requires:
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/isnum.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
%       cd/error_unrecognized.m
%       cd/get_var_name.m
%       cd/iscellnumericvector.m
%
% Used by:
%       cd/compute_combined_data.m
%       cd/compute_average_trace.m
%       cd/compute_maximum_trace.m
%       cd/compute_minimum_trace.m
%       cd/find_in_strings.m
%       cd/compute_combined_data.m

% File History:
% 2019-01-03 Moved from compute_average_trace
% 2019-01-03 Added 'CombineMethod' as an optional argument
% 2019-01-03 Now allows NaNs
% 2019-01-03 Now uses count_samples.m
% 2019-01-03 Added 'TreatRowAsMatrix' as an optional argument
% 2019-01-04 Added 'all', 'any' as valid combine methods
% 2019-01-04 Now uses isnum.m
% 2019-01-12 Added 'Grouping' as an optional parameter
% 2019-01-12 Added 'first', 'last' as valid combine methods
% 2019-01-12 Added 'bootmeans', 'bootmax', 'bootmin as valid combine methods
% TODO: Make 'Seeds' an optional argument
% 

%% Hard-coded parameters
validAlignMethods = {'leftAdjust', 'rightAdjust', ...
                    'leftAdjustPad', 'rightAdjustPad'};
validCombineMethods = {'average', 'mean', 'maximum', 'minimum', ...
                        'all', 'any', 'first', 'last', ...
                        'bootmean', 'bootmax', 'bootmin'};

%% Default values for optional arguments
nSamplesDefault = [];               % set later
alignMethodDefault = 'leftadjust';  % align to the left by default
treatRowAsMatrixDefault = false;    % treat a row vector as a vector by default
groupingDefault = [];               % no grouping by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'traces', ...                   % vectors
    @(x) assert(isnum(x) || iscell(x), ...
                'traces must be either a numeric array or a cell array!'));
addRequired(iP, 'CombineMethod', ...
    @(x) any(validatestring(x, validCombineMethods)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'NSamples', nSamplesDefault, ...
    @(x) validateattributes(x, {'numeric'}, ...
                                {'nonnegative', 'integer', 'scalar'}));
addParameter(iP, 'AlignMethod', alignMethodDefault, ...
    @(x) any(validatestring(x, validAlignMethods)));
addParameter(iP, 'TreatRowAsMatrix', treatRowAsMatrixDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Grouping', groupingDefault);

% Read from the Input Parser
parse(iP, traces, combineMethod, varargin{:});
nSamples = iP.Results.NSamples;
alignMethod = validatestring(iP.Results.AlignMethod, validAlignMethods);
treatRowAsMatrix = iP.Results.TreatRowAsMatrix;
grouping = iP.Results.Grouping;

% Validate combine method
combineMethod = validatestring(combineMethod, validCombineMethods);

%% Do the job
if iscellnumericvector(traces) || ~iscell(traces)
    % Compute combined trace for a set of vectors
    [combTrace, paramsUsed] = ...
        compute_combined_trace_helper(traces, nSamples, grouping, ...
                                alignMethod, combineMethod, treatRowAsMatrix);
else
    % Compute combined traces for many sets of vectors
    [combTrace, paramsUsed] = ...
        cellfun(@(x) compute_combined_trace_helper(x, nSamples, grouping, ...
                    alignMethod, combineMethod, treatRowAsMatrix), ...
                traces, 'UniformOutput', false);
end

%% Make the output the same data type as the input
if ~iscell(traces) && iscell(combTrace)
    combTrace = horzcat(combTrace{:});
elseif iscellnumericvector(traces) && ~iscellnumericvector(combTrace)
    combTrace = force_column_vector(combTrace);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [combTrace, paramsUsed] = ...
        compute_combined_trace_helper(traces, nSamples, grouping, ...
                                alignMethod, combineMethod, treatRowAsMatrix)

%% Preparation
% Force any row vector to be a column vector
%   but do not transform arrays
if ~treatRowAsMatrix
    traces = force_column_vector(traces, 'IgnoreNonVectors', true);
end

% Compute the number of samples for each trace
nSamplesEachTrace = count_samples(traces, 'TreatRowAsMatrix', treatRowAsMatrix);

%% Do the job
% Force traces as a matrix and align appropriately
% TODO: Restrict the number of samples if provided
tracesMatrix = force_matrix(traces, 'AlignMethod', alignMethod);
% tracesMatrix = force_matrix(traces, 'AlignMethod', alignMethod, ...
%                               'NSamples', nSamples);

% Combine traces
if isempty(grouping)
    % Combine all traces
    [combTrace, seeds] = ...
        compute_single_combined_trace(tracesMatrix, combineMethod);
else
    % Combine all traces from each group separately
    [combTrace, seeds] = ...
        compute_combined_trace_each_group(tracesMatrix, grouping, combineMethod);
end

% Count the number of samples
nSamples = count_samples(combTrace);

%% Output info
paramsUsed.nSamplesEachTrace = nSamplesEachTrace;
paramsUsed.alignMethod = alignMethod;
paramsUsed.combineMethod = combineMethod;
paramsUsed.grouping = grouping;
paramsUsed.seeds = seeds;
paramsUsed.nSamples = nSamples;

%% Make the output the same data type as the input
if ~iscell(traces) && iscell(combTrace)
    combTrace = horzcat(combTrace{:});
elseif iscellnumericvector(traces) && ~iscellnumericvector(combTrace)
    combTrace = force_column_vector(combTrace);
    combTrace = match_dimensions(combTrace, size(traces));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [combTraces, seeds] = ...
                compute_combined_trace_each_group(traces, grouping, combineMethod)
%% Computes a combined trace for each group separately

% Find unique grouping values
groups = unique(grouping, 'sorted');

% Count the number of groups
nGroups = length(groups);

% The length of the grouping vector must match the number of traces
if numel(grouping) ~= count_vectors(traces)
    fprintf(['The length of the grouping vector ', ...
                'does not match the number of traces!!\n']);
    combTraceEachGroup = [];
    return
end

% Combine traces from each group separately
combTraceEachGroup = cell(nGroups, 1);
seeds = struct('Type', '', 'Seed', NaN, 'State', NaN);
parfor iGroup = 1:nGroups
    % Get the current grouping value
    groupValueThis = groups(iGroup);

    % Get all indices with the current grouping value
    if istext(groupValueThis)
        indThisGroup = strcmp(grouping, groupValueThis);
    else
        indThisGroup = grouping == groupValueThis;
    end
    
    % Collect all traces with this grouping value
    tracesThisGroup = traces(:, indThisGroup);

    % Save the current seed of the random number generator
    s = rng;

    % Save in output
    seeds(iGroup, 1) = s;

    % Combine the traces from this group
    combTraceEachGroup{iGroup} = ...
        compute_single_combined_trace(tracesThisGroup, combineMethod, s);
end

% Concatenate into a single matrix
combTraces = horzcat(combTraceEachGroup{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [combTrace, seed] = ...
                compute_single_combined_trace(traces, combineMethod, seed)
%% Computes a combined trace based on the combine method

% Save the current seed of the random number generator
seed = rng;

% Combine traces
switch combineMethod
    case {'average', 'mean'}
        % Take the mean of all columns
        combTrace = nanmean(traces, 2);
    case 'maximum'
        % Take the maximum of all columns
        combTrace = max(traces, [], 2);
    case 'minimum'
        % Take the minimum of all columns
        combTrace = min(traces, [], 2);
    case 'all'
        % Take the logical AND of all columns
        combTrace = all(traces, 2);
    case 'any'
        % Take the logical OR of all columns
        combTrace = any(traces, 2);
    case 'first'
        % Take the first column
        combTrace = traces(:, 1);
    case 'last'
        % Take the last column
        combTrace = traces(:, end);
    case 'bootmean'
        % Take the bootstrapped averages
        combTrace = compute_bootstrapped_combos(traces, 'mean', seed);
    case 'bootmax'
        % Take the bootstrapped maximums
        combTrace = compute_bootstrapped_combos(traces, 'maximum', seed);
    case 'bootmin'
        % Take the bootstrapped minimums
        combTrace = compute_bootstrapped_combos(traces, 'minimum', seed);
    otherwise
        error_unrecognized(get_var_name(combineMethod), ...
                            combineMethod, mfilename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combTraces = compute_bootstrapped_combos (traces, method, seed)

% Count the number of traces
nTraces = size(traces, 2);

% Seed the random number generator
rng(seed);

% Generate nTraces X nTraces samples of trace indices with replacement
selections = randi(nTraces, nTraces);

% Take the bootstrapped averages
combTraceCell = ...
    arrayfun(@(x) compute_single_combined_trace(...
                        traces(:, selections(:, x)), method), ...
                transpose(1:nTraces), 'UniformOutput', false);

% Combine them to one 2-D non-cell array
combTraces = horzcat(combTraceCell{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if iscell(traces) && ~all(cellfun(@iscolumn, traces))
    traces = cellfun(@(x) x(:), traces, 'UniformOutput', false);
end

error('The align method %s is not implemented yet!!', alignMethod);

switch alignMethod
case 'leftadjust'
    % Always start from 1
    if iscell(traces)
        tracesAligned = cellfun(@(x) x(1:nSamples), traces, ...
                                    'UniformOutput', false);
    else
        tracesAligned = traces(1:nSamples, :);
    end
case 'rightadjust'
    % Always end at end
    if iscell(traces)
        tracesAligned = cellfun(@(x) x((end-nSamples+1):end), traces, ...
                                    'UniformOutput', false);
    else
        tracesAligned = traces((end-nSamples+1):end, :);
    end
otherwise
    error_unrecognized(get_var_name(alignMethod), alignMethod, mfilename);
end

validAlignMethods = {'leftadjust', 'rightadjust'};

% If the number of samples for each trace are not all equal,
%   align and truncate traces to the desired number of samples
tracesAligned = extract_subvectors(traces, 'AlignMethod', alignMethod);

% Combine into a numeric array with the columns being vectors to be averaged
if iscell(tracesAligned)
    % Place each column vector into a column of an array
    tracesMatrix = horzcat(tracesAligned{:});
else
    % Already a numeric array with the columns being vectors to be averaged
    tracesMatrix = tracesAligned;
end

% Force any row vector to be a column vector
%   but do not transform arrays
traces = force_column_vector(traces, 'IgnoreNonVectors', true);

if iscell(traces)
    % Apply length() to each element
    nSamplesEachTrace = cellfun(@length, traces);
else
    % Whether multiple vectors or not, nSamplesEachTrace is the number of rows
    nSamplesEachTrace = size(traces, 1);
end

% Set default number of samples for the averaged trace
if isempty(nSamples)
    if iscell(traces)
        % Use the minimum length of all traces
        nSamples = min(nSamplesEachTrace);
    else
        % Use the number of rows
        nSamples = nSamplesEachTrace;
    end
end

% Return if there are no samples
if nSamples == 0
    combTrace = [];
    return
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
