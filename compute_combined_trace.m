function [combTrace, paramsUsed] = ...
                compute_combined_trace (traces, combineMethod, varargin)
%% Computes a combined trace from a set of traces
% Usage: [combTrace, paramsUsed] = ...
%               compute_combined_trace (traces, combineMethod, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       combTrace   - the combined trace
%                   specified as a numeric column vector
% Arguments:    
%       traces          - traces to average
%                       Note: If a non-vector array, each column is a vector
%                       must be a numeric array or a cell array
%       combineMethod   - method for combining traces
%                       must be an unambiguous, case-insensitive match to one of: 
%                           'average'   - take the average
%                           'maximum'   - take the maximum
%                           'minimum'   - take the minimum
%                           'all'       - take the logical AND
%                           'any'       - take the logical OR
%                           'first'     - take the first trace
%                           'last'      - take the last trace
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
%       cd/force_matrix.m
%       cd/error_unrecognized.m
%       cd/get_var_name.m
%       cd/iscellnumericvector.m
%
% Used by:
%       cd/compute_average_trace.m
%       cd/compute_maximum_trace.m
%       cd/compute_minimum_trace.m
%       cd/find_in_strings.m
%       cd/compute_average_data.m

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
% 

%% Hard-coded parameters
validAlignMethods = {'leftAdjust', 'rightAdjust', ...
                    'leftAdjustPad', 'rightAdjustPad'};
validCombineMethods = {'average', 'maximum', 'minimum', 'all', 'any', ...
                        'first', 'last'};

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
    combTrace = compute_single_combined_trace(tracesMatrix, combineMethod);
else
    % Combine all traces from each group separately
    combTrace = ...
        compute_combined_trace_each_group(tracesMatrix, grouping, combineMethod);
end

% Count the number of samples
nSamples = count_samples(combTrace);

%% Output info
paramsUsed.nSamplesEachTrace = nSamplesEachTrace;
paramsUsed.alignMethod = alignMethod;
paramsUsed.combineMethod = combineMethod;
paramsUsed.grouping = grouping;
paramsUsed.nSamples = nSamples;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combTraceEachGroup = ...
                compute_combined_trace_each_group(traces, grouping, combineMethod)
%% Computes a combined trace for each group separately

% Find unique grouping values
groups = unique(grouping, 'sorted');

% Count the number of groups
nGroups = length(groups);

count_vectors(traces);

% The length of the grouping vector must match the number of traces
if numel(grouping) ~= count_vectors(traces)
    fprintf(['The length of the grouping vector ', ...
                'does not match the number of traces!!\n']);
    combTraceEachGroup = [];
    return
end

% Combine traces from each group separately
combTraceEachGroup = cell(nGroups, 1);
for iGroup = 1:nGroups
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

    % Combine the traces from this group
    combTraceThisGroup = ...
        compute_single_combined_trace(tracesThisGroup, combineMethod);

    % Save in arrays
    combTraceEachGroup{iGroup} = combTraceThisGroup;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combTrace = compute_single_combined_trace(traces, combineMethod)
%% Computes a combined trace based on the combine method

% Combine traces
switch combineMethod
    case 'average'
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
    otherwise
        error_unrecognized(get_var_name(combineMethod), ...
                            combineMethod, mfilename);
end

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
