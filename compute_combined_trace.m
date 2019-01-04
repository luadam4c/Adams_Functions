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
%                           'average' - take the average
%                           'maximum' - take the maximum
%                           'minimum' - take the minimum
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
%                   
% Requires:
%       cd/count_samples.m
%       cd/force_matrix.m
%       cd/error_unrecognized.m
%       cd/get_var_name.m
%       cd/iscellnumericvector.m
%
% Used by:
%       cd/compute_average_trace.m
%       cd/compute_maximum_trace.m
%       cd/compute_minimum_trace.m

% File History:
% 2019-01-03 Moved from compute_average_trace
% 2019-01-03 Added 'CombineMethod' as an optional argument
% 2019-01-03 Now allows NaNs
% 2019-01-03 Now uses count_samples.m
% 2019-01-03 Added 'TreatRowAsMatrix' as an optional argument
% 

%% Hard-coded parameters
validAlignMethods = {'leftAdjust', 'rightAdjust', ...
                    'leftAdjustPad', 'rightAdjustPad'};
validCombineMethods = {'average', 'maximum', 'minimum'};

%% Default values for optional arguments
nSamplesDefault = [];               % set later
alignMethodDefault = 'leftadjust';  % align to the left by default
treatRowAsMatrixDefault = false;    % treat a row vector as a vector by default

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
    @(x) assert(isnumeric(x) || iscell(x), ...
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

% Read from the Input Parser
parse(iP, traces, combineMethod, varargin{:});
nSamples = iP.Results.NSamples;
alignMethod = validatestring(iP.Results.AlignMethod, validAlignMethods);
treatRowAsMatrix = iP.Results.TreatRowAsMatrix;

% Validate combine method
combineMethod = validatestring(combineMethod, validCombineMethods);

%% Do the job
if iscellnumericvector(traces) || ~iscell(traces)
    [combTrace, paramsUsed] = ...
        compute_one_combined_trace(traces, nSamples, ...
                                alignMethod, combineMethod, treatRowAsMatrix);
else
    [combTrace, paramsUsed] = ...
        cellfun(@(x) compute_one_combined_trace(x, nSamples, ...
                    alignMethod, combineMethod, treatRowAsMatrix), ...
                traces, 'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [combTrace, paramsUsed] = ...
        compute_one_combined_trace(traces, nSamples, alignMethod, ...
                                    combineMethod, treatRowAsMatrix)

%% Preparation
% Force any row vector to be a column vector
%   but do not transform arrays
if ~treatRowAsMatrix
    traces = force_column_numeric(traces, 'IgnoreNonVectors', true);
end

% Compute the number of samples for each trace
nSamplesEachTrace = count_samples(traces, 'TreatRowAsMatrix', treatRowAsMatrix);

%% Do the job
% Force traces as a matrix and align appropriately
tracesMatrix = force_matrix(traces, 'AlignMethod', alignMethod);

% Combine traces
switch combineMethod
    case 'average'
        % Take the mean of all truncated traces (all columns)
        combTrace = nanmean(tracesMatrix, 2);
    case 'maximum'
        % Take the maximum of all truncated traces (all columns)
        combTrace = max(tracesMatrix, [], 2);
    case 'minimum'
        % Take the minimum of all truncated traces (all columns)
        combTrace = min(tracesMatrix, [], 2);
    otherwise
        error_unrecognized(get_var_name(combineMethod), ...
                            combineMethod, mfilename);
end

% Count the number of samples
nSamples = count_samples(combTrace);

%% Output info
paramsUsed.nSamplesEachTrace = nSamplesEachTrace;
paramsUsed.alignMethod = alignMethod;
paramsUsed.combineMethod = combineMethod;
paramsUsed.nSamples = nSamples;

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
traces = force_column_numeric(traces, 'IgnoreNonVectors', true);

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
