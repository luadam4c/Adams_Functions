function [avgTrace, paramsUsed] = compute_average_trace (traces, varargin)
%% Computes the average of traces that are not necessarily the same length
% Usage: [avgTrace, paramsUsed] = compute_average_trace (traces, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       avgTrace    - the average trace
%                   specified as a numeric column vector
% Arguments:    
%       traces      - traces to average
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%       varargin    - 'NSamples': number of samples in the average trace
%                   must be a nonnegative integer scalar
%                   default == minimum of the lengths of all traces
%                   - 'AlignMethod': method for alignment
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'LeftAdjust'  - Align to the left
%                       'RightAdjust' - Align to the right
%                   default == 'LeftAdjust'
%                   
% Requires:
%       cd/iscellnumeric.m
%
% Used by:
%       cd/compute_and_plot_evoked_LFP.m
%       cd/find_passive_params.m
%       cd/force_column_numeric.m
%       cd/m3ha_import_raw_traces.m

% File History:
% 2018-10-11 Created by Adam Lu
% 

%% Hard-coded parameters
validAlignMethods = {'leftadjust', 'rightadjust'};

%% Default values for optional arguments
nSamplesDefault = [];               % set later
alignMethodDefault  = 'leftadjust'; % align to the left by default

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
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['traces must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'NSamples', nSamplesDefault, ...
    @(x) validateattributes(x, {'numeric'}, ...
                                {'nonnegative', 'integer', 'scalar'}));
addParameter(iP, 'AlignMethod', alignMethodDefault, ...
    @(x) any(validatestring(x, validAlignMethods)));

% Read from the Input Parser
parse(iP, traces, varargin{:});
nSamples = iP.Results.NSamples;
alignMethod = validatestring(iP.Results.AlignMethod, validAlignMethods);

%% Preparation
% Force any row vector to be a column vector
%   but do not transform arrays
traces = force_column_numeric(traces, 'IgnoreNonVectors', true);

% Compute the number of samples for each trace
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

%% Do the job
% Return if there are no samples
if nSamples == 0
    avgTrace = [];
    return
end

% If the number of samples for each trace are not all equal,
%   align and truncate traces to the desired number of samples
if iscell(traces) && length(unique(nSamplesEachTrace)) ~= 1
    switch alignMethod
    case 'leftadjust'
        % Always start from 1
        if iscell(traces)
            tracesTruncated = cellfun(@(x) x(1:nSamples), traces, ...
                                        'UniformOutput', false);
        else
            tracesTruncated = traces(1:nSamples, :);
        end
    case 'rightadjust'
        % Always end at end
        if iscell(traces)
            tracesTruncated = cellfun(@(x) x((end-nSamples+1):end), traces, ...
                                        'UniformOutput', false);
        else
            tracesTruncated = traces((end-nSamples+1):end, :);
        end
    otherwise
        error('The align method %s is not implemented yet!!', alignMethod);
    end
else
    tracesTruncated = traces;
end

% Combine into a numeric array with the columns being vectors to be averaged
if iscell(tracesTruncated)
    % Place each column vector into a column of an array
    tracesArray = horzcat(tracesTruncated{:});
else
    % Already a numeric array with the columns being vectors to be averaged
    tracesArray = tracesTruncated;
end

% The average trace is the mean of all truncated traces (all columns)
avgTrace = mean(tracesArray, 2);

%% Output info
paramsUsed.nSamples = nSamples;
paramsUsed.alignMethod = alignMethod;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if iscell(traces) && ~all(cellfun(@iscolumn, traces))
    traces = cellfun(@(x) x(:), traces, 'UniformOutput', false);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%