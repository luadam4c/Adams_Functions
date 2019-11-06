function spikeFreqs = compute_spike_frequency (indSpikes, siMs, varargin)
%% Computes the spike frequency for sets of spike indices given a sampling interval
% Usage: spikeFreqs = compute_spike_frequency (indSpikes, siMs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       spikeFreqs  - spike frequency in Hz
%                   specified as a TODO
%
% Arguments:
%       indSpikes   - spike indices
%                   must be a numeric array or a cell array of numeric vectors
%       siMs        - sampling interval(s) in milliseconds
%                   must be a positive vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/parse_current_family.m

% File History:
% 2019-11-06 Moved from parse_current_family.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'indSpikes', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['indSpikes must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'siMs', siMsDefault, ...
    @(x) assert(ispositivevector(x), ...
                'siMs must be a positive vector!'));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, indSpikes, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Compute the spike frequency for each set of spike indices
if iscell(indSpikes)
    spikeFreqs = cellfun(@(x) compute_spike_frequency_helper(x, siMs), ...
                        indSpikes);
else
    spikeFreqs = compute_spike_frequency_helper(indSpikes, siMs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spikeFreqHz = compute_spike_frequency_helper(indSpikes, siMs)
%% Computes the spike frequency for one set of spike indices

%% Hard-coded parameters
MS_PER_S = 1000;

% Count the number of spikes
nSpikes = numel(indSpikes);

% If less than two spikes, spike frequency is zero
if nSpikes < 2
    spikeFreqHz = 0;
    return
end

% Otherwise, extract the first and last index
idxFirst = indSpikes(1);
idxLast = indSpikes(end);

% Compute the number of samples between the first and last spike
nSamplesBetweenFirstAndLast = idxLast - idxFirst;

% Compute the time in seconds between the first and last spike
timeDiffSeconds = nSamplesBetweenFirstAndLast * siMs / MS_PER_S;

% Compute the average spike frequency in Hz
spikeFreqHz = (nSpikes - 1) / timeDiffSeconds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
