function samplingIntervals = compute_sampling_interval (timeVecs, varargin)
%% Computes sampling intervals from time vectors
% Usage: samplingIntervals = compute_sampling_interval (timeVecs, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       samplingIntervals   - sampling intervals
%                           specified as a numeric vector
% Arguments:
%       timeVecs    - time vectors
%                   must be a numeric vector or a cell array of numeric vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Used by:
%       cd/compute_all_pulse_responses.m
%       cd/compute_average_pulse_response.m
%       cd/create_average_time_vector.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/match_time_info.m
%       cd/parse_LTS.m

% File History:
% 2018-11-28 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default   = [];                   % default TODO: Description of param1

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
addRequired(iP, 'timeVecs', ...
    @(x) assert(isnumericvector(x) || iscellnumericvector(x), ...
                ['vecs must be either a numeric vector ', ...
                    'or a cell array of numeric vectors!']));

% % Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, timeVecs, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
if iscell(timeVecs)
    % Use cellfun
    samplingIntervals = cellfun(@(x) x(2) - x(1), timeVecs);
else
    % Use the first two points
    samplingIntervals = timeVecs(2) - timeVecs(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%