function [tVecAvg, minNSamples] = create_average_time_vector (tVecs)
%% Creates an average time vector from a set of time vectors
% Usage: [tVecAvg, minNSamples] = create_average_time_vector (tVecs)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       tVecAvg     - average time vector
%                   specified as a numeric vector
%       minNSamples - number of samples used for the average time vector
% Arguments:
%       tVecs       - time vectors to "average"
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%
% Requires:
%       cd/compute_sampling_interval.m
%       cd/count_samples.m
%       cd/create_time_vectors.m
%       cd/iscellnumericvector.m
%
% Used by:
%       cd/compute_average_pulse_response.m

% File History:
% 2018-12-15 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default   = [];                   % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1    % TODO: 1 might need to be changed
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'tVecs', ...
    @(x) assert(isnumeric(x) || iscellnumericvector(x), ...
                ['tVecs must be either a numeric array', ...
                    'or a cell array of numeric vectors!']));

% Read from the Input Parser
parse(iP, tVecs);

%% Preparation
% Force as a cell array
tVecs = force_column_cell(tVecs);

%% Do the job
% Count the number of samples of each time vector
nSamples = count_samples(tVecs);

% Find the minimum number of samples
minNSamples = min(nSamples);

% Get all start times
startTimes = cellfun(@(x) x(1), tVecs);

% Get a sampling interval
siMs = compute_sampling_interval(tVecs);

% Compute the average start time
averageStartTime = mean(startTimes);

% Construct a new time vector starting from the average start time
tVecAvg = create_time_vectors(minNSamples, 'SamplingIntervalMs', siMs, ...
                                'TimeStart', averageStartTime, ...
                                'TimeUnits', 'ms');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

siMs = tVecs{1}(2) - tVecs{1}(1);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%