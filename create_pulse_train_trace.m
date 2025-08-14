function [timeVec, boolVec] = create_pulse_train_trace (timesStart, timesEnd, varargin)
%% Creates a square wave trace based on times start and times end for each square pulse, with amplitude of 1
% Usage: [timeVec, boolVec] = create_pulse_train_trace (timesStart, timesEnd, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [timeVec, boolVec] = create_pulse_train_trace ([1, 4, 6], [3, 5, 9]); plot(timeVec, boolVec); xlim([0, 10]), ylim([-1, 2]);
%       [timeVec, boolVec] = create_pulse_train_trace ([1, 4], [3, 5, 9])
%       [timeVec, boolVec] = create_pulse_train_trace ([1, 4, 6], [3, 5])
%
%       [timesStart, timesEnd] = create_pulse_train_times (70, 700, 5000, 150)
%       [timeVec, boolVec] = create_pulse_train_trace (timesStart, timesEnd); plot(timeVec, boolVec); 
%
% Outputs:
%       timeVec     - time vector
%                   specified as a numeric column vector
%       boolVec     - boolean vector (0 or 1)
%                   specified as a numeric column vector
%
% Arguments:
%       timesStart  - square pulse start times
%                   specified as a numeric vector
%       timesStart  - square pulse end times
%                   specified as a numeric vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_vector.m
%
% Used by:
%       \Shared\Code\vIRt\virt_moore.m

% File History:
% 2025-08-14 Pulled from virt_moore.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'timesStart', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'timesEnd', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, timesStart, timesEnd, varargin{:});
param1 = iP.Results.param1;

%% Preparation
% Make sure inputs are column vectors
timesStart = force_column_vector(timesStart);
timesEnd = force_column_vector(timesEnd);

% Compute the number of pulses
nStarts = length(timesStart);
nEnds = length(timesEnd);

% Check if they are the same
if nStarts ~= nEnds
    fprintf('Error: There are %d starts but %d ends!!\n', nStarts, nEnds);
    timeVec = [];
    boolVec = [];
    return
end

%% Do the job
% Place column vectors side by side
timesTemp = [timesStart, timesStart, timesEnd, timesEnd]';
boolTemp = [zeros(nStarts, 1), ones(nStarts, 1), ...
                ones(nStarts, 1), zeros(nStarts, 1)]';

% Linearize into vectors
timeVec = timesTemp(:);
boolVec = boolTemp(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%