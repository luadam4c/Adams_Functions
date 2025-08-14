function [timesStart, timesEnd] = create_pulse_train_times (pulseWidth, period, tMax, perVar, varargin)
%% Randomly generates start and end times for a square wave pulse train
% Usage: [timesStart, timesEnd] = create_pulse_train_times (pulseWidth, period, tMax, perVar, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [timesStart, timesEnd] = create_pulse_train_times (70, 700, 5000, 150)
%
% Outputs:
%       timesStart  - square pulse start times
%                   specified as a numeric column vector
%       timesEnd    - square pulse end times
%                   specified as a numeric column vector
%
% Arguments:
%       pulseWidth  - square pulse width
%                   specified as a numeric scalar
%       period      - pulse train average period
%                   specified as a numeric scalar
%       tMax        - maxixum time
%                   specified as a numeric scalar
%       perVar      - pulse train period variability
%                   specified as a numeric scalar
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
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
if nargin < 4
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'pulseWidth', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addRequired(iP, 'period', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addRequired(iP, 'tMax', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addRequired(iP, 'perVar', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, pulseWidth, period, tMax, perVar, varargin{:});
param1 = iP.Results.param1;

%% Preparation
% Initialize time to 0
tNow = 0;

%% Do the job
% Collect start times of square pulses
timesStart = [];
while tNow + pulseWidth < tMax
    % Add to time start
    timesStart = [timesStart; tNow];

    % Generate next period (at least pulseWidth * 2)
    periodNext = max(pulseWidth * 2, ...
                        period + (1 - 2 * rand(1)) * (perVar / 2));

    % Update next time start
    tNow = tNow + periodNext;
end

% Compute end times of square pulses
timesEnd = timesStart + pulseWidth;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%