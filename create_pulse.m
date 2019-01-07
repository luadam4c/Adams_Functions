function pulse = create_pulse (varargin)
%% Creates a pulse vector
% Usage: pulse = create_pulse (varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       pulse       - a pulse vector
%                   specified as a column numeric vector
% Arguments:
%       varargin    - 'SamplingInterval': sampling interval in time units
%                   must be a positive scalar
%                   default == 0.1 time units
%                   - 'PulseDelay': pulse delay in time units
%                   must be a nonnegative scalar
%                   default == 0 time units
%                   - 'PulseDuration': pulse duration in time units
%                   must be a nonnegative scalar
%                   default == 10 time units
%                   - 'PulseAmplitude': pulse amplitude in amplitude units
%                   must be a numeric scalar
%                   default == 1 amplitude unit
%                   - 'TotalDuration': total duration in time units
%                   must be a nonnegative scalar
%                   default == pulse delay + pulse duration
%
% Used by:
%       cd/create_pulse_train_series.m
%       cd/xolotl_add_current_pulse.m

% File History:
% 2018-12-13 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
samplingIntervalDefault = 0.1;  % 0.1 time units by default
pulseDelayDefault = 0;          % 0 time units delay by default
pulseDurationDefault = 10;      % 10 time units duration by default
pulseAmplitudeDefault = 1;      % 1 amplitude units by default
totalDurationDefault = [];      % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SamplingInterval', samplingIntervalDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'PulseDelay', pulseDelayDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'PulseDuration', pulseDurationDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'PulseAmplitude', pulseAmplitudeDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'TotalDuration', totalDurationDefault, ...
   @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));

% Read from the Input Parser
parse(iP, varargin{:});
samplingInterval = iP.Results.SamplingInterval;
pulseDelay = iP.Results.PulseDelay;
pulseDuration = iP.Results.PulseDuration;
pulseAmplitude = iP.Results.PulseAmplitude;
totalDuration = iP.Results.TotalDuration;

% Check if the total duration is greater than the pulse delay + duration
if totalDuration < pulseDelay + pulseDuration
    fprintf(['The total duration %g cannot be less than ', ...
            'the pulse delay %g + pulse duration %g!!\n'], ...
            totalDuration, pulseDelay, pulseDuration);
    pulse = [];
    return
end

%% Preparation
% Convert to samples
pulseDelaySamples = floor(pulseDelay/ samplingInterval);
pulseDurationSamples = floor(pulseDuration / samplingInterval);

if isempty(totalDuration)
    % Compute the default total duration in samples
    totalDurationSamples = pulseDelaySamples + pulseDurationSamples;
else
    % Convert the total duration to samples
    totalDurationSamples = floor(totalDuration / samplingInterval);
end

%% Do the job
% Initialize the pulse
pulse = zeros(totalDurationSamples, 1);

% Find the indices of the actual pulse
indPulse = pulseDelaySamples + transpose(1:pulseDurationSamples);

% Create the pulse
pulse(indPulse) = pulseAmplitude;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%