function pulseTrain = create_pulse_train (pulse, frequency, totalDuration, varargin)
%% Create a pulse train from a pulse form, a frequency, and a total duration
% Usage: pulseTrain = create_pulse_train (pulse, frequency, totalDuration, varargin)
% Explanation:
%       Creates a pulse train based on a pulse waveform, the pulse frequency
%           in Hz and the total duration in milliseconds
% Example(s):
%       pulse = ones(10, 1) * 5;
%       pulseTrain = create_pulse_train (pulse, 100, 40)
% Outputs:
%       pulseTrain  - a train of pulses
%                   specified as a numeric column vector
% Arguments:    
%       pulse       - a pulse form
%                   must be a numeric vector
%       frequency   - pulse frequency in Hz
%                   must be a positive scalar
%       totalDuration   - total train duration in ms
%                   must be a positive scalar
%       varargin    - 'SamplingRate': sampling rate in Hz
%                   must be a positive scalar
%                   default == 10000 Hz
%
% Requires:
%
% Used by:    
%       /home/Matlab/Adams_Functions/create_pulse_train_series.m

% File History:
% 2018-08-08 Created by Adam Lu


% Hard-coded constants
MS_PER_S = 1000;                % 1000 ms per second

% Hard-coded parameters

%% Default values for optional arguments
samplingRateDefault = 10000;    % default sampling rate is 10 kHz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'pulse', ...                % a pulse form
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'frequency', ...            % pulse frequency in Hz
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addRequired(iP, 'totalDuration', ...        % total train duration in ms
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SamplingRate', samplingRateDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));

% Read from the Input Parser
parse(iP, pulse, frequency, totalDuration, varargin{:});
samplingRate = iP.Results.SamplingRate;

%% Preparation
% Make sure the pulse is a column
pulse = pulse(:);

% Compute the sampling interval in ms
siMs = (1 / samplingRate) * MS_PER_S;

% Count the pulse duration in samples
pulseDurationSamples = length(pulse);

% Count the total number of samples in the train
totalDurationSamples = totalDuration / siMs;

% Compute the pulse interval in ms
pulseIntervalMs = (1 / frequency) * MS_PER_S;

% Compute the pulse interval in samples
pulseIntervalSamples = floor(pulseIntervalMs / siMs);

%% Create the pulseTrain
% Compute the number of pulses
nPulses = floor(totalDurationSamples / pulseIntervalSamples);

% Create a single pulse
singlePulse = zeros(pulseIntervalSamples, 1);
singlePulse(1:pulseDurationSamples) = pulse;

% Create all pulses
allPulses = repmat(singlePulse, nPulses, 1);

% Find the length of all pulses
allPulsesLength = length(allPulses);

% Create the pulse train
pulseTrain = zeros(totalDurationSamples, 1);
pulseTrain(1:allPulsesLength) = allPulses;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}