function xolotlObject = xolotl_add_current_pulse (xolotlObject, varargin)
%% Adds a current pulse to the first compartment of a xolotl object
% Usage: xolotlObject = xolotl_add_current_pulse (xolotlObject, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       xolotlObject    - the created neuron with current pulse mechanism added
%                       specified as a xolotl object
% Arguments:
%       xolotlObject    - the created neuron with simulation parameters set
%                       must be a xolotl object
%       varargin    - 'Delay': delay in ms
%                   must be a TODO
%                   default == TODO
%                   - 'Duration': duration in ms
%                   must be a TODO
%                   default == TODO
%                   - 'Amplitude': amplitude in nA 
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_pulse.m
%
% Used by:
%       cd/m3ha_xolotl_test.m

% File History:
% 2018-12-12 Created by Adam Lu
% TODO: Make more general by adding a 'Compartments' parameter,
%       with only the first compartment by default
% 

%% Hard-coded parameters

%% Default values for optional arguments
delayDefault = 1100;        % default delay in ms
durationDefault = 10;       % default duration in ms
amplitudeDefault = -0.050;  % default amplitude in nA

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
addRequired(iP, 'xolotlObject');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Delay', delayDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'Duration', durationDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'Amplitude', amplitudeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Read from the Input Parser
parse(iP, xolotlObject, varargin{:});
delay = iP.Results.Delay;
duration = iP.Results.Duration;
amplitude = iP.Results.Amplitude;

% Check relationships between arguments
% TODO

%% Preparation
% Extract the sampling interval in ms
siMs = xolotlObject.dt;

% Extract the total duration of the simulation in ms
totalDuration = xolotlObject.t_end;

% Extract the number of neurons from the size of the default I_ext
nCompartments = size(xolotlObject.I_ext, 2);

% Get the number of samples
nSamples = floor(totalDuration / siMs);

% Create a pulse for the first compartment
pulse = create_pulse('SamplingInterval', siMs, 'PulseDelay', delay, ...
                    'PulseDuration', duration, 'PulseAmplitude', amplitude, ...
                    'TotalDuration', totalDuration);

% Match zeros for the other compartments
pulse = [pulse, zeros(nSamples, nCompartments - 1)];

%% Do the job
xolotlObject.I_ext = pulse;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%{
OLD CODE:

% Compute the sampling rate in Hz
samplingRateHz = compute_sampling_rate(siMs);

function samplingRateHz = compute_sampling_rate(siMs)
% TODO: Pull out to its own function

MS_PER_S = 1000;

siSeconds = siMs / MS_PER_S;

samplingRateHz = 1 / siSeconds;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%