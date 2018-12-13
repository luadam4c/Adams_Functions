function xolotlObject = xolotl_add_current_pulse (xolotlObject, varargin)
%% Adds a current pulse to the first compartment of a xolotl object
% Usage: xolotlObject = xolotl_add_current_pulse (xolotlObject, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       xolotlObject    - a created neuron with simulation parameters set
%                       specified as a xolotl object
% Arguments:
%       xolotlObject    - a created neuron with simulation parameters set
%                       must be a xolotl object
%       varargin    - 'Delay': delay in ms
%                   must be a nonnegative scalar
%                   default == 1100 ms
%                   - 'Duration': duration in ms
%                   must be a nonnegative scalar
%                   default == 10 ms
%                   - 'Amplitude': amplitude in nA 
%                   must be a numeric scalar
%                   default == -0.050 nA
%
% Requires:
%       cd/create_pulse.m
%       cd/parse_xolotl_object.m
%
% Used by:
%       cd/m3ha_xolotl_test.m

% File History:
% 2018-12-12 Created by Adam Lu
% 2018-12-12 Now builds upon previous I_ext
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
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'Duration', durationDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'Amplitude', amplitudeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Read from the Input Parser
parse(iP, xolotlObject, varargin{:});
delay = iP.Results.Delay;
duration = iP.Results.Duration;
amplitude = iP.Results.Amplitude;

%% Preparation
% Parse the xolotl object
parsedParams = parse_xolotl_object(xolotlObject);

% Extract parameters
nSamples = parsedParams.nSamples;
nCompartments = parsedParams.nCompartments;
siMs = parsedParams.siMs;
totalDuration = parsedParams.totalDuration;
previousCurrentInjections = parsedParams.externalCurrents;

%% Create pulse vector(s)
% Create a pulse for the first compartment
pulse = create_pulse('SamplingInterval', siMs, 'PulseDelay', delay, ...
                    'PulseDuration', duration, 'PulseAmplitude', amplitude, ...
                    'TotalDuration', totalDuration);

% Match zeros for the other compartments
pulse = [pulse, zeros(nSamples, nCompartments - 1)];

%% Add the pulse vectors to the previously set current injections
xolotlObject.I_ext = previousCurrentInjections + pulse;

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

compartments = xolotlObject.find('compartment');
nCompartments = length(compartments);

xolotlObject.I_ext = pulse;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%