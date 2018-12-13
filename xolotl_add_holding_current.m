function xolotlObject = xolotl_add_holding_current (xolotlObject, varargin)
%% Adds a holding current to the first compartment of a xolotl object
% Usage: xolotlObject = xolotl_add_holding_current (xolotlObject, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       xolotlObject    - a created neuron with simulation parameters
%                       specified as a xolotl object
% Arguments:
%       xolotlObject    - a created neuron with simulation parameters
%                       must be a xolotl object
%       varargin    - 'Amplitude': amplitude in nA 
%                   must be a numeric scalar
%                   default == 0 nA
%
% Requires:
%       cd/parse_xolotl_object.m
%
% Used by:
%       cd/m3ha_xolotl_test.m

% File History:
% 2018-12-12 Modified from xolotl_add_current_pulse.m
% TODO: Make more general by adding a 'Compartments' parameter,
%       with only the first compartment by default
% 

%% Hard-coded parameters

%% Default values for optional arguments
amplitudeDefault = 0;       % default amplitude in nA

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
addParameter(iP, 'Amplitude', amplitudeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Read from the Input Parser
parse(iP, xolotlObject, varargin{:});
amplitude = iP.Results.Amplitude;

%% Preparation
% Parse the xolotl object
parsedParams = parse_xolotl_object(xolotlObject);

% Extract parameters
nSamples = parsedParams.nSamples;
nCompartments = parsedParams.nCompartments;
previousCurrentInjections = parsedParams.externalCurrents;

%% Create holding current(s)
% Create holding current vector for the first compartment
holdingCurrent = amplitude * ones(nSamples, 1);

% Match zeros for the other compartments
holdingCurrent = [holdingCurrent, zeros(nSamples, nCompartments - 1)];

%% Add the holding current vectors to the previously set current injections
xolotlObject.I_ext = previousCurrentInjections + holdingCurrent;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%