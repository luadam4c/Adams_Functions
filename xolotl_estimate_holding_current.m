function [holdingCurrent, testObject] = ...
                xolotl_estimate_holding_current (xolotlObject, varargin)
%% Estimates the holding current necessary to match a certain holding potential
% Usage: [holdingCurrent, testObject] = ...
%               xolotl_estimate_holding_current (xolotlObject, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       holdingCurrent  - the estimated holding current
%                       specified as a numeric scalar
% Arguments:
%       xolotlObject    - a created neuron with simulation parameters
%                       specified as a xolotl object
%       holdingPotential- (opt) required holding potential (mV)
%                       must be a numeric scalar
%                       default == the first leak reversal
%       varargin    - 'CompToPatch': compartment name to patch
%                   must be a string scalar or a character vector
%                   default == set in xolotl_compartment_index.m
%                   - 'TimeToStabilize': time to stabilize in ms
%                   must be a nonnegative scalar
%                   default == 1000 ms
%
% Requires:
%       cd/xolotl_set_simparams.m
%       cd/xolotl_add_voltage_clamp.m
%       cd/xolotl_simulate.m
%       cd/xolotl_compartment_index.m
%
% Used by:
%       cd/xolotlObject_xolotl_test.m

% File History:
% 2018-12-13 Created by Adam Lu
% 2018-12-14 Now returns testObject as the second output
% 

%% Hard-coded parameters

%% Default values for optional arguments
holdingPotentialDefault = [];   % set later
compToPatchDefault = '';        % set later
timeToStabilizeDefault = 1000;  % simulate for 1000 ms by default

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

% Add optional inputs to the Input Parser
addOptional(iP, 'holdingPotential', holdingPotentialDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'CompToPatch', compToPatchDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TimeToStabilize', timeToStabilizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));

% Read from the Input Parser
parse(iP, xolotlObject, varargin{:});
holdingPotential = iP.Results.holdingPotential;
compToPatch = iP.Results.CompToPatch;
timeToStabilize = iP.Results.TimeToStabilize;

%% Preparation
% Find the index for the compartment to patch
[idxCompToPatch, compToPatch] = ...
    xolotl_compartment_index(xolotlObject, compToPatch);

% Get the leak reversal potential for the compartment to patch
leakReversal = xolotlObject.(compToPatch).get('Leak.E');

% Decide on the initial voltage
if ~isempty(leakReversal)
    % Use the leak channel reversal potential
    initialVoltage = leakReversal;
else
    % Otherwise, initialize at the holding potential
    initialVoltage = holdingPotential;
end

% Initialize a test object
testObject = xolotlObject;
% testObject = copy(xolotlObject);

%% Do the job
% Set simulation parameters for voltage clamp
testObject = ...
    xolotl_set_simparams(testObject, 'InitialVoltage', initialVoltage, ...
                                        'TimeEnd', timeToStabilize);

% Add a voltage clamp
testObject = ...
    xolotl_add_voltage_clamp(testObject, 'Compartment', compToPatch, ...
                                        'Amplitude', holdingPotential);

% Simulate the voltage clamp experiment
holdingCurrentsTest = xolotl_simulate(testObject, 'OutputType', 0);


% Find the holding current necessary to match the holding potential
holdingCurrent = holdingCurrentsTest(end, idxCompToPatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

leakReversals = xolotlObject.get('*Leak.E*']);
if ~isempty(leakReversals)
% Use the first leak channel reversal potential if any
initialVoltage = leakReversals(1);

leakReversal = xolotlObject.(compToPatch).get('conductance.E');

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%