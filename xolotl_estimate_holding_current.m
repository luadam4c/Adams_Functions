function [holdingCurrent, xolotlObject] = ...
                xolotl_estimate_holding_current (xolotlObject, varargin)
%% Estimates the holding current necessary to match a certain holding potential
% Usage: [holdingCurrent, xolotlObject] = ...
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
%       cd/m3ha_xolotl_test.m

% File History:
% 2018-12-13 Created by Adam Lu
% 2018-12-14 Now returns xolotlObject as the second output
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
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'xolotlObject');

% Add optional inputs to the Input Parser
addOptional(iP, 'holdingPotential', holdingPotentialDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

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
idxCompToPatch = xolotl_compartment_index(xolotlObject, compToPatch);

% Get the leak reversal potential (mV) for all compartments
leakReversalAll = xolotlObject.get('*Leak.E');

%{
% Get the holding potential for the compartment to patch
holdingPotentialToPatch = holdingPotential(idxCompToPatch);

% Get the leak conductance (uS/mm^2) for all compartments
leakConductanceAll = xolotlObject.get('*Leak.gbar');

% Get the radius (mm) for all compartments
radiusAll = xolotlObject.get('*radius');

% Get the length (mm) for all compartments
lengthAll = xolotlObject.get('*len');

%% Estimate holding current (nA) by calculation
% Compute the surface area for each compartment
surfaceAreaAll = 2 * pi * radiusAll .* lengthAll;

% Compute the total leak conductance for each compartment
totalLeakConductanceAll = leakConductanceAll .* surfaceAreaAll;

% Compute the leak currents in all compartments
%   I = g (V - E)
leakCurrents = totalLeakConductanceAll .* ...
                    (holdingPotentialToPatch - leakReversalAll);

% The holding current should compensate for all leak currents
holdingCurrent = sum(leakCurrents);
%}

%% Save original parameters
% Retrieve original simulation parameters
vClampOrig = xolotlObject.V_clamp;
externalCurrentOrig = xolotlObject.I_ext;
timeEndOrig = xolotlObject.t_end;

% Extract all compartments
compartments = xolotlObject.Children;

% Count the number of compartments
nCompartments = numel(compartments);

% Retrieve original initial voltages
initialVoltageOrig = zeros(nCompartments, 1);
for iCompartment = 1:nCompartments
    % Get the current compartment name
    compartmentName = compartments{iCompartment};

    % Set the initial voltage for the current compartment
    initialVoltageOrig(iCompartment) = xolotlObject.(compartmentName).V;
end

%% Estimate holding current with a voltage clamp simulation
% Decide on the initial voltage
if ~isempty(leakReversalAll(idxCompToPatch))
    % Use the leak channel reversal potential
    initialVoltage = leakReversalAll(idxCompToPatch);
else
    % Otherwise, initialize at the holding potential
    initialVoltage = holdingPotential;
end

% Set simulation parameters for voltage clamp
xolotl_set_simparams(xolotlObject, 'InitialVoltage', initialVoltage, ...
                                    'TimeEnd', timeToStabilize);

% Add a voltage clamp
%   Note: this will set any external current to be zero
xolotl_add_voltage_clamp(xolotlObject, 'Compartment', compToPatch, ...
                                    'Amplitude', holdingPotential);

% Simulate the voltage clamp experiment
holdingCurrentsTest = xolotl_simulate(xolotlObject, 'OutputType', 0);

% Find the holding current necessary to match the holding potential
holdingCurrent = holdingCurrentsTest(end, idxCompToPatch);

%% Restore
% Restore simulation parameters
xolotl_set_simparams(xolotlObject, 'InitialVoltage', initialVoltageOrig, ...
                                    'TimeEnd', timeEndOrig);

% Restore the original voltage clamp or external current, 
%   based on whether it was set or not
if all(all(isnan(vClampOrig)))
    xolotlObject.I_ext = externalCurrentOrig;
else
    xolotlObject.V_clamp = vClampOrig;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

leakReversals = xolotlObject.get('*Leak.E*']);
if ~isempty(leakReversals)
% Use the first leak channel reversal potential if any
initialVoltage = leakReversals(1);

leakReversal = xolotlObject.(compToPatch).get('conductance.E');

% Initialize a test object
testObject = xolotlObject;
% testObject = copy(xolotlObject);

% Parse the xolotl object
parsedParams = parse_xolotl_object(xolotlObject);

% Extract the compartments
compartments = parsedParams.compartments;
nCompartments = parsedParams.nCompartments;

% Get the leak reversal potential (mV) for the compartment to patch
leakReversal = xolotlObject.(compToPatch).get('Leak.E');

% Get the leak conductance (uS/mm^2) for the compartment to patch
leakConductance = xolotlObject.(compToPatch).get('Leak.gbar');

% Get the surface area (mm^2) for the compartment to patch
surfaceArea = xolotlObject.(compToPatch).get('A');

holdingCurrent = leakConductance * surfaceArea * ...
                    (holdingPotentialToPatch - leakReversal);

holdingCurrent = leakCurrents(idxCompToPatch);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%