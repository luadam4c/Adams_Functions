function [xolotlObject, clampedVoltage] = ...
                xolotl_add_voltage_clamp (xolotlObject, varargin)
%% Adds a voltage clamp to the first compartment of a xolotl object
% Usage: [xolotlObject, clampedVoltage] = ...
%               xolotl_add_voltage_clamp (xolotlObject, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       xolotlObject    - a created neuron with simulation parameters
%                       specified as a xolotl object
%       clampedVoltage  - clamped voltage traces added
%                       specified as a numeric array
% Arguments:
%       xolotlObject    - a created neuron with simulation parameters
%                       must be a xolotl object
%       varargin    - 'Compartment': the compartment name to add clamp
%                   must be a string scalar or a character vector
%                   default == set in xolotl_compartment_index.m
%                   - 'Amplitude': amplitude in mV 
%                   must be a numeric scalar
%                   default == -70 mV
%
% Requires:
%       cd/parse_xolotl_object.m
%       cd/xolotl_compartment_index.m
%
% Used by:
%       cd/m3ha_xolotl_test.m

% File History:
% 2018-12-13 Modified from xolotl_add_voltage_clamp.m
% TODO: Make more general by adding a 'Compartment' parameter,
%       with only the first compartment by default
% TODO: vcRs = 0.01             // Voltage clamp series resistance in MOhm
% 

%% Hard-coded parameters

%% Default values for optional arguments
compartmentDefault = '';    % set later
amplitudeDefault = -70;     % default amplitude in mV

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
addParameter(iP, 'Compartment', compartmentDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Amplitude', amplitudeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Read from the Input Parser
parse(iP, xolotlObject, varargin{:});
compartment = iP.Results.Compartment;
amplitude = iP.Results.Amplitude;

%% Preparation
% Parse the xolotl object
parsedParams = parse_xolotl_object(xolotlObject);

% Extract parameters
nSamples = parsedParams.nSamples;
nCompartments = parsedParams.nCompartments;
prevClampedVoltages = parsedParams.clampedVoltages;

% Find the index of the compartment to patch
idxCompartment = xolotl_compartment_index(xolotlObject, compartment);

% Find the number of rows for prevClampedVoltage
nRowsPrev = size(prevClampedVoltages, 1);

%% Create clamped voltage(s)
% Initialize with prevClampedVoltage
newClampedVoltages = prevClampedVoltages;

% Decide on the clamped voltage for the compartment to patch
clampedVoltage = amplitude;

% Match the row count with prevClampedVoltages
clampedVoltage = match_row_count(clampedVoltage, nRowsPrev);

% Replace the corresponding column in prevClampedVoltage 
newClampedVoltages(:, idxCompartment) = clampedVoltage;

%% Replace voltage clamp with the new voltage clamp traces
xolotlObject.V_clamp = newClampedVoltages;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%{
OLD CODE:

% Create clamped voltage vector for the compartment to patch
clampedVoltage = amplitude * ones(nSamples, 1);

% Match NaNs for the other compartments
clampedVoltage = [clampedVoltage, NaN(nSamples, nCompartments - 1)];

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%