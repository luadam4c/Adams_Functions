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
%       varargin    - 'Compartment': the compartment name to add clamp
%                   must be a string scalar or a character vector
%                   default == set in xolotl_compartment_index.m
%                   - 'Amplitude': amplitude in nA 
%                   must be a numeric scalar
%                   default == 0 nA
%
% Requires:
%       cd/parse_xolotl_object.m
%       cd/xolotl_compartment_index.m
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
compartmentDefault = '';    % set later
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
previousCurrentInjections = parsedParams.externalCurrents;

% Find the index of the compartment to patch
idxCompartment = xolotl_compartment_index(xolotlObject, compartment);

% Find the number of rows for previousCurrentInjections
nRowsPrev = size(previousCurrentInjections, 1);

%% Create holding current(s)
% Initialize as zeros
newHoldingCurrents = zeros(nRowsPrev, nCompartments);

% Decide on a holding current for the compartment to patch
holdingCurrent = amplitude;

% Match the row count with previousCurrentInjections
holdingCurrent = match_row_count(holdingCurrent, nRowsPrev);

% Replace the corresponding column in prevClampedVoltage 
newHoldingCurrents(:, idxCompartment) = holdingCurrent;

%% Add the holding current vectors to the previously set current injections
xolotlObject.I_ext = previousCurrentInjections + newHoldingCurrents;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%{
OLD CODE:

% Create holding current vector for the compartment to patch
holdingCurrent = amplitude * ones(nSamples, 1);

% Match zeros for the other compartments
holdingCurrent = [holdingCurrent, zeros(nSamples, nCompartments - 1)];

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%