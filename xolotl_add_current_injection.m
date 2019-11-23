function xolotlObject = xolotl_add_current_injection (xolotlObject, varargin)
%% Adds a current injection to a xolotl object, just the first compartment by default
% Usage: xolotlObject = xolotl_add_current_injection (xolotlObject, varargin)
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
%       varargin    - 'Compartment': the compartment name to add clamp
%                   must be a string scalar or a character vector
%                   default == set in xolotl_compartment_index.m
%                   - 'CurrentVector': current vector(s)
%                   must be a numeric array
%                   default == []
%                   - 'InjectionType': TODO
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/argfun.m
%       cd/match_row_count.m
%       cd/parse_xolotl_object.m
%       cd/xolotl_compartment_index.m
%
% Used by:
%       cd/m3ha_xolotl_test.m
%       cd/xolotl_add_current_pulse.m
%       cd/xolotl_add_holding_current.m

% File History:
% 2019-08-15 Modified from xolotl_add_current_pulse.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
compartmentDefault = '';    % set later
currentVectorDefault = [];  % set later

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

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Compartment', compartmentDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'CurrentVector', currentVectorDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, xolotlObject, varargin{:});
compartment = iP.Results.Compartment;
currentVector = iP.Results.CurrentVector;

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

% Find the number of rows needed to accommodate both previous and new currents
nRowsNeeded = max(nRowsPrev, nSamples);

% Match the row count of current arrays with nRowsNeeded
[previousCurrentInjections, currentVector] = ...
    argfun(@(x) match_current_injection(x, nRowsNeeded), ...
            previousCurrentInjections, currentVector);

%% Add the current vector to the previously set current injections
% Initialize as zeros
newCurrentInjections = zeros(nRowsNeeded, nCompartments);

% Replace the corresponding column in prevClampedVoltage 
newCurrentInjections(:, idxCompartment) = currentVector;

% Set new external current
xolotlObject.I_ext = previousCurrentInjections + newCurrentInjections;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

function currents = match_current_injection (currents, nRowsNeeded)

% Decide on the expansion method
if isscalar(currents)
    expansionMethod = 'repeat';
else
    expansionMethod = 'patchZeros';
end

% Match the row counts
currents = match_row_count(currents, nRowsNeeded, ...
                            'ExpansionMethod', expansionMethod);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%