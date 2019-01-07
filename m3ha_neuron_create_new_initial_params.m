function newTable = m3ha_neuron_create_new_initial_params (prevTable, varargin)
%% Creates a new set of NEURON parameters based on information in the previous parameters table
% Usage: newTable = m3ha_neuron_create_new_initial_params (prevTable, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       newTable    - the new NEURON parameters table
%                   specified as a TODO
% Arguments:
%       prevTable   - the previous NEURON parameters table
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/choose_random_values.m
%
% Used by:
%       ~/m3ha/optimizer4gabab/optimizer_4compgabab.m

% File History:
% 2018-12-11 Created by Adam Lu
% 

%% Hard-coded parameters
valueStr = 'Value';
initValueStr = 'InitValue';
upperBoundStr = 'UpperBound';
lowerBoundStr = 'LowerBound';
isLogStr = 'IsLog';
inUseStr = 'InUse';

%% Default values for optional arguments
% param1Default   = [];                   % default TODO: Description of param1

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
addRequired(iP, 'prevTable', ...
    @(x) validateattributes(x, {'table', 'cell'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, prevTable, varargin{:});
% param1 = iP.Results.param1;

% Check relationships between arguments
% TODO

%% Choose new parameters
% Initialize with the previous set of parameters
newTable = prevTable;

% Find the row indices of the parameters that need to be changed
indInUse = find(prevTable{:, inUseStr});

% Get the names of all parameters in use
paramsInUseNames = prevTable.Properties.RowNames(indInUse);

% Retrieve info for parameters in use
paramsInUseIsLog = newTable{indInUse, isLogStr};
paramsInUseUpperBound = newTable{indInUse, upperBoundStr};
paramsInUseLowerBound = newTable{indInUse, lowerBoundStr};

% Choose new values for the parameters in use
paramsInUseNewValue = choose_random_values(paramsInUseLowerBound, ...
                        paramsInUseUpperBound, 'IsLog', paramsInUseIsLog);

% Update the neuronParamsTable with new initial values
newTable{indInUse, valueStr} = paramsInUseNewValue;
newTable{indInUse, initValueStr} = paramsInUseNewValue;

%% Constrain parameters
%{
% Constrain gpas to match RinEstimate
% TODO: make a function m3ha_constrain_neuronparams.m

% Extract the new geometric parameters
diamSoma = newTable{'diamSoma', 'Value'};
LDend = newTable{'LDend', 'Value'};
diamDend = newTable{'diamDend', 'Value'};

% Compute the surface area for each cell
surfaceArea = compute_surface_area([diamSoma, LDend], [diamSoma, diamDend]);

% Estimate the passive conductances from the input resistances
gpas = compute_gpas(RinEstimate, surfaceArea);

% Use this estimate for gpas
newTable{'gpas', 'Value'} = gpas;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
