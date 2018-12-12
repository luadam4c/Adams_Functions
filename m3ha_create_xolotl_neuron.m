function m3ha = m3ha_create_xolotl_neuron (neuronParamsTableOrFile, varargin)
%% Creates a xolotl object for a neuron based on a parameters table
% Usage: m3ha = m3ha_create_xolotl_neuron (neuronParamsTableOrFile, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       m3ha        - the created neuron
%                   specified as a xolotl object
% Arguments:
%       neuronParamsTable   
%                   - table(s) of single neuron parameters with 
%                       parameter names as 'RowNames' and with variables:
%                       'Value': value of the parameter
%                       'LowerBound': lower bound of the parameter
%                       'UpperBound': upper bound of the parameter
%                       'JitterPercentage': jitter percentage of the parameter
%                       'IsLog': whether the parameter is 
%                                   to be varied on a log scale
%                   must be a 2d table or a cell array of 2d tables
%       neuronParamsFile
%                   - file(s) containing single neuron parameter table(s)
%                   must be a character array, a string array 
%                       or a cell array of character arrays
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/load_params.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2018-12-12 Created by Adam Lu
% 

%% Hard-coded constants
UM_PER_MM = 1e3;
UF_PER_NF = 1e-3;
CM_PER_MM = 1e-1;
OHM_PER_MOHM = 1e3;
S_PER_US = 1e-6;

%% Hard-coded parameters
valueStr = 'Value';

% Shell depth in mm
shellDepth = 1e-4;      % 0.1 um used in TC3.tem

% Concentrations in uM
caIn = 2.4e-1;          % 2.4e-4 mM used in TC3.tem
caOut = 2000;           % 2 mM used in NEURON by default

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
addRequired(iP, 'neuronParamsTableOrFile', ...
    @(x) validateattributes(x, {'table', 'cell', 'string', 'char'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, neuronParamsTableOrFile, varargin{:});
% param1 = iP.Results.param1;

% Check relationships between arguments
% TODO

%% Preparation
% Decipher the first argument

%% Load parameters from NEURON parameter files
if ischar(neuronParamsTableOrFile) || iscell(neuronParamsTableOrFile) || ...
    iscellstr(neuronParamsTableOrFile)
    % Load the parameters table
    neuronParamsTable = load_params(neuronParamsTableOrFile);
else
    % Already the parameters table
    neuronParamsTable = neuronParamsTableOrFile;
end

% Extract the geometric parameters in um
diamSoma = neuronParamsTable{'diamSoma', valueStr};
LDend = neuronParamsTable{'LDend', valueStr};
diamDend = neuronParamsTable{'diamDend', valueStr};

% Extract the specific membrane capacitance in uF/cm^2
cm = neuronParamsTable{'cm', valueStr};

% Extract the axial resistivity in Ω cm
Ra = neuronParamsTable{'Ra', valueStr};

% Extract the leak conductance in S/cm^2
gpas = neuronParamsTable{'gpas', valueStr};

% Extract the leak reversal potential in mV
epas = neuronParamsTable{'epas', valueStr};

%% Convert to xolotl units
% Compute the lengths and radii in mm
radiusSoma = (diamSoma / 2) / UM_PER_MM;
lengthSoma = diamSoma / UM_PER_MM;
radiusDendrite = (diamDend / 2) / UM_PER_MM;
lengthDendrite = LDend / UM_PER_MM;

% Compute the specific membrane capacitance in nF/mm^2
specificMembraneCapacitance = cm / (UF_PER_NF / (CM_PER_MM)^2);

% Compute the axial resistivity in MΩ mm
axialResistivity = Ra / (OHM_PER_MOHM * CM_PER_MM);

% Compute the leak conductance in uS/mm^2
gLeak = gpas / (S_PER_US / (CM_PER_MM)^2);

% The leak reversal potential is still in mV
eLeak = epas;

%% Do the job
% Create a xolotl object
m3ha = xolotl;

% Add soma
m3ha.add('compartment', 'soma', ...
            'radius', radiusSoma, 'len', lengthSoma, ...
            'Cm', specificMembraneCapacitance, ...
            'shell_thickness', shellDepth, ...
            'Ca', caIn, 'Ca_out', caOut);

% Add dend1
m3ha.add('compartment', 'dend1', ...
            'radius', radiusDendrite, 'len', lengthDendrite/2, ...
            'Cm', specificMembraneCapacitance, ...
            'shell_thickness', shellDepth, ...
            'Ca', caIn, 'Ca_out', caOut);

% Add dend2
m3ha.add('compartment', 'dend2', ...
            'radius', radiusDendrite, 'len', lengthDendrite/2, ...
            'Cm', specificMembraneCapacitance, ...
            'shell_thickness', shellDepth, ...
            'Ca', caIn, 'Ca_out', caOut);

% Connect the compartments
% TODO: set axial resistivity
m3ha.connect('soma', 'dend1');
% m3ha.connect('soma', 'dend1', 'gbar', axialConductance);
m3ha.connect('dend1', 'dend2');

% Add passive conductances
m3ha.soma.add('conductance', 'Leak', 'gbar', gLeak, 'E', eLeak);
m3ha.dend1.add('conductance', 'Leak', 'gbar', gLeak, 'E', eLeak);
m3ha.dend2.add('conductance', 'Leak', 'gbar', gLeak, 'E', eLeak);

%% Output results
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

neuronParamsFile = '/media/adamX/m3ha/optimizer4gabab/initial_params/initial_params_D091710.csv';

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%