function xolotlObject = m3ha_xolotl_create_neuron (neuronParamsTableOrFile, varargin)
%% Creates a xolotl object for a 3-compartment neuron based on a parameters table
% Usage: xolotlObject = m3ha_xolotl_create_neuron (neuronParamsTableOrFile, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       xolotlObject    - a created neuron
%                       specified as a xolotl object
% Arguments:
%       neuronParamsTable   
%                   - table(s) of single neuron parameters with 
%                       parameter names as row names and with 
%                       at least the variable:
%                       'Value': value of the parameter
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
%       cd/create_error_for_nargin.m
%       cd/load_params.m
%       cd/print_cellstr.m
%
% Used by:
%       cd/m3ha_xolotl_test.m

% File History:
% 2018-12-12 Created by Adam Lu
% 2019-08-08 Added back tree_idx
%   TODO: Implement a cell array of tables or files
% 

%% Hard-coded constants
UM_PER_MM = 1e3;
UF_PER_NF = 1e-3;
CM_PER_MM = 1e-1;
OHM_PER_MOHM = 1e3;
S_PER_US = 1e-6;

%% Hard-coded parameters
valueStr = 'Value';
diamSomaStr = 'diamSoma';
LDendStr = 'LDend';
diamDendStr = 'diamDend';
cmStr = 'cm';
RaStr = 'Ra';
gpasStr = 'gpas';
epasStr = 'epas';

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
    error(create_error_for_nargin(mfilename));
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

%% Preparation
% Decipher the first argument

%% Load parameters from NEURON parameter files
if ischar(neuronParamsTableOrFile) || isstring(neuronParamsTableOrFile) || ...
    iscellstr(neuronParamsTableOrFile)
    % Check if the file exists
    if ~isfile(neuronParamsTableOrFile)
        fprintf('The file %s does not exist!!\n', neuronParamsTableOrFile);
        xolotlObject = [];
        return
    end

    % Load the parameters table
    neuronParamsTable = load_params(neuronParamsTableOrFile);
else
    % Already the parameters table
    neuronParamsTable = neuronParamsTableOrFile;
end

% Extract the geometric parameters in um
diamSoma = neuronParamsTable{diamSomaStr, valueStr};
LDend = neuronParamsTable{LDendStr, valueStr};
diamDend = neuronParamsTable{diamDendStr, valueStr};

% Extract the specific membrane capacitance in uF/cm^2
cm = neuronParamsTable{cmStr, valueStr};

% Extract the axial resistivity in Ω cm
Ra = neuronParamsTable{RaStr, valueStr};

% Extract the leak conductance in S/cm^2
gpas = neuronParamsTable{gpasStr, valueStr};

% Extract the leak reversal potential in mV
epas = neuronParamsTable{epasStr, valueStr};

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
xolotlObject = xolotl;

% Add soma
%   Note: tree_idx must be set for xolotl to consider
%           the model as multi-compartment
xolotlObject.add('compartment', 'soma', ...
            'radius', radiusSoma, 'len', lengthSoma, ...
            'Cm', specificMembraneCapacitance, ...
            'shell_thickness', shellDepth, ...
            'Ca', caIn, 'Ca_out', caOut, ...
            'tree_idx', 0);

% Add dend1
%   Note: tree_idx does not have to be set here as long as
%           other compartments are connected with 'Axial' synapses
xolotlObject.add('compartment', 'dend1', ...
            'radius', radiusDendrite, 'len', lengthDendrite/2, ...
            'Cm', specificMembraneCapacitance, ...
            'shell_thickness', shellDepth, ...
            'Ca', caIn, 'Ca_out', caOut);

% Add dend2
xolotlObject.add('compartment', 'dend2', ...
            'radius', radiusDendrite, 'len', lengthDendrite/2, ...
            'Cm', specificMembraneCapacitance, ...
            'shell_thickness', shellDepth, ...
            'Ca', caIn, 'Ca_out', caOut);

% Connect the soma with dend1, dend1 with dend2
%   Note: this connects the compartment with 'Axial' synapses
xolotlObject.connect('soma', 'dend1', axialResistivity);
xolotlObject.connect('dend1', 'dend2', axialResistivity);

% Add passive conductances
xolotlObject.soma.add('Leak', 'gbar', gLeak, 'E', eLeak);
xolotlObject.dend1.add('Leak', 'gbar', gLeak, 'E', eLeak);
xolotlObject.dend2.add('Leak', 'gbar', gLeak, 'E', eLeak);

%% Print to standard output
% This returns a column cell array of compartments in alphabetical order
compartments = xolotlObject.find('compartment');

% Print the
disp('These compartments have been made:');
print_cellstr(compartments, 'OmitBraces', true, 'Delimiter', ',');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

xolotlObject.connect('soma', 'dend1', 'gbar', axialConductance);

xolotlObject.soma.add('conductance', 'Leak', 'gbar', gLeak, 'E', eLeak);
xolotlObject.dend1.add('conductance', 'Leak', 'gbar', gLeak, 'E', eLeak);
xolotlObject.dend2.add('conductance', 'Leak', 'gbar', gLeak, 'E', eLeak);

xolotlObject.connect('soma', 'dend1', 'Axial', ...
                        'resistivity', axialResistivity);
xolotlObject.connect('dend1', 'dend2', 'Axial', ...
                        'resistivity', axialResistivity);
xolotlObject.connect('dend1', 'soma', 'Axial', ...
                        'resistivity', axialResistivity);
xolotlObject.connect('dend2', 'dend1', 'Axial', ...
                        'resistivity', axialResistivity);

% Add dend
xolotlObject.add('compartment', 'dend', ...
            'radius', radiusDendrite, 'len', lengthDendrite, ...
            'Cm', specificMembraneCapacitance, ...
            'shell_thickness', shellDepth, ...
            'Ca', caIn, 'Ca_out', caOut);

% Slice the dendrite into two compartments
xolotlObject.slice('dend', 2);
xolotlObject.connect('soma', 'dend1');

xolotlObject.connect('dend1', 'soma', axialResistivity);
xolotlObject.connect('dend2', 'dend1', axialResistivity);


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%