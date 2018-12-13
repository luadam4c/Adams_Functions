function [parsedParams] = parse_xolotl_object (xolotlObject)
%% Parses a xolotl object
% Usage: [parsedParams] = parse_xolotl_object (xolotlObject)
% Explanation:
%       TODO
% Example(s):
%       parsedParams = parse_xolotl_object(xolotlObject);
% Outputs:
%       parsedParams    - parsed parameters, including:
%                           nSamples: number of samples
%                           nCompartments: number of compartments
%                           compartments: compartment names
%                           siMs: sampling interval in ms
%                           totalDuration: total duration in ms
%                           clampedVoltages: clamped voltages
%                           externalCurrents: injected currents
%                       TODO
%                       specified as a structure
% Arguments:
%       xolotlObject    - a created neuron with simulation parameters
%                       must be a xolotl object
%
% Requires:
%       cd/match_row_count.m
%
% Used by:
%       cd/xolotl_add_current_pulse.m
%       cd/xolotl_add_holding_current.m
%       cd/xolotl_add_voltage_clamp.m
%       cd/xolotl_set_simparams.m

% File History:
% 2018-12-13 Created by Adam Lu
% 2018-12-13 Now does not match row counts here
% 

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

% Read from the Input Parser
parse(iP, xolotlObject);

%% Do the job
% Extract the compartments
% Note: This returns a row cell array of compartments in alphabetical order
compartments = xolotlObject.Children;

% Note: This returns a column cell array of compartments in alphabetical order
compartments = xolotlObject.find('compartment');

% Compute the number of compartments
nCompartments = numel(compartments);

% Extract the sampling interval in ms
siMs = xolotlObject.dt;

% Extract the total duration of the simulation in ms
totalDuration = xolotlObject.t_end;

% Extract any set voltage clamps
clampedVoltages = xolotlObject.V_clamp;

% Extract any set injected currents
externalCurrents = xolotlObject.I_ext;

% Get the number of samples
nSamples = floor(totalDuration / siMs);

%% Save in output
parsedParams.nSamples = nSamples;
parsedParams.nCompartments = nCompartments;
parsedParams.compartments = compartments;
parsedParams.siMs = siMs;
parsedParams.totalDuration = totalDuration;
parsedParams.clampedVoltages = clampedVoltages;
parsedParams.externalCurrents = externalCurrents;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Extract the number of compartments from the size of externalCurrents
nCompartments = size(externalCurrents, 2);

% Match the number of rows with nSamples
if size(externalCurrents, 1) == 1
    externalCurrents = ...
        match_row_count(externalCurrents, nSamples);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%