function [idxCompartment, compartment] = xolotl_compartment_index (xolotlObject, varargin)
%% Returns the index of the compartment or the default
% Usage: [idxCompartment, compartment] = xolotl_compartment_index (xolotlObject, varargin)
% Explanation:
%       TODO
% Example(s):
%       idxCompartment = xolotl_compartment_index(xolotlObject);
%       idxCompartment = xolotl_compartment_index(xolotlObject, compartment);
% Outputs:
%       idxCompartment  - index of the compartment
%                       specified as a positive integer scalar
% Arguments:
%       xolotlObject    - a created neuron with simulation parameters
%                       specified as a xolotl object
%       compartment - (opt) compartment name
%                   must be a string scalar or a character vector
%                   default == whatever the first compartment is
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/find_ind_str_in_cell.m
%
% Used by:
%       cd/xolotl_add_current_pulse.m
%       cd/xolotl_add_holding_current.m
%       cd/xolotl_add_voltage_clamp.m
%       cd/xolotl_estimate_holding_current.m

% File History:
% 2018-12-13 Created by Adam Lu
% 2018-12-13 Now returns the compartment name
% 

%% Hard-coded parameters

%% Default values for optional arguments
compartmentDefault = [];        % set later
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
addRequired(iP, 'xolotlObject');

% Add optional inputs to the Input Parser
addOptional(iP, 'compartment', compartmentDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, xolotlObject, varargin{:});
compartment = iP.Results.compartment;
% param1 = iP.Results.param1;

%% Preparation
% Extract the compartments
allCompartments = xolotlObject.Children;

%% Do the job
if isempty(compartment)
    % Match the first (in alphabetical order) compartment by default
    idxCompartment = 1;

    % Return this compartment name
    compartment = allCompartments{1};
else
    % Find the index for the compartment with given name
    idxCompartment = ...
        find_ind_str_in_cell(compartment, allCompartments, ...
                            'SearchMode', 'substrings', 'IgnoreCase', true, ...
                            'MaxNum', 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%