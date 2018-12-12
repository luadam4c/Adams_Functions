%function [output1] = m3ha_create_xolotl_neuron (neuronParamsTable, varargin)
%% TODO: A summary of what the function does (must be a single unbreaked line)
% Usage: [output1] = m3ha_create_xolotl_neuron (neuronParamsTable, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
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
% 201X-XX-XX Created by TODO or Adapted from TODO
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default   = [];                   % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Deal with arguments
% % Check number of required arguments
% if nargin < 1    % TODO: 1 might need to be changed
%     error(['Not enough input arguments, ', ...
%             'type ''help %s'' for usage'], mfilename);
% end

% % Set up Input Parser Scheme
% iP = inputParser;
% iP.FunctionName = mfilename;

% % Add required inputs to the Input Parser
% addRequired(iP, 'neuronParamsTable', ...                  % TODO: Description of neuronParamsTable
%     % TODO: validation function %);

% % Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% % Read from the Input Parser
% parse(iP, neuronParamsTable, varargin{:});
% param1 = iP.Results.param1;

% Check relationships between arguments
% TODO

%% Preparation
% TODO

neuronParamsTable = load_params('')

%% Do the job
% Create a xolotl object
m3ha = xolotl;

% 
m3ha.add('compartment', 'soma', 'radius', .01);

%% Output results
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%