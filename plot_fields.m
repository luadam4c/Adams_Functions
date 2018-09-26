function h = plot_fields (structure, varargin)
%% TODO: A summary of what the function does (must be a single unbreaked line)
% Usage: h = plot_fields (structure, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       h           - TODO: Description of h
%                   specified as a TODO
% Arguments:    
%       structure     - TODO: Description of structure
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/plot_tuning_curve.m
%
% Used by:    
%       /TODO:dir/TODO:file

% File History:
% 2018-09-26 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default   = [];                   % default TODO: Description of param1

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
addRequired(iP, 'structure', ...                  % TODO: Description of structure
    % TODO: validation function %);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default, ...
    % TODO: validation function %);

% Read from the Input Parser
parse(iP, structure, varargin{:});
param1 = iP.Results.param1;

% Check relationships between arguments
% TODO

%% Preparation
% Get all the fields of the structure as a cell array
allFields = fieldnames(structure);

% 

%% Plot all fields


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

