function [normalizedValue, initValue] = ...
                normalize_by_initial_value (newValue, initValue, varargin)
%% Normalize value(s) by the initial value(s) provided, or store initial value(s) if the latter is empty
% Usage: [normalizedValue, initValue] = ...
%               normalize_by_initial_value (newValue, initValue, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       normalizedValue - normalized value(s)
%                       specified as a numeric vector
% Arguments:
%       newValue    - new value(s) to normalize
%                   must be a numeric vector
%       initValue   - initial value(s)
%                   must be a numeric vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/match_dimensions.m
%
% Used by:    
%       cd/compute_lts_errors.m
%       cd/compute_sweep_errors.m

% File History:
% 2018-10-26 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'newValue', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'initValue', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Read from the Input Parser
parse(iP, newValue, initValue, varargin{:});

%% Preparation
% Match the dimensions of initValue to that of newValue
initValue = match_dimensions(initValue, size(newValue));

%% Do the job
if isempty(initValue);
    initValue = newValue;
    normalizedValue = ones(size(newValue));
else
    normalizedValue = newValue ./ initValue;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%