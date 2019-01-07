function valueRel = compute_relative_value(value, limits, varargin)
%% Computes value(s) relative to limits
% Usage: valueRel = compute_relative_value(value, limits, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       valueRel    - value(s) relative to limits
%                   specified as a TODO
% Arguments:
%       value       - value(s)
%                   must be a TODO
%       limits      - value limits
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/compute_relative_time.m
%       cd/plot_pulse_response_with_stimulus.m

% File History:
% 2018-12-28 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'value'); %, ...
    % TODO: validation function %);
addRequired(iP, 'limits'); %, ...
    % TODO: validation function %);

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, value, limits, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Compute the range of the limits
valueRange = range(limits);

% Compute the minimum limit
valueMin = min(limits);

% Compute the value relative to the minimum over the range
valueRel = (value - valueMin) ./ valueRange;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%