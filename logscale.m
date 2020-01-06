function values = logscale (base, target, scaleFactors, varargin)
%% Creates scaled values between base and target based on a log scale
% Usage: values = logscale (base, target, scaleFactors, varargin)
% Explanation:
%       Computes a vectors of values that are spaced on a log scale
%           where the base value is scaleFactor == 0
%           and the target value is scaleFactor == 1
%       For example, logscale(1, 2, 0:0.5:1) == [1, sqrt(2), 2]
%
% Example(s):
%       logscale(1, 2, -1:5)
%       logscale(5, 12, 0:0.1:1)
%
% Outputs:
%       values     - TODO: Description of values
%                   specified as a TODO
%
% Arguments:
%       reqarg1     - TODO: Description of reqarg1
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/linscale.m
%
% Used by:
%       cd/m3ha_compute_gabab_ipsc.m

% File History:
% 2020-01-04 Created by Adam Lu
% TODO for SHINSHIN: Make linscale.m and use it here

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'base');
addRequired(iP, 'target');
addRequired(iP, 'scaleFactors');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, base, target, scaleFactors, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
values = exp(linscale(log(base), log(target), scaleFactors));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
