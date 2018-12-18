function [xLimits, xRange] = compute_xlimits (xVec, varargin)
%% Computes y-axis limits from an x vector (could be endpoints in actual x units)
% Usage: [xLimits, xRange] = compute_xlimits (xVec, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       xLimits     - computed x axis limits
%                   specified as a 2-element numeric vector
%       xRange      - computed x axis range
%                   specified as a numeric scalar
% Arguments:
%       xVec        - x vector
%                   must be a numeric vector
%       varargin    - 'Coverage': percent coverage of x axis
%                   must be a numeric scalar between 0 and 100
%                   default == 100%
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/plot_pulse_response_with_stimulus.m

% File History:
% 2018-12-17 Modified from compute_xlimits.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
coverageDefault = 100;      % 100% coverage of x axis by default

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
addRequired(iP, 'xVec', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Coverage', coverageDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 100}));

% Read from the Input Parser
parse(iP, xVec, varargin{:});
coverage = iP.Results.Coverage;

%% Preparation
% Compute the minimum and maximum values
minValue = min(xVec);
maxValue = max(xVec);

% Check minimum and maximum values
if minValue > maxValue
    error('minimum value of xVec is greater than maximum value!!');
end

%% Do the job
% Compute the range
xDataRange = maxValue - minValue;

% Compute the xLimits range
xRange = xDataRange / (coverage / 100);

% Compute the padding size
xPadSize = (xRange - xDataRange) / 2;

% Compute the lower and upper x limits
xLower = minValue - xPadSize;
xUpper = maxValue + xPadSize;

% Set the x axis limits
xLimits = [xLower, xUpper];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%