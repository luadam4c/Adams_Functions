function [yLimits, yRange] = compute_ylimits (minValue, maxValue, varargin)
%% Computes y-axis limits from a minimum and maximum value
% Usage: [yLimits, yRange] = compute_ylimits (minValue, maxValue, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       yLimits     - computed y axis limits
%                   specified as a 2-element numeric vector
%                       or a cell array of them
%       yRange      - computed y axis range
%                   specified as a numeric vector
% Arguments:    
%       minValue    - minimum value of important data
%                   must be a numeric vector
%       maxValue    - maximum value of important data
%                   must be a numeric vector
%       varargin    - 'Coverage': percent coverage of y axis
%                   must be a numeric scalar between 0 and 100
%                   default == 80%
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/plot_cfit_pulse_response.m
%       cd/plot_pulse.m
%       cd/plot_pulse_response.m
%       cd/plot_pulse_response_with_stimulus.m

% File History:
% 2018-10-12 Created by Adam Lu
% 2018-12-18 Now accepts vectors as arguments
% 

%% Hard-coded parameters

%% Default values for optional arguments
coverageDefault = 80;           % data has 80% coverage of y axis by default

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
addRequired(iP, 'minValue', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addRequired(iP, 'maxValue', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Coverage', coverageDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 100}));

% Read from the Input Parser
parse(iP, minValue, maxValue, varargin{:});
coverage = iP.Results.Coverage;

% Check relationships between arguments
if any(minValue > maxValue)
    error('minValue can''t be greater than maxValue!!');
end

%% Do the job
% Compute the range
dataRange = maxValue - minValue;

% Compute the yLimits range
yRange = dataRange ./ (coverage ./ 100);

% Compute the padding size
yPadSize = (yRange - dataRange) / 2;

% Compute the lower and upper y limits
yLower = minValue - yPadSize;
yUpper = maxValue + yPadSize;

% Set the y axis limits
if isscalar(yLower) && isscalar(yUpper)
    yLimits = [yLower, yUpper];
else
    yLimits = arrayfun(@(x, y) [x, y], yLower, yUpper, 'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%