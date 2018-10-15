function ylimits = compute_ylimits(minValue, maxValue, varargin)
%% Computes y-axis limits from a minimum and maximum value
% Usage: ylimits = compute_ylimits(minValue, maxValue, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       ylimits     - computed y axis limits
%                   specified as a 2-element numeric vector
% Arguments:    
%       minValue    - minimum value of important data
%                   must be a numeric scalar
%       maxValue    - maximum value of important data
%                   must be a numeric scalar
%       varargin    - 'Coverage': percent coverage of y axis
%                   must be a numeric scalar between 0 and 100
%                   default == 80%
%
% Used by:
%       cd/plot_cfit_pulse_response.m
%       cd/plot_pulse.m
%       cd/plot_pulse_response.m

% File History:
% 2018-10-12 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
coverageDefault = 80;           % data has 80% coverage of y axis by default

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
if minValue > maxValue
    error('minValue can''t be greater than maxValue!!');
end

%% Do the job
% Compute the range
dataRange = maxValue - minValue;

% Compute the ylimits range
yRange = dataRange / (coverage / 100);

% Compute the padding size
yPadSize = (yRange - dataRange) / 2;

% Compute the lower and upper y limits
yLower = minValue - yPadSize;
yUpper = maxValue + yPadSize;

% Set the y axis limits
ylimits = [yLower, yUpper];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%