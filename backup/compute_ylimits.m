function [yLimits, yRange] = compute_ylimits (yData, varargin)
%% Computes y-axis limits from y data (works also for a y range [minY, maxY])
% Usage: [yLimits, yRange] = compute_ylimits (yData, varargin)
% Explanation:
%       TODO
% Example(s):
%       minY = apply_iteratively(@min, data);
%       maxY = apply_iteratively(@max, data);
%       yLimits = compute_ylimits(minY, maxY, 'Coverage', 80);
% Outputs:
%       yLimits     - computed y axis limits
%                   specified as a 2-element numeric vector
%                       or a cell array of them
%       yRange      - computed y axis range
%                   specified as a numeric vector
% Arguments:
%       yData       - y data or y range
%                   must be a numeric array or a cell array of numeric arrays
%       varargin    - 'Coverage': percent coverage of y axis
%                   must be a numeric scalar between 0 and 100
%                   default == 80%
%                   - 'Separately': whether to compute y limits separately
%                                       for each vector
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/apply_iteratively.m
%       cd/create_error_for_nargin.m
%       cd/iscellnumeric.m
%
% Used by:
%       cd/m3ha_xolotl_plot.m
%       cd/plot_cfit_pulse_response.m
%       cd/plot_pulse.m
%       cd/plot_pulse_response.m
%       cd/plot_pulse_response_with_stimulus.m
%       cd/plot_traces.m

% File History:
% 2018-10-12 Created by Adam Lu
% 2018-12-18 Now accepts vectors as arguments
% 2018-12-19 Now returns an empty array if there is no range
% 2018-12-19 Now uses apply_iteratively.m

%% Hard-coded parameters

%% Default values for optional arguments
coverageDefault = 80;       % data has 80% coverage of y axis by default
separatelyDefault = false;  % Compute a single set of y axis limits by default

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
addRequired(iP, 'yData', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['yData must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Coverage', coverageDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 100}));
addParameter(iP, 'Separately', separatelyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, yData, varargin{:});
coverage = iP.Results.Coverage;
separately = iP.Results.Separately;

% Check relationships between arguments
if any(minValue > maxValue)
    error('minValue can''t be greater than maxValue!!');
end

%% Preparation
% Compute the minimum and maximum values
minValue = apply_iteratively(@min, yData);
maxValue = apply_iteratively(@max, yData);

% Check minimum and maximum values
if minValue > maxValue
    error('minimum value of xData is greater than maximum value!!');
end

%% Do the job
% Return an empty matrix if there is no range
%   Note: this makes functions like plot_traces.m not use ylim()
if minValue == maxValue
    yLimits = []
    return
end

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

function [yLimits, yRange] = compute_ylimits (minValue, maxValue, varargin)
%% Computes y-axis limits from a minimum and maximum value

%       minValue    - minimum value of important data
%                   must be a numeric vector
%       maxValue    - maximum value of important data
%                   must be a numeric vector

addRequired(iP, 'minValue', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addRequired(iP, 'maxValue', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
parse(iP, minValue, maxValue, varargin{:});


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%