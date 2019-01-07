function [xLimits, xRange] = compute_xlimits (xData, varargin)
%% Computes x-axis limits from x data (works also for an x range [minX, maxX]
% Usage: [xLimits, xRange] = compute_xlimits (xData, varargin)
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
%       xData       - x data or x range
%                   must be a numeric array or a cell array of numeric arrays
%       varargin    - 'Coverage': percent coverage of x axis
%                   must be a numeric scalar between 0 and 100
%                   default == 100%
%                   - 'Separately': whether to compute x limits separately
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
%       cd/plot_pulse_response_with_stimulus.m
%       cd/plot_traces.m

% File History:
% 2018-12-17 Modified from compute_xlimits.m
% 2018-12-19 Now returns an empty array if there is no range
% 2018-12-19 Now uses apply_iteratively.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
coverageDefault = 100;      % 100% coverage of x axis by default
separatelyDefault = false;  % Compute a single set of x axis limits by default

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
addRequired(iP, 'xData', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['xData must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Coverage', coverageDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 100}));
addParameter(iP, 'Separately', separatelyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, xData, varargin{:});
coverage = iP.Results.Coverage;
separately = iP.Results.Separately;

%% Preparation
% Compute the minimum and maximum values
minValue = apply_iteratively(@min, xData);
maxValue = apply_iteratively(@max, xData);

% Check minimum and maximum values
if minValue > maxValue
    error('minimum value of xData is greater than maximum value!!');
end

%% Do the job
% Return an empty matrix if there is no range
%   Note: this makes functions like plot_traces.m not use xlim()
if minValue == maxValue
    xLimits = []
    return
end

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

%       xData       - x vector
%                   must be a numeric vector
@(x) validateattributes(x, {'numeric'}, {'vector'}));

minValue = min(xData);
maxValue = max(xData);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%