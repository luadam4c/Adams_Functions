function h = plot_vertical_line (xValue, varargin)
%% Plots vertical line(s)
% Usage: h = plot_vertical_line (xValue, varargin)
% Explanation:
%       TODO
% Example(s):
%       h = plot_vertical_line(xValue)
%       h = plot_vertical_line(xValue, 'YLimits', yLimits)
% Outputs:
%       h           - handle to the line object(s) created
%                   specified as a primitive line object handle
% Arguments:
%       xValue      - the x value(s) for the vertical line(s)
%                   must be a numeric, datetime or duration array
%       varargin    - 'YLimits': y value limits for the line(s)
%                   must be empty or a numeric vector of 2 elements
%                       or an array of 2 rows
%                   default == get(gca, 'YLim')
%                   - Any other parameter-value pair for the line() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/isnum.m
%       cd/match_format_vector_sets.m
%
% Used by:
%       cd/plot_error_bar.m
%       cd/plot_pulse_response_with_stimulus.m
%       cd/plot_struct.m
%       cd/plot_swd_histogram.m
%       cd/plot_window_boundaries.m

% File History:
% 2018-12-19 Created by Adam Lu
% 2018-12-27 Now allows xValue to be an array
% 2018-12-27 Now accepts datetime and duration arrays
% 2019-01-24 Now accepts multiple y limits
% 

%% Hard-coded parameters

%% Default values for optional arguments
yLimitsDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'xValue', ...
    @(x) validateattributes(x, {'numeric', 'datetime', 'duration'}, {'3d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'YLimits', yLimitsDefault);

% Read from the Input Parser
parse(iP, xValue, varargin{:});
yLimits = iP.Results.YLimits;

% Keep unmatched arguments for the line() function
otherArguments = iP.Unmatched;

%% Preparation
% Set default y value limits
if isempty(yLimits)
    yLimits = get(gca, 'YLim');
end

% Force as a cell array of column vectors and match vectors
[xValue, yLimits] = match_format_vector_sets(num2cell(xValue), yLimits);

%% Do the job
h = cellfun(@(x, y) line(repmat(x, size(y)), y, otherArguments), ...
            xValue, yLimits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

h = arrayfun(@(x) line(repmat(x, size(yLimits)), yLimits, otherArguments), ...
            xValue);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || isnum(x) && isvector(x) && length(x) == 2);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%