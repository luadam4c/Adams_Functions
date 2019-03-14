function h = plot_horizontal_line (yValue, varargin)
%% Plots a horizontal line
% Usage: h = plot_horizontal_line (yValue, varargin)
% Explanation:
%       TODO
% Example(s):
%       h = plot_horizontal_line(yValue)
%       h = plot_horizontal_line(yValue, 'XLimits', xLimits)
% Outputs:
%       h           - handle to the line object created
%                   specified as a primitive line object handle
% Arguments:
%       yValue      - the y value for the horizontal line
%                   must be a numeric, datetime or duration array
%       varargin    - 'XLimits': x value limits for the line
%                   must be empty or a numeric vector of 2 elements
%                       or an array of 2 columns
%                   default == get(gca, 'XLim')
%                   - Any other parameter-value pair for the line() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/isnum.m
%       cd/match_format_vector_sets.m
%
% Used by:
%       cd/parse_multiunit.m
%       cd/plot_pulse_response_with_stimulus.m

% File History:
% 2018-12-19 Created by Adam Lu
% 2018-12-27 Now allows yValue to be an array
% 2018-12-27 Now accepts datetime and duration arrays
% 2019-01-24 Now accepts multiple x limits
% 

%% Hard-coded parameters

%% Default values for optional arguments
xLimitsDefault = [];

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
addRequired(iP, 'yValue', ...
    @(x) validateattributes(x, {'numeric', 'datetime', 'duration'}, {'3d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'XLimits', xLimitsDefault);

% Read from the Input Parser
parse(iP, yValue, varargin{:});
xLimits = iP.Results.XLimits;

% Keep unmatched arguments for the line() function
otherArguments = iP.Unmatched;

%% Preparation
% Set default x value limits
if isempty(xLimits)
    xLimits = get(gca, 'XLim');
end

% Force as a cell array of column vectors and match vectors
[yValue, xLimits] = match_format_vector_sets(num2cell(yValue), xLimits);

%% Do the job
h = cellfun(@(y, x) line(x, repmat(y, size(x)), otherArguments), ...
            yValue, xLimits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

@(x) isempty(x) || isnumeric(x) && isvector(x) && length(x) == 2);

h = arrayfun(@(y) line(xLimits, repmat(y, size(xLimits)), otherArguments), ...
            yValue);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%