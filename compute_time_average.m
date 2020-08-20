function avgValues = compute_time_average (timeVec, valueVec, varargin)
%% Computes the time average of a value vector over time window(s)
% Usage: avgValues = compute_time_average (timeVec, valueVec, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       compute_time_average(1:5, sin(1:5))
%       compute_time_average(1:5, sin(1:5), 'IsRegular', true)
%       compute_time_average(1:100, (sin(1:100)).^2)
%       compute_time_average(1:100, (sin(1:100)).^2, 'IsRegular', true)
%
% Outputs:
%       avgValues   - average values
%                   specified as a numeric vector
%
% Arguments:
%       timeVec     - time vector
%                   must be a a numeric vector
%       valueVec    - value vector
%                   must be a a numeric vector
%       varargin    - 'TimeWindows': time windows to average
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == []
%                   - 'IsRegularSamples': whether the sampling interval 
%                                           is regular
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_subvectors.m
%       cd/find_window_endpoints.m
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%       cd/isnum.m
%
% Used by:
%       cd/plot_relative_events.m

% File History:
% 2020-08-19 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
timeWindowsDefault = [];
isRegularSamplesDefault = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'timeVec', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'valueVec', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TimeWindows', timeWindowsDefault, ...
    @(x) assert(isnum(x) || iscellnumeric(x), ...
                ['TimeWindows must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'IsRegularSamples', isRegularSamplesDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, timeVec, valueVec, varargin{:});
timeWindows = iP.Results.TimeWindows;
isRegularSamples = iP.Results.IsRegularSamples;

%% Do the job
% Find window endpoints
endPoints = find_window_endpoints(timeWindows, timeVec);

% Extract the subvectors
valueSubVecs = extract_subvectors(valueVec, 'Endpoints', endPoints);
        
% Compute the averages of the subvectors
if isRegularSamples
    % Compute the arithmetic mean
    avgValues = mean(valueSubVecs, 1);
else
    % Extract the time subvectors
    timeSubVecs = extract_subvectors(timeVec, 'Endpoints', endPoints);

    % Compute integral means
    if ndims(timeSubVecs) > 1
        % Force as column cell arrays
        [timeSubVecs, valueSubVecs] = ...
            argfun(@force_column_cell, timeSubVecs, valueSubVecs);

        % Compute integral means
        avgValues = cellfun(@compute_integral_mean, ...
                            timeSubVecs, valueSubVecs);
    else
        % Compute integral means
        avgValues = compute_integral_mean(timeSubVecs, valueSubVecs);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function avgValue = compute_integral_mean (x, y)
%% Computes the numerical approximation to the integral mean

avgValue = trapz(x, y) / (x(end) - x(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%