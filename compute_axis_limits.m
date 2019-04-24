function [limits, axisRange] = compute_axis_limits (dataOrRange, axisType, varargin)
%% Computes x or y axis limits from data (works also for a range [min(data), max(data)])
% Usage: [limits, axisRange] = compute_axis_limits (dataOrRange, axisType, varargin)
% Explanation:
%       Computes axis limits from data by using the entire range 
%           for the x axis and expanding by 10% on each side for the y axis
%           (the default 'Coverage' parameter is 100% or 80%, respectively)
%       If 'Separately' is set to true and there are multiple data vectors
%           or ranges, a separate axis limits is computed for each vector
%           (each element of a cell array or each column of a numeric array)
%           and returned as a cell array
%       
% Example(s):
%       limits = compute_axis_limits(xData, 'x');
%       limits = compute_axis_limits(xRange, 'x');
%       limits = compute_axis_limits(yData, 'y', 'Coverage', 70);
% Outputs:
%       limits     - computed y axis limits
%                   specified as a 2-element numeric vector
%                       or a cell array of them
%       axisRange   - computed axis range
%                   specified as a numeric vector
% Arguments:
%       dataOrRange - data for this axis or range along this axis
%                   must be a numeric array or a cell array of numeric arrays
%       axisType    - axis type
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'x'- x axis
%                       'y'- y axis
%       varargin    - 'Coverage': percent coverage of axis
%                   must be empty or a numeric scalar between 0 and 100
%                   default == 100% for x axis and 80% for y axis
%                   - 'Separately': whether to compute axis limits separately
%                                       for each vector
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'AutoZoom': whether to zoom in on the axis 
%                                   to within 3 SDs of the mean of all values
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/apply_iteratively.m
%       cd/create_error_for_nargin.m
%       cd/isnum.m
%
% Used by:
%       cd/m3ha_xolotl_plot.m
%       cd/parse_multiunit.m
%       cd/plot_cfit_pulse_response.m
%       cd/plot_pulse.m
%       cd/plot_pulse_response.m
%       cd/plot_pulse_response_with_stimulus.m
%       cd/plot_traces.m

% File History:
% 2018-12-19 Combined compute_xlimits.m and compute_ylimits.m
% 2019-04-24 Added 'AutoZoom' as an optional argument

%% Hard-coded parameters
validAxisTypes = {'x', 'y'};
xCoverageDefault = 100;     % 100% coverage of x axis by default
yCoverageDefault = 80;      % 80% coverage of x axis by default
maxNStdBeyondMean = 3;      % maximum number of standard deviations 
                            %   beyond the mean

%% Default values for optional arguments
coverageDefault = [];       % set later
separatelyDefault = false;  % compute a single set of axis limits by default
autoZoomDefault = false;    % don't zoom in on axis by default

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
addRequired(iP, 'dataOrRange', ...
    @(x) assert(isnum(x) || iscell(x), ...
                'dataOrRange must be either a numeric array or a cell array!'));
addRequired(iP, 'axisType', ...
    @(x) any(validatestring(x, validAxisTypes)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Coverage', coverageDefault, ...
    @(x) isempty(x) || isnumeric(x) && isscalar(x) && x >= 0 && x <= 100);
addParameter(iP, 'Separately', separatelyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'AutoZoom', autoZoomDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, dataOrRange, axisType, varargin{:});
coverage = iP.Results.Coverage;
separately = iP.Results.Separately;
autoZoom = iP.Results.AutoZoom;

% Validate axisType
axisType = validatestring(axisType, validAxisTypes);

%% Preparation
% Set default coverage
if isempty(coverage)
    switch axisType
        case 'x'
            coverage = xCoverageDefault;
        case 'y'
            coverage = yCoverageDefault;
        otherwise
            error('Code logic error!');
    end
end

% Compute the minimum and maximum values
if autoZoom
    % Compute the mean and standard deviation of all values
    if separately
        if iscell(dataOrRange)
            meanValue = cellfun(@mean, dataOrRange);
            stdValue = cellfun(@std, dataOrRange);
        elseif isnumeric(dataOrRange)
            meanValue = mean(dataOrRange);
            stdValue = std(dataOrRange);
        end
    else
        dataLin = force_column_vector(dataOrRange, 'ToLinearize', true);
        meanValue = mean(dataLin);
        stdValue = std(dataLin);
    end

    % Compute the "minimum" and "maximum" values
    %   based on the mean and standard deviation of all values
    minValue = meanValue - maxNStdBeyondMean * stdValue;
    maxValue = meanValue + maxNStdBeyondMean * stdValue;
else
    % Compute the minimum and maximum values
    if separately
        if iscell(dataOrRange)
            minValue = cellfun(@min, dataOrRange);
            maxValue = cellfun(@max, dataOrRange);
        elseif isnumeric(dataOrRange)
            minValue = min(dataOrRange, [], 1);
            maxValue = max(dataOrRange, [], 1);
        end
    else
        minValue = apply_iteratively(@min, dataOrRange);
        maxValue = apply_iteratively(@max, dataOrRange);
    end
end

% Check minimum and maximum values
if any(minValue > maxValue)
    error('minimum value of axis data is greater than maximum value!!');
end

%% Do the job
% Return an empty matrix if there is no range
%   Note: this makes functions like plot_traces.m not use xlim() or ylim()
if minValue == maxValue
    limits = [];
    return
end

% Compute the data range
dataRange = maxValue - minValue;

% Compute the limits range
axisRange = dataRange ./ (coverage ./ 100);

% Compute the padding size
padSize = (axisRange - dataRange) / 2;

% Compute the lower and upper axis limits
lowerLimit = minValue - padSize;
upperLimit = maxValue + padSize;

% Set the y axis limits
if isscalar(lowerLimit) && isscalar(upperLimit)
    limits = [lowerLimit, upperLimit];
else
    limits = arrayfun(@(x, y) [x, y], lowerLimit, upperLimit, ...
                        'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%