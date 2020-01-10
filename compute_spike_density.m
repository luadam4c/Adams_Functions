function spikeDensityHz = compute_spike_density (spikeTimes, varargin)
%% Computes the spike density from spike times and overlapping bins
% Usage: spikeDensityHz = compute_spike_density (spikeTimes, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       spikeDensityHz = compute_spike_density([1, 1, 2, 4, 8, 9, 9, 10])
%
% Outputs:
%       spikeDensityHz  - spike density in Hz
%                       specified as a numeric vector
%
% Arguments:
%       spikeTimes      - spike times
%                       must be a numeric vector
%       varargin    - 'TimeWindow': time window
%                   must be a 2-element numeric vector
%                   default == [0; max(spikeTimes)]
%                   - 'Resolution': time difference between each data point
%                   must be a positive scalar
%                   default == range(timeWindow) / nBinsDefault
%                   - 'BinWidth': bin width for spike counts
%                   must be a positive scalar
%                   default == resolution * binWidth2ResolutionDefault
%                   - 'TimeUnits': time units used
%                   must be an unambiguous, case-insensitive match to one of: 
%                       's'     - seconds
%                       'ms'    - milliseconds
%                       'us'    - microseconds
%                   default == 's'
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/parse_multiunit.m

% File History:
% 2019-05-15 Created by Adam Lu
% 2019-05-16 Now allows spikeTimes to be empty
% 

%% Hard-coded parameters
MS_PER_S = 1e3;
US_PER_S = 1e6;
validTimeUnits = {'s', 'ms', 'us'};

% TODO: Make these optional parameters
nBinsDefault = 200;
binWidth2ResolutionDefault = 2;

%% Default values for optional arguments
timeWindowDefault = [];     % set later
resolutionDefault = [];     % set later
binWidthDefault = [];       % set later
timeUnitsDefault = 'ms';    % time in ms by default

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
addRequired(iP, 'spikeTimes', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TimeWindow', timeWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'Resolution', resolutionDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'BinWidth', binWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'TimeUnits', timeUnitsDefault, ...
    @(x) any(validatestring(x, validTimeUnits)));

% Read from the Input Parser
parse(iP, spikeTimes, varargin{:});
timeWindow = iP.Results.TimeWindow;
resolution = iP.Results.Resolution;
binWidth = iP.Results.BinWidth;
timeUnits = validatestring(iP.Results.TimeUnits, validTimeUnits);

%% Preparation
% Set default time window
if isempty(timeWindow)
    timeWindow = [0; max(spikeTimes)];
end

% Set default resolution
if isempty(resolution)
    resolution = range(timeWindow) / nBinsDefault;
end

% Set default bin width
if isempty(binWidth)
    binWidth = resolution * binWidth2ResolutionDefault;
end

% Compute the desired conversion factor
switch timeUnits
    case 's'
        conversion = 1;
    case 'ms'
        conversion = 1/MS_PER_S;
    case 'us'
        conversion = 1/US_PER_S;
    otherwise
        error(['TimeUnits %s unrecognized!\n', ...
                'Type ''help %s'' for usage'], timeUnits, mfilename);
end

%% Do the job
% Compute the time window length
windowLength = range(timeWindow);

% Compute the number of bins
nBins = floor((windowLength - binWidth) / resolution) + 1;

% Compute the bin left end points
binLeft = transpose(0:nBins-1) * resolution;

% Compute the bin right end points
binRight = binLeft + binWidth;

% Compute the number of spikes in each bin
binCounts = arrayfun(@(x, y) sum(spikeTimes >= x & spikeTimes < y), ...
                        binLeft, binRight);

% Compute the spike density in Hz
spikeDensityHz = binCounts / (binWidth * conversion);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%