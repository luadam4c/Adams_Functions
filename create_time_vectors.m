function timeVecs = create_time_vectors (nSamples, varargin)
%% Creates time vector(s) in seconds from number(s) of samples and other optional arguments
% Usage: timeVecs = create_time_vectors (nSamples, varargin)
% Explanation:
%       TODO
% Example(s):
%       timeVecs1 = create_time_vectors(10, 'SamplingRateHz', 20)
%       timeVecs2 = create_time_vectors(10, 'SamplingIntervalUs', 100)
%       timeVecs3 = create_time_vectors(6, 'BoundaryMode', 'leftadjust')
%       timeVecs4 = create_time_vectors([5, 10], 'BoundaryMode', 'span')
%       create_time_vectors(3, 'TimeUnits', 'ms', 'SamplingIntervalUs', 100)
% Outputs:
%       timeVecs    - created time vectors(s)
%                   specified as a column vector
%                       or a cell array of column vectors
% Arguments:    
%       nSamples    - number of sample points for each vector to build
%                   must be a positive integer vector
%       varargin    - 'BoundaryMode': boundary mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'span'        - both endpoints are included
%                       'leftadjust'  - the left endpoint is included,
%                                           but not the right endpoint
%                       'rightadjust' - the right endpoint is included
%                                           but not the left endpoint
%                   default == 'restrictive'
%                   - 'TimeUnits': output time units
%                   must be an unambiguous, case-insensitive match to one of: 
%                       's'     - seconds
%                       'ms'    - milliseconds
%                       'us'    - microseconds
%                   default == 's'
%                   - 'SamplingRateHz': sampling rate in Hz
%                   must be a positive vector
%                   default == 1 Hz
%                   - 'SamplingIntervalUs': sampling interval in us
%                   must be a positive vector
%                   default == 1e6 us (1 second)
%                   - 'SamplingIntervalMs': sampling interval in ms
%                   must be a positive vector
%                   default == 1e3 ms (1 second)
%                   - 'TimeStart': start time(s) in ms
%                   must be a positive vector
%                   default == 0
%
% Requires:
%       cd/match_dimensions.m
%       cd/match_reciprocals.m
%
% Used by:
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/create_average_time_vector.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_xolotl_plot.m

% File History:
% 2018-10-25 Adapted from make_time_column.m in Marks_Functions
% 2018-11-28 Added 'TimeUnits' and 'SamplingIntervalMs' as optional arguments
% 2018-12-15 Added 'TimeStart' as an optional argument
% 

%% Hard-coded constants
US_PER_S = 1e6;
MS_PER_S = 1e3;

%% Hard-coded parameters
validBoundaryModes = {'span', 'leftadjust', 'rightadjust'};
validTimeUnits = {'s', 'ms', 'us'};

%% Default values for optional arguments
boundaryModeDefault = 'rightadjust';    % don't start from zero by default
timeUnitsDefault = 's';                 % return seconds by default
samplingRateHzDefault = [];             % set by match_reciprocals.m
samplingIntervalUsDefault = [];         % set by match_reciprocals.m
samplingIntervalMsDefault = [];         % set by match_reciprocals.m
timeStartDefault = 0;                   % start at 0 by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'nSamples', ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive', 'integer'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BoundaryMode', boundaryModeDefault, ...
    @(x) any(validatestring(x, validBoundaryModes)));
addParameter(iP, 'TimeUnits', timeUnitsDefault, ...
    @(x) any(validatestring(x, validTimeUnits)));
addParameter(iP, 'SamplingRateHz', samplingRateHzDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive'}));
addParameter(iP, 'SamplingIntervalUs', samplingIntervalUsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive'}));
addParameter(iP, 'SamplingIntervalMs', samplingIntervalMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive'}));
addParameter(iP, 'TimeStart', timeStartDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Read from the Input Parser
parse(iP, nSamples, varargin{:});
boundaryMode = validatestring(iP.Results.BoundaryMode, validBoundaryModes);
timeUnits = validatestring(iP.Results.TimeUnits, validTimeUnits);
samplingRateHz = iP.Results.SamplingRateHz;
samplingIntervalUs = iP.Results.SamplingIntervalUs;
samplingIntervalMs = iP.Results.SamplingIntervalMs;
tStart = iP.Results.TimeStart;

%% Preparation
% Count the number of vectors to build
nVectors = length(nSamples);

% Query the dimensions of nSamples
dimOutput = size(nSamples);

% Convert sampling interval(s) to seconds
if ~isempty(samplingIntervalUs)
    siSeconds = samplingIntervalUs / US_PER_S;
elseif ~isempty(samplingIntervalMs)
    siSeconds = samplingIntervalMs / MS_PER_S;

    % TODO: Display warning if samplingIntervalUs also provided
else
    siSeconds = [];
end

% Decide on the sampling interval(s)
%   Note: match_reciprocals returns 1 if both arguments are empty
siSeconds = match_reciprocals(siSeconds, samplingRateHz);

% Match the dimensions of siSeconds to that of nSamples
siSeconds = match_dimensions(siSeconds, dimOutput);

% Match the dimensions of tStart to that of nSamples
tStart = match_dimensions(tStart, dimOutput);

%% Do the job
if nVectors == 1
    timeVecs = create_time_vectors_helper(nSamples, siSeconds, tStart, ...
                                        boundaryMode, timeUnits, mfilename);
else
    timeVecs = ...
        arrayfun(@(x, y, z) create_time_vectors_helper(x, y, z, ...
                                        boundaryMode, timeUnits, mfilename), ...
                nSamples, siSeconds, tStart, 'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tVec = create_time_vectors_helper (nSamples, siSeconds, tStart, ...
                                    boundaryMode, timeUnits, functionName)
%% Create a time vector with target units

%% Hard-coded constants
US_PER_S = 1e6;
MS_PER_S = 1e3;

% Compute the sampling interval in the correct units
switch timeUnits
    case 's'
        siUnits = siSeconds;
    case 'ms'
        siUnits = siSeconds * MS_PER_S;
    case 'us'
        siUnits = siSeconds * US_PER_S;
    otherwise
        error(['TimeUnits %s unrecognized!\n', ...
                'Type ''help %s'' for usage'], timeUnits, functionName);
end

% Create the time vector based the boundary mode
switch boundaryMode
    case 'span'
        % Include both endpoints
        tVec = tStart + transpose(linspace(0, nSamples * siUnits, nSamples));
    case 'leftadjust'
        % Include the left but not right endpoint
        tVec = tStart + transpose(0:(nSamples - 1)) * siUnits;
    case 'rightadjust'
        % Include the right but not left endpoint
        tVec = tStart + transpose(1:nSamples) * siUnits;
    otherwise
        error(['BoundaryMode %s unrecognized!\n', ...
                'Type ''help %s'' for usage'], boundaryMode, functionName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Match sampling rates and intervals to those of nSamples
[samplingRateHz, siSeconds] = ...
    argfun(@(x) match_dimensions_if_not_empty(x, dimOutput), ...
                samplingRateHz, siSeconds);


function x = match_dimensions_if_not_empty(x, varargin)
%% Apply match_dimensions only if first argument not empty

if ~isempty(x)
    x = match_dimensions(x, varargin{:});
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%