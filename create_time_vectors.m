function timeVecs = create_time_vectors (nSamples, varargin)
%% Creates time vector(s) in seconds from number(s) of samples and other optional arguments
% Usage: timeVecs = create_time_vectors (nSamples, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       timeVecs1 = create_time_vectors(10, 'SamplingRateHz', 20)
%       timeVecs2 = create_time_vectors(10, 'SamplingIntervalUs', 100)
%       timeVecs3 = create_time_vectors(6, 'BoundaryMode', 'leftadjust')
%       timeVecs4 = create_time_vectors([5, 10], 'BoundaryMode', 'span')
%       create_time_vectors(3, 'TimeUnits', 'ms', 'SamplingIntervalUs', 100)
%
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
%                   default == 'rightadjust'
%                   - 'TimeUnits': output time units
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'min'   - minutes
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
%                   - 'SamplingIntervalSeconds': sampling interval in seconds
%                   must be a positive vector
%                   default == 1 second
%                   - 'TimeStart': start time(s) in output time units
%                   must be a numeric vector
%                   default == 0
%                   - 'ForceCellOutput': whether to force output as a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   
% Requires:
%       cd/create_error_for_nargin.m
%       cd/match_format_vectors.m
%       cd/match_reciprocals.m
%
% Used by:
%       cd/atfwrite.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/create_average_time_vector.m
%       cd/create_pleth_EEG_movies.m
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/find_pulse_response_endpoints.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_xolotl_plot.m
%       cd/match_time_info.m
%       cd/parse_gas_trace.m
%       cd/parse_multiunit.m
%       cd/plot_traces_spike2_mat.m

% File History:
% 2018-10-25 Adapted from make_time_column.m in Marks_Functions
% 2018-11-28 Added 'TimeUnits' and 'SamplingIntervalMs' as optional arguments
% 2018-12-15 Added 'TimeStart' as an optional argument
% 2018-12-17 Now uses match_format_vectors.m
% 2019-01-01 Added 'ForceCellOutput' as an optional argument
% 2019-03-14 Added 'SamplingIntervalSeconds' as an optional argument
% 2019-09-11 Added 'min' as a valid time unit
% 

%% Hard-coded constants
S_PER_MIN = 60;
US_PER_S = 1e6;
MS_PER_S = 1e3;

%% Hard-coded parameters
validBoundaryModes = {'span', 'leftadjust', 'rightadjust'};
validTimeUnits = {'min', 's', 'ms', 'us'};

%% Default values for optional arguments
boundaryModeDefault = 'rightadjust';    % don't start from zero by default
timeUnitsDefault = 's';                 % return seconds by default
samplingRateHzDefault = [];             % set by match_reciprocals.m
samplingIntervalUsDefault = [];         % set by match_reciprocals.m
samplingIntervalMsDefault = [];         % set by match_reciprocals.m
samplingIntervalSecDefault = [];        % set by match_reciprocals.m
timeStartDefault = 0;                   % start at 0 by default
forceCellOutputDefault = false; % don't force output as a cell array by default

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
addParameter(iP, 'SamplingIntervalSeconds', samplingIntervalSecDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive'}));
addParameter(iP, 'TimeStart', timeStartDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'ForceCellOutput', forceCellOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, nSamples, varargin{:});
boundaryMode = validatestring(iP.Results.BoundaryMode, validBoundaryModes);
timeUnits = validatestring(iP.Results.TimeUnits, validTimeUnits);
samplingRateHz = iP.Results.SamplingRateHz;
samplingIntervalUs = iP.Results.SamplingIntervalUs;
samplingIntervalMs = iP.Results.SamplingIntervalMs;
samplingIntervalSec = iP.Results.SamplingIntervalSeconds;
tStart = iP.Results.TimeStart;
forceCellOutput = iP.Results.ForceCellOutput;

%% Preparation
% TODO: Display warning if more than one sampling rate/interval is provided
%   esp. that samplingIntervalUs overrides samplingIntervalMs, etc.

% Convert any provided sampling interval(s) to seconds
%   Note: if not provided, keep empty as default 
%           will be set with match_reciprocals.m
if ~isempty(samplingIntervalUs)
    siSeconds = samplingIntervalUs / US_PER_S;
elseif ~isempty(samplingIntervalMs)
    siSeconds = samplingIntervalMs / MS_PER_S;
elseif ~isempty(samplingIntervalSec)
    siSeconds = samplingIntervalSec;
else
    siSeconds = [];
end

% Decide on the sampling interval(s) by either using the 
%   provided sampling interval(s) or sampling rate(s)
%   Note: match_reciprocals returns 1 if both arguments are empty
%           Therefore, the default sampling interval is 1 second (1 Hz)
siSeconds = match_reciprocals(siSeconds, samplingRateHz);

% Compute the sampling interval(s) in the desired time units
switch timeUnits
    case 'min'
        siUnits = siSeconds / S_PER_MIN;
    case 's'
        siUnits = siSeconds;
    case 'ms'
        siUnits = siSeconds * MS_PER_S;
    case 'us'
        siUnits = siSeconds * US_PER_S;
    otherwise
        error(['TimeUnits %s unrecognized!\n', ...
                'Type ''help %s'' for usage'], timeUnits, mfilename);
end

% Match the formats of vectors
[nSamples, siUnits, tStart] = match_format_vectors(nSamples, siUnits, tStart);

% Count the number of vectors
nVectors = numel(nSamples);

%% Create the time vector(s)
if nVectors == 1
    timeVecs = create_time_vector(nSamples, siUnits, tStart, ...
                                        boundaryMode, mfilename);
else
    timeVecs = ...
        arrayfun(@(x, y, z) create_time_vector(x, y, z, ...
                                        boundaryMode, mfilename), ...
                nSamples, siUnits, tStart, 'UniformOutput', false);
end

%% Output
% Force as a cell array if requested
if forceCellOutput && ~iscell(timeVecs)
    timeVecs = {timeVecs};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tVec = create_time_vector (nSamples, siUnits, tStart, ...
                                    boundaryMode, functionName)
%% Create a time vector with target units

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

nVectors = length(nSamples);

% Force as column vectors
[nSamples, siUnits, tStart] = ...
    argfun(@force_column_vector, nSamples, siUnits, tStart);

% Query the number of values provided for each parameter
nRowsNSamples = length(nSamples);
nRowsSiUnits = length(siUnits);
nRowsTStart = length(tStart);

% Decide on the number of vectors to build
nVectors = max([nRowsNSamples, nRowsSiUnits, nRowsTStart]);

% Match the number of rows of each parameter to nVectors
[nSamples, siUnits, tStart] = ...
    argfun(@(x) match_row_count(x, nVectors), nSamples, siUnits, tStart);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
