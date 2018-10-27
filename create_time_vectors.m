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
%                   - 'SamplingRateHz': sampling rate in Hz
%                   must be a positive vector
%                   default == 1 Hz
%                   - 'SamplingIntervalUs': sampling interval in us
%                   must be a positive vector
%                   default == 1e6 us
%
% Requires:
%       cd/match_dimensions.m
%       cd/match_reciprocals.m
%
% Used by:    
%       cd/compute_single_neuron_errors.m

% File History:
% 2018-10-25 Adapted from make_time_column.m in Marks_Functions
% 

%% Hard-coded constants
US_PER_S = 1e6;

%% Hard-coded parameters
validBoundaryModes = {'span', 'leftadjust', 'rightadjust'};

%% Default values for optional arguments
boundaryModeDefault = 'rightadjust';    % don't start from zero by default
samplingRateHzDefault = [];             % set later
samplingIntervalUsDefault = [];         % set later

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
addParameter(iP, 'SamplingRateHz', samplingRateHzDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive'}));
addParameter(iP, 'SamplingIntervalUs', samplingIntervalUsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive'}));

% Read from the Input Parser
parse(iP, nSamples, varargin{:});
boundaryMode = validatestring(iP.Results.BoundaryMode, validBoundaryModes);
samplingRateHz = iP.Results.SamplingRateHz;
samplingIntervalUs = iP.Results.SamplingIntervalUs;

%% Preparation
% Count the number of vectors to build
nVectors = length(nSamples);

% Query the dimensions of nSamples
dimOutput = size(nSamples);

% Convert sampling interval(s) to seconds
siSeconds = samplingIntervalUs / US_PER_S;

% Decide on the sampling interval(s)
siSeconds = match_reciprocals(siSeconds, samplingRateHz);

% Match the dimensions of siSeconds to that of nSamples
siSeconds = match_dimensions(siSeconds, dimOutput);

%% Do the job
if nVectors == 1
    timeVecs = create_time_vectors_helper(nSamples, boundaryMode, siSeconds);
else
    timeVecs = ...
        arrayfun(@(x, y) create_time_vectors_helper(x, boundaryMode, y), ...
                nSamples, siSeconds, 'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tVec = create_time_vectors_helper (nSamples, boundaryMode, siSeconds)
%% Create a time vector

% Choose based on boundary mode
switch boundaryMode
    case 'span'
        % Include both endpoints
        tVec = transpose(linspace(0, nSamples * siSeconds, nSamples));
    case 'leftadjust'
        % Include the left but not right endpoint
        tVec = transpose(0:(nSamples - 1)) * siSeconds;
    case 'rightadjust'
        % Include the right but not left endpoint
        tVec = transpose(1:nSamples) * siSeconds;
    otherwise
        error(['BoundaryMode unrecognized!\n', ...
                'Type ''help %s'' for usage'], mfilename);
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