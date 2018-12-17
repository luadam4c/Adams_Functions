function [parsedParams, parsedData] = ...
                parse_pulse_response (vectors, siMs, varargin)
%% Parses pulse response widths, endpoints, amplitudes for vector(s) containing a pulse response
% Usage: [parsedParams, parsedData] = ...
%               parse_pulse_response (vectors, siMs, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       parsedParams    - a table containing the parsed parameters, 
%                           each row corresponding to a vector, with fields:
%                           nSamples
%                           responseWidthSamples
%                           responseWidthMs
%                           nSamplesRising
%                           nSamplesFalling
%                           idxResponseStart
%                           idxResponseEnd
%                           idxResponseMid
%                           idxBaseStart
%                           idxBaseEnd
%                           idxSteadyStart
%                           idxSteadyEnd
%                           baseValue
%                           steadyValue
%                           steadyAmplitude
%                           minValue
%                           maxValue
%                           idxPeakSearchBegin
%                           idxPeak
%                           peakValue
%                           peakAmplitude
%                           hasJump
%                           siMs
%                       specified as a table
%       parsedData      - a table containing the parsed data, with fields:
%                           vectors
%                           indBase
%                           indSteady
%                           indRising
%                           indFalling
%                           tvecsRising
%                           vvecsRising
%                           tvecsFalling
%                           vvecsFalling
%                       specified as a table
% Arguments:
%       vectors     - vectors containing a pulse response
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%       siMs        - sampling interval in ms
%                   must be a positive vector
%       varargin    - 'PulseVectors': vector that contains the pulse itself
%                   must be a numeric vector
%                   default == [] (not used)
%                   - 'SameAsPulse': whether always the same as 
%                                       the current pulse endpoints
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'MeanValueWindowMs': window in ms for 
%                                           calculating mean values
%                   must be a positive scalar
%                   default == 0.5 ms
%                   - 'MinPeakDelayMs': minimum peak delay (ms)
%                               after the end of the pulse
%                   must be a positive scalar
%                   default == 0 ms
%
% Requires:
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/find_pulse_response_endpoints.m
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%       cd/match_dimensions.m
%
% Used by:
%       cd/compute_all_pulse_responses.m
%       cd/compute_average_pulse_response.m
%       cd/find_passive_params.m
%       cd/plot_pulse_response.m

% File History:
% 2018-10-10 Adapted from parse_pulse.m
% 2018-10-11 Fixed tvecRising so that it starts from 0
% 2018-11-13 Added 'MeanValueWindowMs' as an optional argument
% 2018-12-15 Added 'peakValue', 'idxPeak', 'peakAmplitude' in parsedParams
% 2018-12-15 Added 'MinPeakDelayMs' as an optional argument
% 2018-12-17 Now allows siMs to be a vector
% TODO: Compute peakDelaySamples and peakDelayMs
% TODO: Compute maxSlope, halfWidthSamples, halfWidthMs, 
% TODO: Compute timeConstantSamples, timeConstantMs

%% Hard-coded parameters

%% Default values for optional arguments
pulseVectorsDefault = [];       % don't use pulse vectors by default
sameAsPulseDefault = true;      % use pulse endpoints by default
meanValueWindowMsDefault = 0.5; % calculating mean values over 0.5 ms by default
minPeakDelayMsDefault = 0;      % no minimum peak delay by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vectors', ...                   % vectors
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'siMs', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PulseVectors', pulseVectorsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['PulseVectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'SameAsPulse', sameAsPulseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MeanValueWindowMs', meanValueWindowMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'scalar'}));
addParameter(iP, 'MinPeakDelayMs', minPeakDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'scalar'}));

% Read from the Input Parser
parse(iP, vectors, siMs, varargin{:});
pulseVectors = iP.Results.PulseVectors;
sameAsPulse = iP.Results.SameAsPulse;
meanValueWindowMs = iP.Results.MeanValueWindowMs;
minPeakDelayMs = iP.Results.MinPeakDelayMs;

%% Preparation
% Force vectors to be a column cell array
vectors = force_column_cell(vectors);

% Count the number of samples for each vector
nSamples = count_samples(vectors);

% Count the number of vectors
nVectors = count_vectors(vectors);

% Make sure siMs is a column vector
siMs = match_dimensions(siMs, [nVectors, 1]);

% Convert time parameters to samples
meanValueWindowSamples = round(meanValueWindowMs ./ siMs);
minPeakDelaySamples = round(minPeakDelayMs ./ siMs);

% For clarity, define a column vector of ones
allOnes = ones(nVectors, 1);

%% Do the job
% Find indices for all the pulse response endpoints
[idxResponseStart, idxResponseEnd, hasJump, idxPulseStart, idxPulseEnd] = ...
    find_pulse_response_endpoints(vectors, siMs, ...
                                    'PulseVectors', pulseVectors, ...
                                    'SameAsPulse', sameAsPulse, ...
                                    'ResponseLengthMs', 0, ...
                                    'BaselineLengthMs', 0);

% Find indices for all the pulse response midpoints
idxResponseMid = round((idxResponseStart + idxResponseEnd) ./ 2);

% Find all pulse widths in samples
responseWidthSamples = idxResponseEnd - idxResponseStart;

% Convert pulse widths to milliseconds
responseWidthMs = responseWidthSamples .* siMs;

% Compute the endpoints of the baseline
%   Note: this cannot be smaller than 1
idxBaseStart = max([allOnes, idxResponseStart - meanValueWindowSamples], [], 2);
idxBaseEnd = idxResponseStart - 1;

% Compute the endpoints of the steady state
%   Note: this cannot be smaller than idxResponseStart
% TODO: Use argfun.m
allRespStart = idxResponseStart .* allOnes;
idxSteadyStartIdeal = idxResponseEnd - meanValueWindowSamples;
idxSteadyEndIdeal = idxResponseEnd - 1;
idxSteadyStart = max([allRespStart, idxSteadyStartIdeal], [], 2);
idxSteadyEnd = max([allRespStart, idxSteadyStart], [], 2);

% Construct a vector that goes from -n:-1
indBefore = (-1) * transpose(fliplr(1:meanValueWindowSamples));

% Find the baseline indices using a window before response start
indBase = arrayfun(@(x) x + indBefore, ...
                    idxResponseStart, 'UniformOutput', false);

% Find the steady state indices using a window before response end
indSteady = arrayfun(@(x) x + indBefore, ...
                    idxResponseEnd, 'UniformOutput', false);

% Make sure nothing is out of bounds
indBase = cellfun(@(x) x(x >= 1), indBase, 'UniformOutput', false);
indSteady = cellfun(@(x, y) x(x <= y), indSteady, num2cell(nSamples), ...
                    'UniformOutput', false);

% TODO: Use argfun.m
% Find the average baseline value
baseValue = cellfun(@(x, y) mean(x(y)), vectors, indBase);

% Find the average steady state value
steadyValue = cellfun(@(x, y) mean(x(y)), vectors, indSteady);

% Find the steady state amplitudes
steadyAmplitude = steadyValue - baseValue;

% Find the minimum and maximum values
% TODO: Use extract_elements(vectors, 'min', 'MaxNum', 1)
% TODO: Use extract_elements(vectors, 'max', 'MaxNum', 1)
minValue = cellfun(@min, vectors);
maxValue = cellfun(@max, vectors);

% Find the index to begin searching for the peak
%   Note: this cannot be greater than nSamples
idxPeakSearchBegin = arrayfun(@(x, y, z) min([x + y, z]), ...
                                idxPulseEnd, minPeakDelaySamples, nSamples);

% Find the minimum and maximum values after the peak search begins
% TODO: Use argfun.m
[minValueAfterPulse, idxMinValueAfterPulseRel] = ...
    cellfun(@(x, y) min(x(y:end)), vectors, num2cell(idxPeakSearchBegin));
[maxValueAfterPulse, idxMaxValueAfterPulseRel] = ...
    cellfun(@(x, y) max(x(y:end)), vectors, num2cell(idxPeakSearchBegin));

% Record the corresponding indices
idxMinValueAfterPulse = idxMinValueAfterPulseRel + (idxPeakSearchBegin - 1);
idxMaxValueAfterPulse = idxMaxValueAfterPulseRel + (idxPeakSearchBegin - 1);

% Find the peak values (the minimum or maximum value after pulse end 
%   + minPeakDelaySamples with largest magnitude)
[peakValue, idxPeak] = ...
    arrayfun(@(x, y, z, w) choose_peak_value(x, y, z, w), ...
            minValueAfterPulse, maxValueAfterPulse, ...
            idxMinValueAfterPulse, idxMaxValueAfterPulse);

% Compute the relative peak amplitude
peakAmplitude = peakValue - baseValue;

% Find the indices for the rising and falling phases, respectively
% TODO: Use argfun and make a function create_indices.m
%   indices = create_indices(idxStart, idxEnd)
indRising = arrayfun(@(x, y) transpose(x:y), ...
                    idxResponseStart, idxResponseEnd, ...
                    'UniformOutput', false);
indFalling = arrayfun(@(x, y) transpose(x:y), ...
                    idxResponseEnd, nSamples, ...
                    'UniformOutput', false);
indCombined = arrayfun(@(x, y) transpose(x:y), ...
                    idxResponseStart, nSamples, ...
                    'UniformOutput', false);

% Count the number of samples in the rising and falling phases, respectively
[nSamplesRising, nSamplesFalling, nSamplesCombined] = ...
    argfun(@count_samples, indRising, indFalling, indCombined);

% Convert base values to a cell array
baseValueCell = num2cell(baseValue);

% Generate shifted rising/falling phase vectors so that time starts at zero
%   and steady state value is zero
% Note: This will make curve fitting easier
% TODO: Use argfun
tvecsRising = arrayfun(@(x, y) transpose((1:x) - 1) * y, ...
                        nSamplesRising, siMs, 'UniformOutput', false);
vvecsRising = cellfun(@(x, y, z) x(y) - z, ...
                        vectors, indRising, baseValueCell, ...
                        'UniformOutput', false);
tvecsFalling = arrayfun(@(x, y) transpose((1:x) - 1) * y, ...
                        nSamplesRising, siMs, 'UniformOutput', false);
vvecsFalling = cellfun(@(x, y, z) x(y) - z, ...
                        vectors, indFalling, baseValueCell, ...
                        'UniformOutput', false);

% Generate shifted pulse response vectors so that time starts at zero
%   and steady state value is zero
tvecsCombined = arrayfun(@(x, y) transpose((1:x) - 1) * y, ...
                        nSamplesCombined, siMs, 'UniformOutput', false);
vvecsCombined = cellfun(@(x, y, z) x(y) - z, ...
                        vectors, indCombined, baseValueCell, ...
                        'UniformOutput', false);

%% Store results in output
parsedParams = table(nSamples, siMs, responseWidthSamples, responseWidthMs, ...
                        nSamplesRising, nSamplesFalling, nSamplesCombined, ...
                        idxResponseStart, idxResponseEnd, idxResponseMid, ...
                        idxPulseStart, idxPulseEnd, ...
                        idxBaseStart, idxBaseEnd, ...
                        idxSteadyStart, idxSteadyEnd, ...
                        baseValue, steadyValue, steadyAmplitude, ...
                        minValue, maxValue, ...
                        minValueAfterPulse, idxMinValueAfterPulse, ...
                        maxValueAfterPulse, idxMaxValueAfterPulse, ...
                        idxPeakSearchBegin, idxPeak, ...
                        peakValue, peakAmplitude, hasJump);
parsedData = table(vectors, indBase, indSteady, ...
                    indRising, indFalling, indCombined, ...
                    tvecsRising, vvecsRising, tvecsFalling, vvecsFalling, ...
                    tvecsCombined, vvecsCombined);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [peakValue, idxPeak] = ...
                choose_peak_value (minValueAfterPulse, maxValueAfterPulse, ...
                                idxMinValueAfterPulse, idxMaxValueAfterPulse)
%% Chooses the peak value from minimum and maximum

% Find the larger magnitude of the two extrema
[~, iPeak] = max(abs([minValueAfterPulse, maxValueAfterPulse]));

% Get the actual index and value for the peak
if iPeak == 1
    idxPeak = idxMinValueAfterPulse;
    peakValue = minValueAfterPulse;
elseif iPeak == 2
    idxPeak = idxMaxValueAfterPulse;
    peakValue = maxValueAfterPulse;
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

tvecsRising = arrayfun(@(x) transpose(1:x) * siMs, nSamplesRising, ...
                        'UniformOutput', false);
tvecsFalling = arrayfun(@(x) transpose(1:x) * siMs, nSamplesFalling, ...
                        'UniformOutput', false);

idxBaseStart = idxResponseStart - meanValueWindowSamples;
idxSteadyStart = idxResponseEnd - meanValueWindowSamples;

% Find the index for each vector after pulse ends
idxAfterPulseEnd = arrayfun(@(x, y) min([x + 1, y]), idxPulseEnd, nSamples);

nSamplesRising = count_samples(indRising);
nSamplesFalling = count_samples(indFalling);
nSamplesCombined = count_samples(indCombined);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%