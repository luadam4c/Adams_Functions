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
%                           siMs
%                           meanValueWindowMs
%                           minPeakDelayMs
%                           meanValueWindowSamples
%                           minPeakDelaySamples
%                           nSamples
%                           responseWidthSamples
%                           responseWidthMs
%                           nSamplesRising
%                           nSamplesFalling
%                           idxResponseStart
%                           idxResponseEnd
%                           idxResponseMid
%                           idxBeforePulseStart
%                           idxBeforePulseEnd
%                           idxBaseStart
%                           idxBaseEnd
%                           idxSteadyStart
%                           idxSteadyEnd
%                           baseValue
%                           steadyValue
%                           steadyAmplitude
%                           minValue
%                           maxValue
%                           minValueAfterMinDelay
%                           idxMinValueAfterMinDelay
%                           maxValueAfterMinDelay
%                           idxMaxValueAfterMinDelay
%                           idxMinPeakTime
%                           idxPeak
%                           peakValue
%                           halfPeakValue
%                           peakAmplitude
%                           peakDelaySamples
%                           peakDelayMs
%                           hasJump
%                           nSamplesPulse
%                           idxAfterPulseStart
%                           idxAfterPulseEnd
%                           idxPulseMidpoint
%                           pulseBaseValue
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
%                           pulseVectors
%                           indPulseBase
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
%                   - 'ResponseLengthMs': length of the pulse response
%                                           after pulse endpoint in ms
%                   must be a nonnegative scalar
%                   default = 20 ms
%                   - 'BaselineLengthMs': length of the pulse response
%                                           before pulse start in ms
%                   must be a nonnegative scalar
%                   default = 0 ms
%
% Requires:
%       cd/argfun.m
%       cd/compute_means.m
%       cd/convert_to_samples.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/create_error_for_nargin.m
%       cd/create_indices.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/find_pulse_response_endpoints.m
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%       cd/match_dimensions.m
%       cd/parse_pulse.m
%       cd/renamevars.m
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
% 2018-12-17 Now computes peakDelaySamples and peakDelayMs
% 2018-12-17 Added all detected results of parse_pulse.m
% 2018-12-18 Added halfPeakValue
% TODO: Compute maxSlope, halfWidthSamples, halfWidthMs, 
% TODO: Compute timeConstantSamples, timeConstantMs

%% Hard-coded parameters
% Pulse parameter names to change or remove
%   Note: must be consistent with parse_pulse.m
pulseParamsToRemove = {'idxBeforeStart', 'idxBeforeEnd'};
pulseParamNamesOld = {'nSamples', 'idxAfterStart', ...
                        'idxAfterEnd', 'idxMidpoint', 'baseValue'};
pulseParamNamesNew = {'nSamplesPulse', 'idxAfterPulseStart', ...
                    'idxAfterPulseEnd', 'idxPulseMidpoint', 'pulseBaseValue'};
pulseDataNamesOld = {'vectors', 'indBase'};
pulseDataNamesNew = {'pulseVectors', 'indPulseBase'};

%% Default values for optional arguments
pulseVectorsDefault = [];       % don't use pulse vectors by default
sameAsPulseDefault = true;      % use pulse endpoints by default
meanValueWindowMsDefault = 0.5; % calculating mean values over 0.5 ms by default
minPeakDelayMsDefault = 0;      % no minimum peak delay by default
responseLengthMsDefault = 0;    % 0 ms after pulse end by default
baselineLengthMsDefault = 0;    % don't include a baseline by default

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
addParameter(iP, 'ResponseLengthMs', responseLengthMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'BaselineLengthMs', baselineLengthMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));

% Read from the Input Parser
parse(iP, vectors, siMs, varargin{:});
pulseVectors = iP.Results.PulseVectors;
sameAsPulse = iP.Results.SameAsPulse;
meanValueWindowMs = iP.Results.MeanValueWindowMs;
minPeakDelayMs = iP.Results.MinPeakDelayMs;
responseLengthMs = iP.Results.ResponseLengthMs;
baselineLengthMs = iP.Results.BaselineLengthMs;

%% Preparation
% Force vectors to be a column cell array
vectors = force_column_cell(vectors);

% Count the number of samples for each vector
nSamples = count_samples(vectors);

% Count the number of vectors
nVectors = count_vectors(vectors);

% Convert time parameters to samples
[meanValueWindowSamples, minPeakDelaySamples] = ...
    argfun(@(x) convert_to_samples(x, siMs), meanValueWindowMs, minPeakDelayMs);

% Make sure parameters are column vectors
[siMs, meanValueWindowMs, minPeakDelayMs, ...
    meanValueWindowSamples, minPeakDelaySamples] = ...
    argfun(@(x) match_dimensions(x, [nVectors, 1]), ...
            siMs, meanValueWindowMs, minPeakDelayMs, ...
            meanValueWindowSamples, minPeakDelaySamples);

% For clarity, define a column vector of ones
allOnes = ones(nVectors, 1);

%% Do the job
% If not empty, parse the pulse vectors with a given sampling interval
if ~isempty(pulseVectors)
    [pulseParams, pulseData] = ...
        parse_pulse(pulseVectors, 'SamplingIntervalMs', siMs);

    % Remove variables in pulseParams
    pulseParams = removevars(pulseParams, pulseParamsToRemove);

    % Rename variables in pulseParams
    pulseParams = renamevars(pulseParams, pulseParamNamesOld, pulseParamNamesNew);

    % Rename variables in pulseData
    pulseData = renamevars(pulseData, pulseDataNamesOld, pulseDataNamesNew);
end

% Find indices for all the pulse response endpoints
[idxResponseStart, idxResponseEnd, hasJump, ...
    idxBeforePulseStart, idxBeforePulseEnd] = ...
    find_pulse_response_endpoints(vectors, siMs, ...
                'PulseVectors', pulseVectors, 'SameAsPulse', sameAsPulse, ...
                'ResponseLengthMs', responseLengthMs, ...
                'BaselineLengthMs', baselineLengthMs);

% Find indices for all the pulse response midpoints
idxResponseMid = round((idxResponseStart + idxResponseEnd) ./ 2);

% Find all response widths in samples
responseWidthSamples = idxResponseEnd - idxResponseStart;

% Convert response widths to milliseconds
responseWidthMs = responseWidthSamples .* siMs;

% Compute the endpoints of the baseline
%   Note: this cannot be smaller than 1
[idxBaseStart, idxBaseEnd] = ...
    argfun(@(x) max([allOnes, idxResponseStart - x], [], 2), ...
            meanValueWindowSamples, 1);

% Compute the endpoints of the steady state
%   Note: this cannot be smaller than idxResponseStart
allRespStart = idxResponseStart .* allOnes;
[idxSteadyStart, idxSteadyEnd] = ...
    argfun(@(x) max([allRespStart, idxResponseEnd - x], [], 2), ...
            meanValueWindowSamples, 1);

% Construct a vector that goes from -n:-1
indBefore = transpose(-meanValueWindowSamples:-1);

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

% Find the average baseline and steady state value
[baseValue, steadyValue] = ...
    argfun(@(x) compute_means(vectors, 'Indices', x), indBase, indSteady);

% Find the steady state amplitudes
steadyAmplitude = steadyValue - baseValue;

% Find the minimum and maximum values
[minValue, maxValue] = argfun(@(x) extract_elements(vectors, x), 'min', 'max');

% Find the index to begin searching for the peak (the minimum peak time)
%   Note: this cannot be greater than nSamples
idxMinPeakTime = arrayfun(@(x, y, z) min([x + y, z]), ...
                            idxBeforePulseEnd, minPeakDelaySamples, nSamples);

% Find the endpoints for the peak search
endPointsPeakSearch = transpose([idxMinPeakTime, nSamples]);

% Extract the parts of the vector after the minimum peak time(s)
vecsPeakSearch = extract_subvectors(vectors, 'Endpoints', endPointsPeakSearch);

% Find the minimum and maximum values after the minimum peak time(s)\
[minValueAfterMinDelay, idxMinValueAfterMinDelayRel] = ...
    extract_elements(vecsPeakSearch, 'min');
[maxValueAfterMinDelay, idxMaxValueAfterMinDelayRel] = ...
    extract_elements(vecsPeakSearch, 'max');

% Record the corresponding indices
[idxMinValueAfterMinDelay, idxMaxValueAfterMinDelay] = ...
    argfun(@(x) x + idxMinPeakTime - 1, ...
            idxMinValueAfterMinDelayRel, idxMaxValueAfterMinDelayRel);

% Find the peak values (the minimum or maximum value after pulse end 
%   + minPeakDelaySamples with largest magnitude)
[peakValue, idxPeak] = ...
    arrayfun(@(x, y, z, w) choose_peak_value(x, y, z, w), ...
            minValueAfterMinDelay, maxValueAfterMinDelay, ...
            idxMinValueAfterMinDelay, idxMaxValueAfterMinDelay);

% Compute the half peak value
halfPeakValue = mean([peakValue, baseValue], 2);

% Compute the relative peak amplitude
peakAmplitude = peakValue - baseValue;

% Compute the peak delay (pulse end to peak) in samples
peakDelaySamples = idxPeak - idxBeforePulseEnd;

% Compute the peak delay (pulse end to peak) in ms
% Use convert_to_time.m
peakDelayMs = peakDelaySamples .* siMs;

% Find the endpoint for the different phases
endPointsRising = transpose([idxResponseStart, idxResponseEnd]);
endPointsFalling = transpose([idxResponseEnd, nSamples]);
endPointsCombined = transpose([idxResponseStart, nSamples]);

% Construct the corresponding indices
% TODO: Use argfun and make a function create_indices.m
[indRising, indFalling, indCombined] = ...
    argfun(@(x) create_indices(x, 'ForceCellOutput', true), ...
            endPointsRising, endPointsFalling, endPointsCombined);

% Count the number of samples in the rising and falling phases, respectively
[nSamplesRising, nSamplesFalling, nSamplesCombined] = ...
    argfun(@count_samples, indRising, indFalling, indCombined);

% Convert base values to a cell array
baseValueCell = num2cell(baseValue);

% Generate shifted rising/falling phase vectors so that time starts at zero
%   and steady state value is zero
% Note: This will make curve fitting easier
% TODO: Use argfun.m and create_time_vectors.m
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
% Put together the pulse response parameters
parsedParams = table(siMs, meanValueWindowMs, minPeakDelayMs, ...
                        meanValueWindowSamples, minPeakDelaySamples, ...
                        nSamples, responseWidthSamples, responseWidthMs, ...
                        nSamplesRising, nSamplesFalling, nSamplesCombined, ...
                        idxResponseStart, idxResponseEnd, idxResponseMid, ...
                        idxBeforePulseStart, idxBeforePulseEnd, ...
                        idxBaseStart, idxBaseEnd, ...
                        idxSteadyStart, idxSteadyEnd, ...
                        baseValue, steadyValue, steadyAmplitude, ...
                        minValue, maxValue, ...
                        minValueAfterMinDelay, idxMinValueAfterMinDelay, ...
                        maxValueAfterMinDelay, idxMaxValueAfterMinDelay, ...
                        idxMinPeakTime, idxPeak, ...
                        peakValue, halfPeakValue, peakAmplitude, ...
                        peakDelaySamples, peakDelayMs, hasJump);

% Put together the pulse response data
parsedData = table(vectors, indBase, indSteady, ...
                    indRising, indFalling, indCombined, ...
                    vecsPeakSearch, ...
                    tvecsRising, vvecsRising, tvecsFalling, vvecsFalling, ...
                    tvecsCombined, vvecsCombined);

% Append pulse params and data to the table if provided
if ~isempty(pulseVectors)
    parsedParams = horzcat(parsedParams, pulseParams);
    parsedData = horzcat(parsedData, pulseData);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [peakValue, idxPeak] = ...
                choose_peak_value (minValueAfterMinDelay, maxValueAfterMinDelay, ...
                                idxMinValueAfterMinDelay, idxMaxValueAfterMinDelay)
%% Chooses the peak value from minimum and maximum

% Find the larger magnitude of the two extrema
[~, iPeak] = max(abs([minValueAfterMinDelay, maxValueAfterMinDelay]));

% Get the actual index and value for the peak
if iPeak == 1
    idxPeak = idxMinValueAfterMinDelay;
    peakValue = minValueAfterMinDelay;
elseif iPeak == 2
    idxPeak = idxMaxValueAfterMinDelay;
    peakValue = maxValueAfterMinDelay;
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
idxAfterMinDelayEnd = arrayfun(@(x, y) min([x + 1, y]), idxBeforePulseEnd, nSamples);

nSamplesRising = count_samples(indRising);
nSamplesFalling = count_samples(indFalling);
nSamplesCombined = count_samples(indCombined);

siMs = match_dimensions(siMs, [nVectors, 1]);

meanValueWindowSamples = round(meanValueWindowMs ./ siMs);
minPeakDelaySamples = round(minPeakDelayMs ./ siMs);

indBefore = (-1) * transpose(fliplr(1:meanValueWindowSamples));

minValue = cellfun(@min, vectors);
maxValue = cellfun(@max, vectors);

idxBaseStart = max([allOnes, idxResponseStart - meanValueWindowSamples], [], 2);
idxBaseEnd = max([allOnes, idxResponseStart - 1], [], 2);
idxSteadyStartIdeal = idxResponseEnd - meanValueWindowSamples;
idxSteadyEndIdeal = idxResponseEnd - 1;
idxSteadyStart = max([allRespStart, idxSteadyStartIdeal], [], 2);
idxSteadyEnd = max([allRespStart, idxSteadyEndIdeal], [], 2);

idxMinValueAfterMinDelay = idxMinValueAfterMinDelayRel + (idxMinPeakTime - 1);
idxMaxValueAfterMinDelay = idxMaxValueAfterMinDelayRel + (idxMinPeakTime - 1);

[minValueAfterMinDelay, idxMinValueAfterMinDelayRel] = ...
    cellfun(@(x, y) min(x(y:end)), vectors, num2cell(idxMinPeakTime));
[maxValueAfterMinDelay, idxMaxValueAfterMinDelayRel] = ...
    cellfun(@(x, y) max(x(y:end)), vectors, num2cell(idxMinPeakTime));

% Find the average baseline value
baseValue = cellfun(@(x, y) mean(x(y)), vectors, indBase);
% Find the average steady state value
steadyValue = cellfun(@(x, y) mean(x(y)), vectors, indSteady);

indRising = arrayfun(@(x, y) transpose(x:y), ...
                    idxResponseStart, idxResponseEnd, ...
                    'UniformOutput', false);
indFalling = arrayfun(@(x, y) transpose(x:y), ...
                    idxResponseEnd, nSamples, ...
                    'UniformOutput', false);
indCombined = arrayfun(@(x, y) transpose(x:y), ...
                    idxResponseStart, nSamples, ...
                    'UniformOutput', false);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%