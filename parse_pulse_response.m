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
%                           hasJump
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
%
% Requires:
%       cd/count_samples.m
%       cd/find_pulse_response_endpoints.m
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%
% Used by:
%       cd/compute_average_pulse_response.m
%       cd/find_passive_params.m
%       cd/plot_pulse_response.m

% File History:
% 2018-10-10 Adapted from parse_pulse.m
% 2018-10-11 Fixed tvecRising so that it starts from 0
% 2018-11-13 Added 'MeanValueWindowMs' as an optional argument
% 2018-12-15 Added 'peakValue', 'idxPeak', 'peakAmplitude' in parsedParams
% TODO: 2018-11-28 Now allows siMs to be a vector
% TODO: Compute the slope of the peak

%% Hard-coded parameters

%% Default values for optional arguments
pulseVectorsDefault = [];       % don't use pulse vectors by default
sameAsPulseDefault = true;      % use pulse endpoints by default
meanValueWindowMsDefault = 0.5; % calculating mean values over 0.5 ms by default

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

% Read from the Input Parser
parse(iP, vectors, siMs, varargin{:});
pulseVectors = iP.Results.PulseVectors;
sameAsPulse = iP.Results.SameAsPulse;
meanValueWindowMs = iP.Results.MeanValueWindowMs;

%% Preparation
% Force vectors to be a column cell array
vectors = force_column_cell(vectors);

% Count the number of samples for each vector
nSamples = count_samples(vectors);

% Calculate the mean value calculation window in samples
meanValueWindowSamples = round(meanValueWindowMs ./ siMs);

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
responseWidthMs = responseWidthSamples * siMs;

% Compute the endpoints of the baseline
%   Note: this cannot be smaller than 1
idxBaseStart = max([1, idxResponseStart - meanValueWindowSamples]);
idxBaseEnd = idxResponseStart - 1;

% Compute the endpoints of the steady state
%   Note: this cannot be smaller than idxResponseStart
idxSteadyStart = max([idxResponseStart, idxResponseEnd - meanValueWindowSamples]);
idxSteadyEnd = max([idxResponseStart, idxResponseEnd - 1]);

% Construct a vector that goes from -n:-1
indBefore = (-1) * transpose(fliplr(1:meanValueWindowSamples));

% Find the baseline indices using a window before response start
indBase = arrayfun(@(x) x + indBefore, ...
                    idxResponseStart, 'UniformOutput', false);

% Find the average baseline value
baseValue = cellfun(@(x, y) mean(x(y)), vectors, indBase);

% Find the steady state indices using a window before response end
indSteady = arrayfun(@(x) x + indBefore, ...
                    idxResponseEnd, 'UniformOutput', false);

% Find the average steady state value
steadyValue = cellfun(@(x, y) mean(x(y)), vectors, indSteady);

% Find the steady state amplitudes
steadyAmplitude = steadyValue - baseValue;

% Find the minimum and maximum values
minValue = cellfun(@min, vectors);
maxValue = cellfun(@max, vectors);

% Find the index for each vector after pulse ends
%   Note: this cannot be greater than nSamples
idxAfterPulseEnd = arrayfun(@(x) min([idxPulseEnd + 1, x]), nSamples);

% Find the minimum and maximum values after the pulse ends
[minValueAfterPulse, idxMinValueAfterPulseRel] = ...
    cellfun(@(x, y) min(x(y:end)), vectors, num2cell(idxAfterPulseEnd));
[maxValueAfterPulse, idxMaxValueAfterPulseRel] = ...
    cellfun(@(x, y) max(x(y:end)), vectors, num2cell(idxAfterPulseEnd));

% Record the corresponding indices
idxMinValueAfterPulse = idxMinValueAfterPulseRel + (idxAfterPulseEnd - 1);
idxMaxValueAfterPulse = idxMaxValueAfterPulseRel + (idxAfterPulseEnd - 1);

% Find the peak values (the minimum or maximum value after pulse end
%   with largest magnitude)
[peakValue, idxPeak] = ...
    arrayfun(@(x, y, z, w) choose_peak_value(x, y, z, w), ...
            minValueAfterPulse, maxValueAfterPulse, ...
            idxMinValueAfterPulse, idxMaxValueAfterPulse);

% Compute the relative peak amplitude
peakAmplitude = peakValue - baseValue;

% Find the indices for the rising and falling phases, respectively
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
nSamplesRising = count_samples(indRising);
nSamplesFalling = count_samples(indFalling);
nSamplesCombined = count_samples(indCombined);

% Convert base values to a cell array
baseValueCell = num2cell(baseValue);

% Generate shifted rising/falling phase vectors so that time starts at zero
%   and steady state value is zero
% Note: This will make curve fitting easier
tvecsRising = arrayfun(@(x) transpose((1:x) - 1) * siMs, nSamplesRising, ...
                        'UniformOutput', false);
vvecsRising = cellfun(@(x, y, z) x(y) - z, ...
                        vectors, indRising, baseValueCell, ...
                        'UniformOutput', false);
tvecsFalling = arrayfun(@(x) transpose((1:x) - 1) * siMs, nSamplesFalling, ...
                        'UniformOutput', false);
vvecsFalling = cellfun(@(x, y, z) x(y) - z, ...
                        vectors, indFalling, baseValueCell, ...
                        'UniformOutput', false);

% Generate shifted pulse response vectors so that time starts at zero
%   and steady state value is zero
tvecsCombined = arrayfun(@(x) transpose((1:x) - 1) * siMs, nSamplesCombined, ...
                        'UniformOutput', false);
vvecsCombined = cellfun(@(x, y, z) x(y) - z, ...
                        vectors, indCombined, baseValueCell, ...
                        'UniformOutput', false);

%% Store results in output
parsedParams = table(nSamples, responseWidthSamples, responseWidthMs, ...
                        nSamplesRising, nSamplesFalling, ...
                        idxResponseStart, idxResponseEnd, idxResponseMid, ...
                        idxPulseStart, idxPulseEnd, ...
                        idxBaseStart, idxBaseEnd, ...
                        idxSteadyStart, idxSteadyEnd, ...
                        baseValue, steadyValue, steadyAmplitude, ...
                        minValue, maxValue, ...
                        minValueAfterPulse, idxMinValueAfterPulse, ...
                        maxValueAfterPulse, idxMaxValueAfterPulse, ...
                        peakValue, idxPeak, peakAmplitude, hasJump);
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


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%