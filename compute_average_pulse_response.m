function [tVecAvg, respAvg, stimAvg, featuresAvg] = ...
                compute_average_pulse_response (fileName, responseType, varargin)
%% Computes an average pulse response as well as its features
% Usage: [tVecAvg, respAvg, stimAvg, featuresAvg] = ...
%               compute_average_pulse_response (fileName, responseType, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       tVecAvg     - time vector for average response
%                   specified as a numeric column vector
%       respAvg     - average pulse response vector
%                   specified as a numeric column vector
%       stimAvg     - average stimulation pulse vector
%                   specified as a numeric column vector
%       featuresAvg - computed features of the average pulse response, 
%                       a table with variables:
%                           those returned by parse_pulse_response.m, 
%                           labels: labels for respAvg and stimAvg, respectively
%                   specified as a table of 1 row
% Arguments:
%       fileName    - file name could be either the full path or 
%                       a relative path in current directory
%                       .abf is not needed (e.g. 'B20160908_0004')
%                   must be a string scalar or a character vector
%       responseType- the channel type for the pulse response
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'Voltage'       - voltage
%                       'Current'       - current
%                       'Conductance'   - conductance
%                       'Other'         - other un-identified types
%       varargin    - 'LowPassFrequency': frequency of lowpass filter in Hz
%                   must be a nonnegative scalar
%                   default = []
%                   - 'ResponseLengthMs': length of the pulse response
%                                           after pulse endpoint in ms
%                   must be a nonnegative scalar
%                   default = 20 ms
%                   - 'BaselineLengthMs': length of the pulse response
%                                           before pulse start in ms
%                   must be a nonnegative scalar
%                   default = 5 ms
%                   - 'MinPeakDelayMs': minimum peak delay (ms)
%                               after the end of the pulse
%                   must be a positive scalar
%                   default == 0 ms
%                   - 'ChannelTypes': the channel types
%                   must be a cellstr with nChannels elements
%                       each being one of the following:
%                           'Voltage'
%                           'Current'
%                           'Conductance'
%                           'Other'
%                   default == detected with identify_channels()
%                   - 'ParsedParams': parsed parameters returned by parse_abf.m
%                   must be a scalar structure
%                   default == what the file provides
%                   - 'ParsedData': parsed data returned by parse_abf.m
%                   must be a scalar structure
%                   default == what the file provides
%
% Requires:
%       cd/argfun.m
%       cd/compute_average_trace.m
%       cd/create_average_time_vector.m
%       cd/parse_pulse_response.m
%
% Used by:
%       cd/compute_and_plot_average_response.m

% File History:
% 2018-12-15 Moved from compute_and_plot_evoked_LFP.m
% 2018-12-15 Made baselineLengthMs and responseLengthMs optional parameters
% 2018-12-15 Now lowpass filters each trace first
% 

%% Hard-coded parameters
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Other'};

%% Default values for optional arguments
lowPassFrequencyDefault = [];   % do not lowpass filter by default
responseLengthMsDefault = 20;   % a response of 20 ms by default
baselineLengthMsDefault = 5;    % a baseline of 5 ms by default
minPeakDelayMsDefault = 0;      % no minimum peak delay by default
channelTypesDefault = {};       % set later
channelUnitsDefault = {};       % set later
channelLabelsDefault = {};      % set later
parsedParamsDefault = [];       % set later
parsedDataDefault = [];         % set later

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
addRequired(iP, 'fileName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'responseType', ...
    @(x) any(validatestring(x, validChannelTypes)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'LowPassFrequency', lowPassFrequencyDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'ResponseLengthMs', responseLengthMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'BaselineLengthMs', baselineLengthMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'MinPeakDelayMs', minPeakDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'scalar'}));
addParameter(iP, 'ChannelTypes', channelTypesDefault, ...
    @(x) isempty(x) || isstring(x) || iscellstr(x));
addParameter(iP, 'ChannelUnits', channelUnitsDefault, ...
    @(x) isempty(x) || isstring(x) || iscellstr(x));
addParameter(iP, 'ChannelLabels', channelLabelsDefault, ...
    @(x) isempty(x) || isstring(x) || iscellstr(x));
addParameter(iP, 'ParsedParams', parsedParamsDefault, ...
    @(x) isempty(x) || isstruct(x));
addParameter(iP, 'ParsedData', parsedDataDefault, ...
    @(x) isempty(x) || isstruct(x));

% Read from the Input Parser
parse(iP, fileName, responseType, varargin{:});
lowPassFrequency = iP.Results.LowPassFrequency;
responseLengthMs = iP.Results.ResponseLengthMs;
baselineLengthMs = iP.Results.BaselineLengthMs;
minPeakDelayMs = iP.Results.MinPeakDelayMs;
channelTypes = iP.Results.ChannelTypes;
channelUnits = iP.Results.ChannelUnits;
channelLabels = iP.Results.ChannelLabels;
parsedParams = iP.Results.ParsedParams;
parsedData = iP.Results.ParsedData;

%% Filter and extract pulse response(s)
[tVecsResponse, respVecsResponse, stimVecsResponse, labels] = ...
    filter_and_extract_pulse_response(fileName, responseType, varargin{:});

%% Average the pulse response
% Create a new time vector starting from the average start time
[tVecAvg, nSamplesAvg] = create_average_time_vector(tVecsResponse);

% Average the stimulation pulses and pulse responses
[respAvg, stimAvg] = ...
    argfun(@(x) compute_average_trace(x, 'AlignMethod', 'LeftAdjust', ...
                                        'NSamples', nSamplesAvg), ...
            respVecsResponse, stimVecsResponse);


%% Compute features of the average pulse response
% Compute the sampling interval
siMs = tVecAvg(2) - tVecAvg(1);

% Parse the average pulse response vector
featuresAvg = parse_pulse_response(respAvg, siMs, ...
                                'PulseVectors', stimAvg, ...
                                'SameAsPulse', true, ...
                                'MeanValueWindowMs', baselineLengthMs, ...
                                'MinPeakDelayMs', minPeakDelayMs);

% Add labels to features
featuresAvg.labels = labels;

% Use the file name as the row name
featuresAvg.Properties.RowNames = {fileName};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Get the minimum number of samples
nSamples = count_samples(tVecsResponse)

% Find the minimum number of samples
idxStartAveraged = max(1, round(mean(idxResponseStarts)));
idxEndAveraged = (idxStartAveraged - 1) + minNSamples;
tVecResponse = tVec(idxStartAveraged:idxEndAveraged);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%