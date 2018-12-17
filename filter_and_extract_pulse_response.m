function [tVecsResponse, respVecsResponse, stimVecsResponse, labels] = ...
            filter_and_extract_pulse_response (fileName, responseType, varargin)
%% Filters and extracts pulse response(s) from a .abf file
% Usage: [tVecsResponse, respVecsResponse, stimVecsResponse] = ...
%           filter_and_extract_pulse_response (fileName, responseType, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       tVecsResponse       - time vector(s) for average response
%                           specified as cell array of numeric column vectors
%       respVecsResponse    - average pulse response vector(s)
%                           specified as cell array of numeric column vectors
%       stimVecsResponse    - average stimulation pulse vector(s)
%                           specified as cell array of numeric column vectors
%       labels              - labels for response and stimulation, respectively
%                           specified as 2-element cell array of char vectors
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
%       cd/choose_stimulation_type.m
%       cd/extract_channel.m
%       cd/extract_subvectors.m
%       cd/find_pulse_response_endpoints.m
%       cd/force_column_cell.m
%       cd/freqfilter.m
%       cd/parse_abf.m
%
% Used by:
%       cd/compute_all_pulse_responses.m
%       cd/compute_average_pulse_response.m

% File History:
% 2018-12-15 Moved from compute_average_pulse_response.m
% 2018-12-15 Fixed order of labels
% 

%% Hard-coded parameters
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Other'};

%% Default values for optional arguments
lowPassFrequencyDefault = [];   % do not lowpass filter by default
responseLengthMsDefault = 20;   % a response of 20 ms by default
baselineLengthMsDefault = 5;    % a baseline of 5 ms by default
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
iP.KeepUnmatched = true;                        % allow extraneous options

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
channelTypes = iP.Results.ChannelTypes;
channelUnits = iP.Results.ChannelUnits;
channelLabels = iP.Results.ChannelLabels;
parsedParams = iP.Results.ParsedParams;
parsedData = iP.Results.ParsedData;

% Validate the channel type
responseType = validatestring(responseType, validChannelTypes);

%% Preparation
% Load and parse the abf file if parsedParams and parsedData not both provided
if isempty(parsedParams) || isempty(parsedData)
    [parsedParams, parsedData] = ...
        parse_abf(fileName, 'Verbose', false, ...
                  'ChannelTypes', channelTypes, ...
                  'ChannelUnits', channelUnits, ...
                  'ChannelLabels', channelLabels, ...
                  'ExtractChannels', true);
end

% Extract the parsed parameters
siMs = parsedParams.siMs;
siSeconds = parsedParams.siSeconds;
channelTypes = parsedParams.channelTypes;
channelUnits = parsedParams.channelUnits;
channelLabels = parsedParams.channelLabels;

% Extract the time vector
tVec = parsedData.tVec;

% Initialize labels
labels = cell(1, 2);

% Decide on the stimulation type based on the response type
stimType = choose_stimulation_type(responseType);

% Extract the vectors containing pulse responses
[respVecs, labels{1}] = ...
    extract_channel(fileName, responseType, 'MaxNChannels', 1, ...
            'ParsedParams', parsedParams, 'ParsedData', parsedData, ...
            'ChannelTypes', channelTypes, 'ChannelUnits', channelUnits, ...
            'ChannelLabels', channelLabels);

% Extract the vectors containing pulses
[stimVecs, labels{2}] = ...
    extract_channel(fileName, stimType, 'MaxNChannels', 1, ...
            'ParsedParams', parsedParams, 'ParsedData', parsedData, ...
            'ChannelTypes', channelTypes, 'ChannelUnits', channelUnits, ...
            'ChannelLabels', channelLabels);

%% Filter and extract the pulse responses
% Low-pass filter if requested
if ~isempty(lowPassFrequency)
    respVecs = freqfilter(respVecs, lowPassFrequency, ...
                            'FilterType', 'low', 'si', siSeconds);
end

% Force as a cell array of vectors
[tVec, respVecs, stimVecs] = ...
    argfun(@force_column_cell, tVec, respVecs, stimVecs);

% Identify the pulse response endpoints
[idxResponseStarts, idxResponseEnds] = ...
    find_pulse_response_endpoints(respVecs, siMs, ...
                                'PulseVectors', stimVecs, ...
                                'BaselineLengthMs', baselineLengthMs, ...
                                'ResponseLengthMs', responseLengthMs);

% Place endpoints together as a matrix, with each column corresponding
%   to each vector
endPointsResponse = transpose([idxResponseStarts, idxResponseEnds]);

% Extract the pulse responses
[tVecsResponse, respVecsResponse, stimVecsResponse] = ...
    argfun(@(x) extract_subvectors(x, 'Endpoints', endPointsResponse), ...
            tVec, respVecs, stimVecs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Extract the vectors containing pulses responses and pulses, respectively
[respVecs, stimVecs] = ...
    argfun(@(x) extract_channel(fileName, x, ...
                'ParsedParams', parsedParams, 'ParsedData', parsedData, ...
                'ChannelTypes', channelTypes), stimType, responseType);

% Compute the baseline length in samples
baselineLengthSamples = floor(baselineLengthMs / siMs);

% Extract the time vector using the average starting index
idxStartAveraged = max(1, round(mean(idxResponseStarts)));
idxEndAveraged = (idxStartAveraged - 1) + nSamplesAvgResponse;
tVecsResponse = tVec(idxStartAveraged:idxEndAveraged);

% Determine number of samples in the shortest pulse response
%   This will be the length of the average response
nSamplesAvgResponse = min(nSamplesEachResponse);

% Compute the number of samples in each pulse response
nSamplesEachResponse = idxResponseEnds - idxResponseStarts + 1;


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%