function [tVecResp, vecResp, vecStim, labels] = ...
                compute_average_pulse_response (fileName, responseType, varargin)
%% Computes an average pulse response
% Usage: [tVecResp, vecResp, vecStim, labels] = ...
%               compute_average_pulse_response (fileName, responseType, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       tVecResp    - time vector for the average response
%                   specified as a numeric column vector
%       vecResp     - vector for the average response
%                   specified as a numeric column vector
%       vecStim     - vector for the average stimulation pulse
%                   specified as a numeric column vector
%       labels      - labels for vecResp and vecStim, respectively
%                   specified as a 2-element cell array of character vectors
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
%       varargin    - 'ChannelTypes': the channel types
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
%       cd/extract_channel.m
%       cd/extract_subvectors.m
%       cd/find_pulse_response_endpoints.m
%       cd/force_column_cell.m
%       cd/parse_abf.m
%
% Used by:
%       cd/compute_and_plot_evoked_LFP.m

% File History:
% 2018-12-15 Moved from compute_and_plot_evoked_LFP.m
% TODO: Make baselineLengthMs and responseLengthMs optional parameters
% 

%% Hard-coded parameters
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Other'};
baselineLengthMs = 5;           % baseline length in ms
responseLengthMs = 20;          % response length in ms

%% Default values for optional arguments
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
addParameter(iP, 'ChannelTypes', channelTypesDefault, ...
    @(x) isempty(x) || isstring(x) || iscellstr(x));
addParameter(iP, 'ChannelUnits', channelUnitsDefault, ...
    @(x) isempty(x) || isstring(x) || iscellstr(x));
addParameter(iP, 'ChannelLabels', channelLabelsDefault, ...
    @(x) isempty(x) || isstring(x) || iscellstr(x));
addParameter(iP, 'ParsedParams', parsedParamsDefault, ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));
addParameter(iP, 'ParsedData', parsedDataDefault, ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));

% Read from the Input Parser
parse(iP, fileName, responseType, varargin{:});
channelTypes = iP.Results.ChannelTypes;
channelUnits = iP.Results.ChannelUnits;
channelLabels = iP.Results.ChannelLabels;
parsedParams = iP.Results.ParsedParams;
parsedData = iP.Results.ParsedData;

% Validate the channel type
responseType = validatestring(responseType, validChannelTypes);

%% Preparation
% Decide on the stimulation type based on the response type
switch responseType
    case 'Voltage'
        stimType = 'Current';
    case {'Current', 'Conductance', 'Other'}
        stimType = 'Voltage';
    otherwise
        error('Code logic error!!');
end

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
channelLabels = parsedParams.channelLabels;

% Extract the time vector
tVec = parsedData.tVec;

% Initialize labels
labels = cell(1, 2);

% Extract the vectors containing pulses
[stimVecs, labels{1}] = ...
    extract_channel(fileName, stimType, ...
                'ParsedParams', parsedParams, 'ParsedData', parsedData, ...
                'ChannelTypes', channelTypes);

% Extract the vectors containing pulse responses
[respVecs, labels{2}] = ...
    extract_channel(fileName, responseType, ...
                'ParsedParams', parsedParams, 'ParsedData', parsedData, ...
                'ChannelTypes', channelTypes);

%% Do the job
% Force as a cell array of vectors
[respVecs, stimVecs] = argfun(@force_column_cell, respVecs, stimVecs);

% Compute the baseline length in samples
baselineLengthSamples = floor(baselineLengthMs / siMs);

% Identify the pulse response endpoints
[idxResponseStarts, idxResponseEnds, ~, ~, idxStimEnds] = ...
    find_pulse_response_endpoints(respVecs, siMs, ...
                                'PulseVectors', stimVecs, ...
                                'BaselineLengthMs', baselineLengthMs, ...
                                'ResponseLengthMs', responseLengthMs);

% Compute the number of samples in each pulse response
nSamplesEachResponse = idxResponseEnds - idxResponseStarts + 1;

% Determine number of samples in the shortest pulse response
%   This will be the length of the average response
nSamplesAvgResponse = min(nSamplesEachResponse);

% Extract the time vector using the average starting index
idxStartAveraged = max(1, round(mean(idxResponseStarts)));
idxEndAveraged = (idxStartAveraged - 1) + nSamplesAvgResponse;
tVecResp = tVec(idxStartAveraged:idxEndAveraged);

% Place endpoints together as a matrix, with each column corresponding
%   to each vector
endPointsResponse = transpose([idxResponseStarts, idxResponseEnds]);

% Extract the pulse responses
[respVecsResponse, stimVecsResponse] = ...
    argfun(@(x) extract_subvectors(x, 'Endpoints', endPointsResponse), ...
            respVecs, stimVecs);

% Average the pulse responses to get the evoked local field potential
%   and average the current pulses to get the stimulation pulse
[vecResp, vecStim] = ...
    argfun(@(x) compute_average_trace(x, 'AlignMethod', 'LeftAdjust'), ...
            respVecsResponse, stimVecsResponse);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Extract the vectors containing pulses responses and pulses, respectively
[respVecs, stimVecs] = ...
    argfun(@(x) extract_channel(fileName, x, ...
                'ParsedParams', parsedParams, 'ParsedData', parsedData, ...
                'ChannelTypes', channelTypes), stimType, responseType);


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%