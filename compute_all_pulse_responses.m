function [tVecAll, respAll, stimAll, featuresAll] = ...
                compute_all_pulse_responses (fileName, responseType, varargin)
%% Filter and extract all pulse response and compute features
% Usage: [tVecAll, respAll, stimAll, featuresAll] = ...
%               compute_all_pulse_responses (fileName, responseType, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       tVecAll     - time vector for each response
%                   specified as a numeric column vector
%       respAll     - pulse response vector(s)
%                   specified as a numeric column vector
%                       or a cell array of column vectors
%       stimAll     - stimulation pulse vector(s)
%                   specified as a numeric column vector
%                       or a cell array of column vectors
%       featuresAll - computed features of the average pulse response, 
%                       a table with variables:
%                           those returned by parse_pulse_response.m, 
%                           labels: labels for respAll and stimAll, respectively
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
%       cd/compute_average_trace.m
%       cd/extract_channel.m
%       cd/extract_subvectors.m
%       cd/find_pulse_response_endpoints.m
%       cd/force_column_cell.m
%       cd/freqfilter.m
%       cd/parse_abf.m
%       cd/parse_pulse_response.m
%
% Used by:
%       cd/compute_and_plot_all_responses.m

% File History:
% 2018-12-15 Modified from compute_average_pulse_response.m
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

%% Filter and extract pulse response(s)
[tVecsResponse, respVecsResponse, stimVecsResponse, labels] = ...
    filter_and_extract_pulse_response(fileName, responseType, varargin{:});

%% Average the pulse response
% Average the stimulation pulses and pulse responses
[respAll, stimAll] = ...
    argfun(@(x) compute_average_trace(x, 'AlignMethod', 'LeftAdjust'), ...
            respVecsResponse, stimVecsResponse);

%% Compute features of the average pulse response
% Compute the sampling interval
siMs = tVecAll(2) - tVecAll(1);

% Parse the average pulse response vector
featuresAll = parse_pulse_response(respAll, siMs, ...
                                'PulseVectors', stimAll, ...
                                'SameAsPulse', true, ...
                                'MeanValueWindowMs', baselineLengthMs);

% Add labels to features
featuresAll.labels = labels;

% Use the file name as the row name
featuresAll.Properties.RowNames = {fileName};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%