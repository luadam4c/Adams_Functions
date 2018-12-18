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
%       varargin    - 'MinRowNumber': minimum row number for featuresAll
%                   must be a positive integer scalar
%                   default = 1
%                   - 'LowPassFrequency': frequency of lowpass filter in Hz
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
%                   - 'SaveFlag': whether to save the pulse train series
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/compute_sampling_interval.m
%       cd/count_vectors.m
%       cd/create_labels_from_numbers.m
%       cd/parse_pulse_response.m
%
% Used by:
%       cd/plot_protocols.m

% File History:
% 2018-12-15 Modified from compute_average_pulse_response.m
% 2018-12-17 Now uses create_labels_from_numbers.m
% TODO: save table if SaveFlag is true
% 

%% Hard-coded parameters
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Other'};

%% Default values for optional arguments
minRowNumberDefault = 1;        % row number to start counting at is 1
lowPassFrequencyDefault = [];   % do not lowpass filter by default
responseLengthMsDefault = 20;   % a response of 20 ms by default
baselineLengthMsDefault = 5;    % a baseline of 5 ms by default
minPeakDelayMsDefault = 0;      % no minimum peak delay by default
channelTypesDefault = {};       % set later
channelUnitsDefault = {};       % set later
channelLabelsDefault = {};      % set later
parsedParamsDefault = [];       % set later
parsedDataDefault = [];         % set later
saveFlagDefault = true;         % save the pulse train series by default

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
addParameter(iP, 'MinRowNumber', minRowNumberDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
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
addParameter(iP, 'SaveFlag', saveFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, fileName, responseType, varargin{:});
minRowNumber = iP.Results.MinRowNumber;
baselineLengthMs = iP.Results.BaselineLengthMs;
minPeakDelayMs = iP.Results.MinPeakDelayMs;
saveFlag = iP.Results.SaveFlag;

%% Filter and extract pulse response(s)
[tVecAll, respAll, stimAll, labels] = ...
    filter_and_extract_pulse_response(fileName, responseType, varargin{:});

%% Compute features of the average pulse response
% Compute the sampling interval(s)
siMs = compute_sampling_interval(tVecAll);

% Count the number of vectors
nVectors = count_vectors(tVecAll);

% Parse the average pulse response vector
featuresAll = parse_pulse_response(respAll, siMs, ...
                                'PulseVectors', stimAll, ...
                                'SameAsPulse', true, ...
                                'MeanValueWindowMs', baselineLengthMs, ...
                                'MinPeakDelayMs', minPeakDelayMs);

% Add labels to each row
labels = repmat({labels}, [nVectors, 1]);
featuresAll = addvars(featuresAll, labels);

% Construct sweep names
sweepName = create_labels_from_numbers(1:nVectors, 'Prefix', 'Swp', ...
                                    'ForceColumnOutput', true);

% Repeat the file names
filePath = repmat({fileName}, [nVectors, 1]);

% Insert the sweep name before the first column of featuresAll
featuresAll = addvars(featuresAll, sweepName, 'Before', 1);

% Insert the file path before the sweep name of featuresAll
featuresAll = addvars(featuresAll, filePath, 'Before', 1);

%% Save the features table
% writetable(featuresAll, sheetPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Insert before the second column of featuresAll the median time

if iscell(tVecAll)
    siMs = tVecAll{1}(2) - tVecAll{1}(1);
else
    siMs = tVecAll(2) - tVecAll(1);
end

sweepName = arrayfun(@(x) ['Swp', num2str(x)], ...
                    transpose(1:nVectors), 'UniformOutput', false);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
