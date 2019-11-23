function [vectors, label] = extract_channel (abfFileName, channelType, varargin)
%% Extracts vectors of a given type from a .abf file
% Usage: [vectors, label] = extract_channel (abfFileName, channelType, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       vectors     - any identified current vector(s) of the given type
%                       (Note: 2nd dimension: sweep; 
%                           optional 3rd dimension: channel)
%                   specified as a numeric array
%       label       - label(s) for the channel to extract
%                   specified as a character vector 
%                       or a cell array of character vectors
% Arguments:
%       abfFileName - file name could be either the full path or 
%                       a relative path in current directory
%                       .abf is not needed (e.g. 'B20160908_0004')
%                   must be a string scalar or a character vector
%       channelType - the channel type to extract
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'Voltage'       - voltage
%                       'Current'       - current
%                       'Conductance'   - conductance
%                       'Other'         - other un-identified types
%       varargin    - 'MaxNChannels': maximum number of channels to extract
%                   must be a positive integer scalar or Inf
%                   default == Inf
%                   - 'MaxNVectors': maximum number of vectors to extract
%                   must be a positive integer scalar or Inf
%                   default == Inf
%                   - 'ChannelTypes': type assigned to each channel, possibly:
%                           'Voltage', 'Current' or 'Conductance'
%                   must as a row cell array with the
%                        number of elements same as the length of the 
%                        2nd dimension of abfdata
%                   - 'ParsedParams': parsed parameters returned by parse_abf.m
%                   must be a scalar structure
%                   default == what the file provides
%                   - 'ParsedData': parsed data returned by parse_abf.m
%                   must be a scalar structure
%                   default == what the file provides
%
% Requires:
%       cd/match_positions.m
%       cd/parse_abf.m
%
% Used by:
%       cd/filter_and_extract_pulse_response.m
%       cd/identify_eLFP_protocol.m
%       cd/identify_gabab_protocol.m

% File History:
% 2018-12-15 Moved from identify_eLFP_protocol.m
% 2018-12-15 Added 'ParsedData' as an optional parameter
% 2018-12-15 Now uses ismatrix() per MATLAB's suggestion
% 2018-12-15 Now returns the corresponding label(s)
% 2018-12-16 Fixed usage of maxNChannels and added maxNVectors
% TODO: Generalize to other file types
% 

%% Hard-coded parameters
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Other'};

%% Default values for optional arguments
maxNChannelsDefault = Inf;      % set later
maxNVectorsDefault = Inf;       % set later
channelTypesDefault = {};       % set later
channelUnitsDefault = {};       % set later
channelLabelsDefault = {};      % set later
parsedParamsDefault = [];       % set later
parsedDataDefault = [];         % set later

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
addRequired(iP, 'abfFileName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'channelType', ...
    @(x) any(validatestring(x, validChannelTypes)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MaxNChannels', maxNChannelsDefault, ...
    @(x) isinf(x) || ispositiveintegerscalar(x));
addParameter(iP, 'MaxNVectors', maxNVectorsDefault, ...
    @(x) isinf(x) || ispositiveintegerscalar(x));
addParameter(iP, 'ChannelTypes', channelTypesDefault, ...
    @(x) isempty(x) || isstring(x) || iscellstr(x));
addParameter(iP, 'ChannelUnits', channelUnitsDefault, ...
    @(x) isempty(x) || isstring(x) || iscellstr(x));
addParameter(iP, 'ChannelLabels', channelLabelsDefault, ...
    @(x) isempty(x) || isstring(x) || iscellstr(x));
addParameter(iP, 'ParsedParams', parsedParamsDefault, ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));
addParameter(iP, 'ParsedData', parsedDataDefault, ...
    @(x) validateattributes(x, {'struct'}, {'2d'}));

% Read from the Input Parser
parse(iP, abfFileName, channelType, varargin{:});
maxNChannels = iP.Results.MaxNChannels;
maxNVectors = iP.Results.MaxNVectors;
channelTypes = iP.Results.ChannelTypes;
channelUnits = iP.Results.ChannelUnits;
channelLabels = iP.Results.ChannelLabels;
parsedParams = iP.Results.ParsedParams;
parsedData = iP.Results.ParsedData;

% Validate the channel type
channelType = validatestring(channelType, validChannelTypes);

%% Preparation
% Parse the abf file if parsedData not provided
if isempty(parsedParams) || isempty(parsedData)
    [parsedParams, parsedData] = parse_abf(abfFileName, 'Verbose', false, ...
                                            'ChannelTypes', channelTypes, ...
                                            'ChannelUnits', channelUnits, ...
                                            'ChannelLabels', channelLabels, ...
                                            'ExtractChannels', true);
end

% Extract the channel types if needed
if isempty(channelTypes)
    channelTypes = parsedParams.channelTypes;
end

% Extract the channel label if needed
if isempty(channelLabels)
    channelLabels = parsedParams.channelLabels;
end

%% Do the job
% Extract the vectors from the channel type of interest
switch channelType
    case 'Voltage'
        vectors = parsedData.vVecs;
    case 'Current'
        vectors = parsedData.iVecs;
    case 'Conductance'
        vectors = parsedData.gVecs;
    case 'Other'
        vectors = parsedData.otherVecs;
    otherwise
        error('Logic error!');
end

% Restrict to the desired maximum number of channels
% If vectors is 3-D, use the first two dimensions 
%   Note: ismatrix is false if one of the dimensions is zero
if ~ismatrix(vectors) && min(size(vectors)) > 0
    vectors = squeeze(vectors(:, :, 1:maxNChannels));
end

% Count the number of vectors
nVectors = count_vectors(vectors);

% Restrict to the desired maximum number of vectors
if ~isinf(maxNVectors) && nVectors > maxNVectors
    if iscell(vectors)
        % If vectors is a cellarray, use the first maxNVectors elements
        vectors = vectors(1:maxNVectors);
    else
        % If vectors is a numeric array, use the first maxNVectors columns
        vectors = vectors(:, 1:maxNVectors);
    end
end

% Find the corresponding channel label(s)
label = match_positions(channelLabels, channelTypes, channelType, ...
                        'maxNChannels', maxNChannels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

[~, ~, ~, ~, iVecs, ~] = ...
    parse_abf(abfFileName, 'Verbose', false, ...
                'ChannelTypes', channelTypes);

if ndims(vectors) > 2

vectors = squeeze(vectors(:, :, 1));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%