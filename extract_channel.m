function vectors = extract_channel (abfFileName, channelType, varargin)
%% Extracts vectors of a given type from a .abf file
% Usage: vectors = extract_channel (abfFileName, channelType, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       vectors     - any identified current vector(s) of the given type
%                       (Note: 2nd dimension: sweep; 
%                           optional 3rd dimension: channel)
%                   specified as a numeric array
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
%       varargin    - 'ChannelTypes': type assigned to each channel, possibly:
%                           'Voltage', 'Current' or 'Conductance'
%                   must as a row cell array with the
%                        number of elements same as the length of the 
%                        2nd dimension of abfdata
%                   - 'ParsedData': parsed data returned by parse_abf.m
%                   must be a scalar structure
%                   default == what the file provides
%
% Requires:
%       cd/parse_abf.m
%
% Used by:
%       cd/identify_eLFP.m

% File History:
% 2018-12-15 Moved from identify_eLFP.m
% 2018-12-15 Added 'ParsedData' as an optional parameter
% 2018-12-15 Now uses ismatrix() per MATLAB's suggestion
% TODO: Generalize to other file types
% 

%% Hard-coded parameters
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Other'};

%% Default values for optional arguments
channelTypesDefault = {};       % set later
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
addRequired(iP, 'abfFileName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'channelType', ...
    @(x) any(validatestring(x, validChannelTypes)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ChannelTypes', channelTypesDefault, ...
    @(x) validateattributes(x, {'cell'}, {'nonempty'}));
addParameter(iP, 'ParsedData', parsedDataDefault, ...
    @(x) validateattributes(x, {'struct'}, {'2d'}));

% Read from the Input Parser
parse(iP, abfFileName, channelType, varargin{:});
channelTypes = iP.Results.ChannelTypes;
parsedData = iP.Results.ParsedData;

% Validate the channel type
channelType = validatestring(channelType, validChannelTypes);

%% Preparation
% Parse the abf file if parsedData not provided
if isempty(parsedData)
    [~, parsedData] = parse_abf(abfFileName, 'Verbose', false, ...
                                'ChannelTypes', channelTypes, ...
                                'ExtractChannels', true);
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

% If vectors is a cellarray, use the first element
if iscell(vectors)
    vectors = vectors{1};
end

% If vectors is 3-D, use the first two dimensions 
if ~ismatrix(vectors)
    vectors = squeeze(vectors(:, :, 1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

[~, ~, ~, ~, iVecs, ~] = ...
    parse_abf(abfFileName, 'Verbose', false, ...
                'ChannelTypes', channelTypes);

if ndims(vectors) > 2

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%