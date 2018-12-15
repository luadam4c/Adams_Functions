function iVecs = extract_current_vectors (abfFileName, varargin)
%% Extracts current vectors from a .abf file
% Usage: iVecs = extract_current_vectors (abfFileName, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       iVecs       - any identified current vector(s)
%                       (Note: 2nd dimension: sweep; 
%                           optional 3rd dimension: channel)
%                   specified as a numeric array
% Arguments:
%       abfFileName - file name could be either the full path or 
%                       a relative path in current directory
%                       .abf is not needed (e.g. 'B20160908_0004')
%                   must be a string scalar or a character vector
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

%% Default values for optional arguments
channelTypesDefault = {};       % set later
parsedDataDefault = [];         % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'abfFileName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ChannelTypes', channelTypesDefault, ...
    @(x) validateattributes(x, {'cell'}, {'nonempty'}));
addParameter(iP, 'ParsedData', parsedDataDefault, ...
    @(x) validateattributes(x, {'struct'}, {'2d'}));

% Read from the Input Parser
parse(iP, abfFileName, varargin{:});
channelTypes = iP.Results.ChannelTypes;
parsedData = iP.Results.ParsedData;

%% Preparation
% Parse the abf file if parsedData not provided
if isempty(parsedData)
    [~, parsedData] = parse_abf(abfFileName, 'Verbose', false, ...
                                'ChannelTypes', channelTypes);
end

%% Do the job
% Extract the current vectors
iVecs = parsedData.iVecs;

% If iVecs is a cellarray, use the first element
if iscell(iVecs)
    iVecs = iVecs{1};
end

% If iVecs is 3-D, use the first two dimensions 
if ~ismatrix(iVecs)
    iVecs = squeeze(iVecs(:, :, 1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

[~, ~, ~, ~, iVecs, ~] = ...
    parse_abf(abfFileName, 'Verbose', false, ...
                'ChannelTypes', channelTypes);

if ndims(iVecs) > 2

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%