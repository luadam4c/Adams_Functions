function isGababProt = identify_gabab_protocol (vVecsORfileName, varargin)
%% Identifies whether a .abf file or a set of voltage vectors follows a GABA-B IPSC protocol
% Usage: isGababProt = identify_gabab_protocol (vVecsORfileName, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       isGababProt - whether the data corresponds to a GABA-B IPSC protocol
%                   specified as a logical scalar
%
% Arguments:    
%       vVecsORfileName
%                   - voltage vector(s) OR
%                       .abf file name (could be either the full path or 
%                       a relative path in current directory
%                       .abf is not needed (e.g. 'B20160908_0004'))
%                   must be a numeric array 
%                       or a string scalar or a character vector
%       varargin    - 'ChannelTypes': type assigned to each channel, possibly:
%                           'Voltage', 'Current', 'Conductance' or 'Other'
%                   must as a row cell array with the
%                        number of elements same as the length of the 
%                        2nd dimension of abfdata
%                   - 'MinSweeps': minimum number of sweeps 
%                   must be a positive integer scalar
%                   default == 2
%
% Requires:
%       cd/extract_channel.m
%       cd/identify_repetitive_pulses.m
%
% Used by:    
%       cd/parse_abf.m

% File History:
% 2018-12-15 - Adapted from identify_eLFP_protocol.m
% 2018-12-15 - Changed minSweepsDefault to 1
% 

%% Hard-coded parameters
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Undefined'};

%% Default values for optional arguments
channelTypesDefault = {};       % set later
minSweepsDefault = 1;           % must have at least 1 sweep by default
                                %   Note: In 2018_12_11_A for instance 
                                %       there is only one sweep

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
addRequired(iP, 'vVecsORfileName');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ChannelTypes', channelTypesDefault, ...
    @(x) validateattributes(x, {'cell'}, {'nonempty'}));
addParameter(iP, 'MinSweeps', minSweepsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));

% Read from the Input Parser
parse(iP, vVecsORfileName, varargin{:});
channelTypes = iP.Results.ChannelTypes;
minSweeps = iP.Results.MinSweeps;

% Validate channel types
if ~isempty(channelTypes)
    channelTypes = cellfun(@(x) validatestring(x, validChannelTypes), ...
                            channelTypes, 'UniformOutput', false);
end

%% Preparation
% Parse the first argument and extract the voltage vectors if needed
if ischar(vVecsORfileName) || isstring(vVecsORfileName)
    % The first argument is a file name
    fileName = vVecsORfileName;

    % Extract the voltage vectors
    vVecs = extract_channel(fileName, 'voltage', 'MaxNum', 1, ...
                            'ChannelTypes', channelTypes);
else
    % The first argument are the voltage vectors
    vVecs = vVecsORfileName;

    % Restrict
    if iscell(vVecs) && numel(vVecs) > 1
        vVecs = vVecs{1};
    else
        vVecs = vVecs(:, 1);
    end
end

%% Do the job
% Identify whether the voltage vectors are repetitive pulses
isGababProt = identify_repetitive_pulses(vVecs, 'MinSweeps', minSweeps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%