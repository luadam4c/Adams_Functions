function isEvokedLfp = identify_eLFP (iVecsORfileName, varargin)
%% Identifies whether an abf file follows an eLFP protocol
% Usage: isEvokedLfp = identify_eLFP (iVecsORfileName, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       isEvokedLfp - whether the data corresponds to an evoked LFP protocol
%                   specified as a logical scalar
%
% Arguments:    
%       iVecsORfileName
%                   - current vector(s) from a .abf file OR
%                       .abf file name (could be either the full path or 
%                       a relative path in current directory
%                       .abf is not needed (e.g. 'B20160908_0004'))
%                   must be a numeric array
%       varargin    - 'ChannelTypes': type assigned to each channel, possibly:
%                           'Voltage', 'Current' or 'Conductance'
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
% 2018-09-17 - Created by Adam Lu
% 2018-09-21 - Considered the case when iVecs is a cellarray
% 2018-09-21 - Considered the case when iVecs is 3-D
% 2018-10-03 - Updated usage of parse_abf.m
% 2018-12-15 - Made 'MinSweeps' an optional parameter with default 2
% 2018-12-15 - Moved code to identify_repetitive_pulses.m
% 2018-12-15 - Moved code to extract_channel.m
% 

%% Hard-coded parameters
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Undefined'};

%% Default values for optional arguments
channelTypesDefault = {};       % set later
minSweepsDefault = 2;           % must have at least 2 sweeps by default

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
addRequired(iP, 'iVecsORfileName');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ChannelTypes', channelTypesDefault, ...
    @(x) validateattributes(x, {'cell'}, {'nonempty'}));
addParameter(iP, 'MinSweeps', minSweepsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));

% Read from the Input Parser
parse(iP, iVecsORfileName, varargin{:});
channelTypes = iP.Results.ChannelTypes;
minSweeps = iP.Results.MinSweeps;

% Validate channel types
if ~isempty(channelTypes)
    channelTypes = cellfun(@(x) validatestring(x, validChannelTypes), ...
                            channelTypes, 'UniformOutput', false);
end

%% Preparation
% Parse the first argument and extract the current vectors if needed
if ischar(iVecsORfileName) || isstring(iVecsORfileName)
    % The first argument is a file name
    fileName = iVecsORfileName;

    % Extract the current vectors
    iVecs = extract_channel(fileName, 'current', 'ChannelTypes', channelTypes);
else
    % The first argument are the current vectors
    iVecs = iVecsORfileName;
end

%% Do the job
% Identify whether the current vectors are repetitive pulses
isEvokedLfp = identify_repetitive_pulses(iVecs, 'MinSweeps', minSweeps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

ampCp = iVecs(idxCpMid, iSwp);

disp('done');

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%