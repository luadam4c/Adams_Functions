function [isEvokedLfp] = identify_eLFP (iVecsORfileName, varargin)
%% Identifies whether an abf file follows an eLFP protocol
% Usage: [isEvokedLfp] = identify_eLFP (iVecsORfileName, varargin)
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
%
% Requires:
%       cd/parse_abf.m
%
% Used by:    
%       cd/plot_all_abfs_dir.m

% File History:
% 2018-09-17 Created by Adam Lu
% 2018-09-21 Considered the case when iVecs is a cellarray
% 2018-09-21 Considered the case when iVecs is 3-D
% 

%% Hard-coded parameters
coeffVarThreshold = 0.01;       % coefficient of variation threshold
minSweeps = 2;                  % must have at least 2 sweeps
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Undefined'};

%% Default values for optional arguments
channelTypesDefault = {};       % set later

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

% Read from the Input Parser
parse(iP, iVecsORfileName, varargin{:});
channelTypes = iP.Results.ChannelTypes;

% Validate channel types
if ~isempty(channelTypes)
    channelTypes = cellfun(@(x) validatestring(x, validChannelTypes), ...
                            channelTypes, 'UniformOutput', false);
end

%% Preparation
if ischar(iVecsORfileName) || isstring(iVecsORfileName)
    % The first arguement is a file name
    fileName = iVecsORfileName;

    % Parse the abf file to get the current vectors
    [~, ~, ~, ~, iVecs, ~] = ...
        parse_abf(fileName, 'Verbose', false, ...
                    'ChannelTypes', channelTypes);

    % If iVecs is a cellarray, use the first element
    if iscell(iVecs)
        iVecs = iVecs{1};
    end

    % If iVecs is 3-D, use the first two dimensions 
    if ndims(iVecs) > 2
        iVecs = squeeze(iVecs(:, :, 1));
    end
else
    % The first arguement are the current vectors
    iVecs = iVecsORfileName;
end

% Count the number of sweeps
nSweeps = size(iVecs, 2);

%% Do the job
% Not an evoked LFP protocol if there are no current vectors recorded
if isempty(iVecs)
    isEvokedLfp = false;
    return;
end

% Not an evoked LFP protocol if there are too few current vectors
if nSweeps < minSweeps
    isEvokedLfp = false;
    return;
end

% Identify the current pulse response endpoints and midpoints
idxCpStarts = zeros(nSweeps, 1);
idxCpEnds = zeros(nSweeps, 1);
idxCpMids = zeros(nSweeps, 1);
ampCps = zeros(nSweeps, 1);
parfor iSwp = 1:nSweeps
%for iSwp = 1:nSweeps
    % Identify the current pulse endpoints
    [idxCpStart, idxCpEnd] = ...
        find_pulse_endpoints(iVecs(:, iSwp));

    if isempty(idxCpStart) || isempty(idxCpEnd)
        idxCpStart = NaN;
        idxCpEnd = NaN;
        idxCpMid = NaN;
        ampCp = NaN;
    else
        % Identify the current pulse midpoints
        idxCpMid = ceil(mean([idxCpStart, idxCpEnd]));

        % Identify the current pulse amplitudes
        ampCp = max(abs(iVecs(:, iSwp)));
    end
    
    % Store in arrays
    idxCpEnds(iSwp) = idxCpEnd;
    idxCpStarts(iSwp) = idxCpStart;
    idxCpMids(iSwp) = idxCpMid;
    ampCps(iSwp) = ampCp;
end

% Check whether the variation of amplitudes and starting indices
%   are small enough
if nanstd(ampCps) / abs(nanmean(ampCps)) < coeffVarThreshold && ...
    nanstd(idxCpStarts) / abs(nanmean(idxCpStarts)) < coeffVarThreshold
    isEvokedLfp = true;
else
    isEvokedLfp = false;
end

disp('done');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

ampCp = iVecs(idxCpMid, iSwp);

%}
