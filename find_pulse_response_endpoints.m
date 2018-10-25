function [idxResponseStart, idxResponseEnd, hasJump, idxPulseStart, idxPulseEnd] = ...
            find_pulse_response_endpoints (vectors, siMs, varargin)
%% Returns the start and end indices of the first pulse response (from pulse start to 20 ms after pulse ends by default) from vector(s)
% Usage: [idxResponseStart, idxResponseEnd, hasJump, idxPulseStart, idxPulseEnd] = ...
%           find_pulse_response_endpoints (vectors, siMs, varargin)
%
% Outputs:
%       idxResponseStart    - index of pulse response start
%                           specified as a positive integer vector
%       idxResponseEnd      - index of pulse response end
%                           specified as a positive integer vector
%       hasJump             - whether there was a significant "jump" detected
%                           specified as a logical vector
%       idxPulseStart       - index of pulse start
%                           specified as a positive integer vector
%       idxPulseEnd         - index of pulse end
%                           specified as a positive integer vector
%
% Arguments:
%       vectors     - vector(s) containing a pulse response
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%       siMs                - sampling interval in ms
%                           must be a positive scalar
%       varargin    - 'PulseVectors': vectors that contains the pulse itself
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%                   default == [] (not used)
%                   - 'SameAsPulse': whether always the same as 
%                                       the current pulse endpoints
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'ResponseLengthMs': length of the pulse response
%                                           after pulse endpoint in ms
%                   must be a nonnegative scalar
%                   default = 20 ms
%                   - 'BaselineLengthMs': length of the pulse response
%                                           before pulse start in ms
%                   must be a nonnegative scalar
%                   default = 0 ms
%
% Requires:
%       cd/find_first_jump.m
%       cd/find_pulse_endpoints.m
%       cd/iscellnumeric.m
%       cd/match_vector_counts.m
%
% Used by:    
%       cd/compute_initial_slopes.m
%       cd/compute_and_plot_evoked_LFP.m
%       cd/parse_pulse_response.m

% File History:
% 2018-08-13 AL - Adapted from compute_average_initial_slopes.m
% 2018-09-17 AL - Changed required arguement tVecCpr to siMs
% 2018-09-17 AL - Added optional parameters SameAsPulse and ResponseLengthMs
% 2018-09-23 AL - Added the optional parameter baselineLengthSamples
% 2018-10-09 AL - Renamed isUnbalanced -> hasJump and improved documentation

%% Hard-coded parameters
nSamplesPerJump = 2;            % number of samples apart for calculating jump
signal2Noise = 10;              % signal to noise ratio for a jump
noiseWindowSize = 5;            % number of samples to measure baseline noise

%% Default values for optional arguments
pulseVectorsDefault = [];       % don't use pulse vectors by default
sameAsPulseDefault = true;      % use pulse endpoints by default
responseLengthMsDefault = 20;   % a response of 20 ms by default
baselineLengthMsDefault = 0;    % don't include a baseline by default

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
addRequired(iP, 'vectors', ...                   % vectors
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'siMs', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PulseVectors', pulseVectorsDefault, ...
    @(x) isnumeric(x) || iscell(x) && all(cellfun(@isnumeric, x)) );
addParameter(iP, 'SameAsPulse', sameAsPulseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ResponseLengthMs', responseLengthMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'BaselineLengthMs', baselineLengthMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));

% Read from the Input Parser
parse(iP, vectors, siMs, varargin{:});
pulseVectors = iP.Results.PulseVectors;
sameAsPulse = iP.Results.SameAsPulse;
responseLengthMs = iP.Results.ResponseLengthMs;
baselineLengthMs = iP.Results.BaselineLengthMs;

%% Preparation
% Match up pulseVectors with vectors and make sure they are both cell arrays
[pulseVectors, vectors] = ...
    match_vector_counts(pulseVectors, vectors, 'ForceCellOutputs', true);

%% Do the job
[idxResponseStart, idxResponseEnd, hasJump, idxPulseStart, idxPulseEnd] = ...
    cellfun(@(x, y) fpre_helper(x, siMs, y, sameAsPulse, responseLengthMs, ...
        baselineLengthMs, nSamplesPerJump, signal2Noise, noiseWindowSize), ...
        vectors, pulseVectors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idxResponseStart, idxResponseEnd, hasJump, ...
            idxPulseStart, idxPulseEnd] = ...
                fpre_helper (vectors, siMs, pulseVector, sameAsPulse, ...
                                responseLengthMs, baselineLengthMs, ...
                                nSamplesPerJump, signal2Noise, noiseWindowSize)

% Compute the length of the pulse response in samples
cprLengthSamples = floor(responseLengthMs / siMs);

% Compute the baseline length in samples
baselineLengthSamples = floor(baselineLengthMs / siMs);

% Find the start and end points of the pulse
if ~isempty(pulseVector)
    % If provided, use the trace containg the pulse
    [idxPulseStart, idxPulseEnd] = find_pulse_endpoints(pulseVector);

    % If an index is not detected, return warning message
    if isempty(idxPulseStart) || isempty(idxPulseEnd)
        fprintf('The pulse could not be detected!\n');
        idxResponseStart = [];
        idxResponseEnd = [];
        hasJump = [];
        idxPulseStart = [];
        idxPulseEnd = [];
        return
    end
else
    % If not, estimate by smoothing the trace, then look for inflection points 
    % TODO
    error('Please provide a pulseVector for now!\n');
end

% Use windows straddling the the start/end points of the current pulse
%   as regions for finding the start/end points of the current pulse response
%   Note: this will always cause idxPulseStart/idxPulseEnd to be the first
%           index for checking, and will check noiseWindowSize more points
idxRegion1Start = idxPulseStart - noiseWindowSize;
idxRegion1End = idxPulseStart + noiseWindowSize;
responseRegion1 = vectors(idxRegion1Start:idxRegion1End);

idxRegion2Start = idxPulseEnd - noiseWindowSize;
idxRegion2End = idxPulseEnd + noiseWindowSize;
responseRegion2 = vectors(idxRegion2Start:idxRegion2End);

% Initialize hasJump as false
hasJump = false;

% Find the start point of the current pulse response
%   by detecting the 'first jump' in the region of interest #1
%   If it doesn't exist, use the start point of the current pulse
[~, idxTemp1] = ...
    find_first_jump(responseRegion1, 'NSamplesPerJump', nSamplesPerJump, ...
                                 'Signal2Noise', signal2Noise, ...
                                 'NoiseWindowSize', noiseWindowSize);
if ~isempty(idxTemp1) && ~sameAsPulse
    idxResponseStart = (idxRegion1Start - 1) + idxTemp1 - baselineLengthSamples;
    hasJump = true;
else
    idxResponseStart = idxPulseStart - baselineLengthSamples;
end

% Find the end point of the current pulse response
%   by detecting the 'first jump' in the region of interest #2
%   If it doesn't exist, use the end point of the current pulse
%       plus cprLengthSamples
[~, idxTemp2] = ...
    find_first_jump(responseRegion2, 'NSamplesPerJump', nSamplesPerJump, ...
                                 'Signal2Noise', signal2Noise, ...
                                 'NoiseWindowSize', noiseWindowSize);
if ~isempty(idxTemp2) && ~sameAsPulse
    idxResponseEnd = (idxRegion2Start - 1) + idxTemp2 + cprLengthSamples;
    hasJump = true;
else
    idxResponseEnd = idxPulseEnd + cprLengthSamples;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Compute slope right after current pulse start
idxFirst1 = idxPulseStart;
idxLast1 = idxPulseStart + nSamples - 1;
startSlope1 = compute_slope(tvecCpr, vectors, idxFirst1, idxLast1);

% Compute slope right after current pulse end
idxFirst2 = idxPulseEnd;
idxLast2 = idxPulseEnd + nSamples - 1;
endSlope2 = compute_slope(tvecCpr, vectors, idxFirst2, idxLast2);

% Compute slope right after current pulse start
idxFirst3 = idxCpStart2;
idxLast3 = idxCpStart2 + nSamples - 1;
startSlope3 = compute_slope(tvecCpr, vectors, idxFirst3, idxLast3);

% Compute slope right after current pulse end
idxFirst4 = idxCpEnd2;
idxLast4 = idxCpEnd2 + nSamples - 1;
endSlope4 = compute_slope(tvecCpr, vectors, idxFirst4, idxLast4);

% Choose the more negative of the start slopes 
%  and the more positive of the end slopes
startSlope = min([startSlope1, startSlope3]);
endSlope = max([endSlope2, endSlope4]);

addRequired(iP, 'nSamples', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
parse(iP, tvecCpr, vectors, pulseVector, nSamples);

function [avgSlope, startSlope, endSlope, indsUsedForPlot] = find_pulse_response_endpoints (tvecCpr, vectors, pulseVector, nSamples)

% Note: function calls are slower
%       /home/Matlab/Adams_Functions/compute_slope.m
startSlope = compute_slope(tvecCpr, vectors, idxFirst1, idxLast1);
endSlope = compute_slope(tvecCpr, vectors, idxFirst2, idxLast2);

% Crop the voltage trace
vvecCropped = vectors((idxResponseStart + 1):end);

%       tvecCpr     - time vector of the current pulse response in ms
%                   must be a numeric vector
addRequired(iP, 'tvecCpr', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
parse(iP, tvecCpr, vectors, varargin{:});

@(x) validateattributes(x, {'numeric'}, {'vector'}));

%       varargin    - 'PulseVector': vector that contains the pulse itself
%                   must be a numeric vector
%                   default == [] (not used)
pulseVectorDefault = [];        % don't use pulse vector by default
addParameter(iP, 'PulseVector', pulseVectorDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
pulseVector = iP.Results.PulseVector;

idxResponseStart = zeros(nVectors, 1);
idxResponseEnd = zeros(nVectors, 1);
hasJump = zeros(nVectors, 1);
idxPulseStart = zeros(nVectors, 1);
idxPulseEnd = zeros(nVectors, 1);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%