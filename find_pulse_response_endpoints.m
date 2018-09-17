function [idxCprStart, idxCprEnd, isUnbalanced, idxCpStart, idxCpEnd] = ...
            find_pulse_response_endpoints (vvecCpr, siMs, varargin)
%% Computes the average initial slope from a current pulse response
% Usage: [idxCprStart, idxCprEnd, isUnbalanced] = ...
%           find_pulse_response_endpoints (vvecCpr, siMs, varargin)
%
% Arguments:    
%       vvecCpr     - voltage vector of the current pulse response
%                   must be a numeric vector
%       siMs        - sampling interval in ms
%                   must be a positive scalar
%       varargin    - 'IvecCpr': current vector of the current pulse response
%                   must be a numeric vector
%                   default == [] (not used)
%                   - 'SameAsIvec': whether always the same as 
%                                       the current pulse endpoints
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'CprLengthMs': length of the current pulse response
%                                       after pulse endpoint in ms
%                   must be a nonnegative scalar
%                   default = 20 ms
%
% Requires:
%       /home/Matlab/Adams_Functions/find_first_jump.m
%       /home/Matlab/Adams_Functions/find_pulse_endpoints.m
%
% Used by:    
%       /home/Matlab/Adams_Functions/compute_average_initial_slopes.m
%       /home/Matlab/Adams_Functions/plot_evoked_LFP.m

% File History:
% 2018-08-13 AL - Adapted from compute_average_initial_slopes.m
% 2018-09-17 AL - Changed required arguement tVecCpr to siMs
% 2018-09-17 AL - Added optional parameters SameAsIvec and CprLengthMs

%% Hard-coded parameters
signal2Noise = 10;
noiseWindowSize = 5;

%% Default values for optional arguments
ivecCprDefault = [];            % don't use current vector by default
sameAsIvecDefault = true;       % use current pulse endpoints by default
cprLengthMsDefault = 20;        % a response of 20 ms by default

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
addRequired(iP, 'vvecCpr', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'siMs', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'IvecCpr', ivecCprDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'SameAsIvec', sameAsIvecDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'CprLengthMs', cprLengthMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));

% Read from the Input Parser
parse(iP, vvecCpr, siMs, varargin{:});
ivecCpr = iP.Results.IvecCpr;
sameAsIvec = iP.Results.SameAsIvec;
cprLengthMs = iP.Results.CprLengthMs;

%% Do the job
% Compute the length of the current pulse response in samples
cprLengthSamples = floor(cprLengthMs / siMs);

% Find the start and end points of the current pulse
if ~isempty(ivecCpr)
    % If provided, use the current trace
    [idxCpStart, idxCpEnd] = find_pulse_endpoints(ivecCpr);

    % If an index is not detected, return warning message
    % TODO: fprintf('The current pulse could not be detected!\n');
else
    % If not, estimate by smoothing the trace, then look for inflection points 
    % TODO
end

% Use windows straddling the the start/end points of the current pulse
%   as regions for finding the start/end points of the current pulse response
%   Note: this will always cause idxCpStart/idxCpEnd to be the first
%           index for checking, and will check noiseWindowSize more points
idxRegion1Start = idxCpStart - noiseWindowSize;
idxRegion1End = idxCpStart + noiseWindowSize;
vvecRegion1 = vvecCpr(idxRegion1Start:idxRegion1End);

idxRegion2Start = idxCpEnd - noiseWindowSize;
idxRegion2End = idxCpEnd + noiseWindowSize;
vvecRegion2 = vvecCpr(idxRegion2Start:idxRegion2End);

% Initialize isUnbalanced as false
isUnbalanced = false;

% Find the start point of the current pulse response
%   by detecting the 'first jump' in the region of interest #1
%   If it doesn't exist, use the start point of the current pulse
[~, idxTemp1] = ...
    find_first_jump(vvecRegion1, 'NSamplesPerJump', 2, ...
                                 'Signal2Noise', signal2Noise, ...
                                 'NoiseWindowSize', noiseWindowSize);
if ~isempty(idxTemp1) && ~sameAsIvec
    idxCprStart = (idxRegion1Start - 1) + idxTemp1;
    isUnbalanced = true;
else
    idxCprStart = idxCpStart;
end

% Find the end point of the current pulse response
%   by detecting the 'first jump' in the region of interest #2
%   If it doesn't exist, use the end point of the current pulse
%       plus cprLengthSamples
[~, idxTemp2] = ...
    find_first_jump(vvecRegion2, 'NSamplesPerJump', 2, ...
                                 'Signal2Noise', signal2Noise, ...
                                 'NoiseWindowSize', noiseWindowSize);
if ~isempty(idxTemp2) && ~sameAsIvec
    idxCprEnd = (idxRegion2Start - 1) + idxTemp2 + cprLengthSamples;
    isUnbalanced = true;
else
    idxCprEnd = idxCpEnd + cprLengthSamples;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Compute slope right after current pulse start
idxFirst1 = idxCpStart;
idxLast1 = idxCpStart + nSamples - 1;
startSlope1 = compute_slope(tvecCpr, vvecCpr, idxFirst1, idxLast1);

% Compute slope right after current pulse end
idxFirst2 = idxCpEnd;
idxLast2 = idxCpEnd + nSamples - 1;
endSlope2 = compute_slope(tvecCpr, vvecCpr, idxFirst2, idxLast2);

% Compute slope right after current pulse start
idxFirst3 = idxCpStart2;
idxLast3 = idxCpStart2 + nSamples - 1;
startSlope3 = compute_slope(tvecCpr, vvecCpr, idxFirst3, idxLast3);

% Compute slope right after current pulse end
idxFirst4 = idxCpEnd2;
idxLast4 = idxCpEnd2 + nSamples - 1;
endSlope4 = compute_slope(tvecCpr, vvecCpr, idxFirst4, idxLast4);

% Choose the more negative of the start slopes 
%  and the more positive of the end slopes
startSlope = min([startSlope1, startSlope3]);
endSlope = max([endSlope2, endSlope4]);

addRequired(iP, 'nSamples', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
parse(iP, tvecCpr, vvecCpr, ivecCpr, nSamples);

function [avgSlope, startSlope, endSlope, indsUsedForPlot] = find_pulse_response_endpoints (tvecCpr, vvecCpr, ivecCpr, nSamples)

% Note: function calls are slower
%       /home/Matlab/Adams_Functions/compute_slope.m
startSlope = compute_slope(tvecCpr, vvecCpr, idxFirst1, idxLast1);
endSlope = compute_slope(tvecCpr, vvecCpr, idxFirst2, idxLast2);

% Crop the voltage trace
vvecCropped = vvecCpr((idxCprStart + 1):end);

%       tvecCpr     - time vector of the current pulse response in ms
%                   must be a numeric vector
addRequired(iP, 'tvecCpr', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
parse(iP, tvecCpr, vvecCpr, varargin{:});

%}

