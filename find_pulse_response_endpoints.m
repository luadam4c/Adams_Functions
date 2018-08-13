function [idxCprStart, idxCprEnd, isUnbalanced] = ...
            find_pulse_response_endpoints (tvecCpr, vvecCpr, varargin)
%% Computes the average initial slope from a current pulse response
% Usage: [idxCprStart, idxCprEnd, isUnbalanced] = ...
%           find_pulse_response_endpoints (tvecCpr, vvecCpr, varargin)
%
% Arguments:    
%       tvecCpr     - time vector of the current pulse response
%                   must be a numeric vector
%       vvecCpr     - voltage vector of the current pulse response
%                   must be a numeric vector
%       varargin    - 'IvecCpr': current vector of the current pulse response
%                   must be a numeric vector
%                   default == [] (not used)
%                   
%
% Requires:
%       /home/Matlab/Adams_Functions/find_first_jump.m
%       /home/Matlab/Adams_Functions/find_pulse_endpoints.m
%
% Used by:    
%       /home/Matlab/Adams_Functions/compute_average_initial_slopes.m

% File History:
% 2018-08-13 AL - Adapted from compute_average_initial_slopes.m

%% Hard-coded parameters
signal2Noise = 10;
noiseWindowSize = 5;

%% Default values for optional arguments
defaultIvecCpr = [];                % don't use current vector by default

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
addRequired(iP, 'tvecCpr', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'vvecCpr', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'IvecCpr', defaultIvecCpr, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Read from the Input Parser
parse(iP, tvecCpr, vvecCpr, varargin{:});
ivecCpr = iP.Results.IvecCpr;

%% Do the job
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

% Find the start/end point of the current pulse response
%   by detecting the 'first jump' in the region of interest
%   If it doesn't exist, use the start/end point of the current pulse
[~, idxTemp1] = ...
    find_first_jump(vvecRegion1, 'NSamplesPerJump', 2, ...
                             'Signal2Noise', signal2Noise, ...
                             'NoiseWindowSize', noiseWindowSize);
if ~isempty(idxTemp1)
    idxCprStart = (idxRegion1Start - 1) + idxTemp1;
    isUnbalanced = true;
else
    idxCprStart = idxCpStart;
end
[~, idxTemp2] = ...
    find_first_jump(vvecRegion2, 'NSamplesPerJump', 2, ...
                                 'Signal2Noise', signal2Noise, ...
                                 'NoiseWindowSize', noiseWindowSize);
if ~isempty(idxTemp2)
    idxCprEnd = (idxRegion2Start - 1) + idxTemp2;
    isUnbalanced = true;
else
    idxCprEnd = idxCpEnd;
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


%}

