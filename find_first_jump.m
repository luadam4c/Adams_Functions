function [valJump, idxJumpStart] = find_first_jump (vector, varargin)
%% Finds the index of the first jump in a time series
% Usage: [valJump, idxJump] = find_first_jump (vector, varargin)
% Explanation:
%       TODO
%
% Examples:
%       TODO
%
% Outputs:
%       valJump     - value of the jump 
%                       (large difference between consecutive points)
%                   specified as a numeric scalar
%       idxJumpStart - starting index of the jump
%                   specified as a positive integer scalar
%
% Arguments:    
%       vector      - a time series
%                   must be a numeric vector
%       varargin    - 'NSamplesPerJump': number of samples per jump
%                   must be a positive integer scalar
%                   default == 2
%                   - 'Signal2Noise': signal to noise ratio
%                   must be a positive scalar
%                   default == 2
%                   - 'NoiseWindowSize': noise window size in samples
%                   must be a positive integer scalar
%                   default == 5 samples
%
%
% Requires:
%       cd/find_first_deviant.m
%       
% Used by:
%       cd/find_pulse_response_endpoints.m
%       cd/parse_lts.m

% File History:
% 2018-08-10 Created by Adam Lu

%% Default values for optional arguments
nSamplesPerJumpDefault = 2;     % default is consecutive samples
signal2noiseDefault = 2;        % default signal to noise ratio
noiseWindowSizeDefault = 5;     % default noise window in samples

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
addRequired(iP, 'vector', ...                  % voltage vector
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'NSamplesPerJump', nSamplesPerJumpDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'Signal2noise', signal2noiseDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'NoiseWindowSize', noiseWindowSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Read from the Input Parser
parse(iP, vector, varargin{:});
nSamplesPerJump = iP.Results.NSamplesPerJump;
signal2noise = iP.Results.Signal2noise;
noiseWindowSize = iP.Results.NoiseWindowSize;

%% Do the job
% Get all differences between consecutive samples
if nSamplesPerJump == 2
    allDiffs = diff(vector);
else
    allDiffs = vector(nSamplesPerJump:end) - vector(1:end-nSamplesPerJump+1);
end

% Find the first deviating difference
[valJump, idxJumpStart] = ...
    find_first_deviant(allDiffs, 'Deviant2Peers', signal2noise, ...
                                'PeersWindowSize', noiseWindowSize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%