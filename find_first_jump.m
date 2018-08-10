function [valJump, idxJump] = find_first_jump (vector, varargin)
%% Finds the index of the first jump in a time series
% Usage: [valJump, idxJump] = find_first_jump (vector, varargin)
% Outputs:
%       valJump     - value of the jump 
%                       (large difference between consecutive points)
%                   specified as a numeric scalar
%       idxJump     - starting index of the jump
%                   specified as a positive integer scalar
% Arguments:    
%       vector      - a time series
%                   must be a numeric vector
%       varargin    - 'Signal2Noise': signal to noise ratio
%                   must be a positive scalar
%                   default == 2
%                   - 'NoiseWindowSize': noise window size in samples
%                   must be a positive integer scalar
%                   default == 5 samples
%
%
% Requires:
%       
% Used by:
%       /home/Matlab/Brians_Functions/compute_average_initial_slopes.m

% File History:
% 2018-08-10 Created by Adam Lu

%% Default values for optional arguments
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
addParameter(iP, 'Signal2noise', signal2noiseDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'NoiseWindowSize', noiseWindowSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Read from the Input Parser
parse(iP, vector, varargin{:});
signal2noise = iP.Results.Signal2noise;
noiseWindowSize = iP.Results.NoiseWindowSize;

%% Do the job
% Get all differences between consecutive samples
allDiffs = diff(vector);

% Get the number of jumps
nJumps = length(allDiffs);

% Get the starting index of the first jump to test
idxFirstJumpToTest = noiseWindowSize + 1;

% Iterate through all noise windows
for idxJump = idxFirstJumpToTest:nJumps
    % Get the indices of the noise window
    indNoise = (idxJump - noiseWindowSize):(idxJump - 1);

    % Compute the root-mean-square averaged 'noise'
    avgNoise = rms(allDiffs(indNoise));

    % Get the value of the jump
    valJump = allDiffs(idxJump);

    % If the jump magnitude is more than signal2noise of avgNoise, stop
    if abs(valJump) > avgNoise
        break
    elseif idxJump == nJumps
        idxJump = [];
        valJump = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:


%}
