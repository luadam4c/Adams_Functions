function isRepetitivePulse = identify_repetitive_pulses (vectors, varargin)
%% Identifies whether a set of vectors are repetitive pulses
% Usage: isRepetitivePulse = identify_repetitive_pulses (vectors, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       isRepetitivePulse   -  whether the vectors consists of repetitive pulses
%                           specified as a logical scalar
% Arguments:
%       vectors     - vectors containing a pulse
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%       varargin    - 'MinSweeps': minimum number of vectors 
%                   must be a positive integer scalar
%                   default == 2
%
% Requires:
%       cd/count_vectors.m
%       cd/parse_pulse.m
%       cd/relative_std.m
%
% Used by:
%       cd/identify_eLFP_protocol.m
%       cd/identify_gabab_protocol.m

% File History:
% 2018-12-15 Moved from identify_eLFP_protocol.m
% 2018-12-15 Now uses count_vectors.m
% 2018-12-15 Now uses a relative standard deviation (%) threshold

%% Hard-coded parameters
rsdThreshold = 1;               % relative standard deviation (%) threshold

%% Default values for optional arguments
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
addRequired(iP, 'vectors', ...                   % vectors
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MinSweeps', minSweepsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));

% Read from the Input Parser
parse(iP, vectors, varargin{:});
minSweeps = iP.Results.MinSweeps;

%% Preparation
% Count the number of vectors
nVectors = count_vectors(vectors);

%% Do the job
% Not repetitive pulses if there are too few vectors
if nVectors < minSweeps
    isRepetitivePulse = false;
    return
end

% Parse the pulse(s)
parsedParams = parse_pulse(vectors);

% Extract the pulse response start indices and amplitudes
idxCpStarts = parsedParams.idxBeforeStart;
ampCps = parsedParams.pulseAmplitude;

% Check whether the variation of amplitudes and starting indices
%   are small enough by the relative standard deviation measure
if relative_std(idxCpStarts) < rsdThreshold && ...
    relative_std(ampCps) < rsdThreshold
    isRepetitivePulse = true;
else
    isRepetitivePulse = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Identify the pulse response endpoints and midpoints
idxCpStarts = zeros(nVectors, 1);
idxCpEnds = zeros(nVectors, 1);
idxCpMids = zeros(nVectors, 1);
ampCps = zeros(nVectors, 1);
parfor iSwp = 1:nVectors
    % Identify the pulse endpoints
    [idxCpStart, idxCpEnd] = ...
        find_pulse_endpoints(iVecs(:, iSwp));

    if isempty(idxCpStart) || isempty(idxCpEnd)
        idxCpStart = NaN;
        idxCpEnd = NaN;
        idxCpMid = NaN;
        ampCp = NaN;
    else
        % Identify the pulse midpoints
        idxCpMid = ceil(mean([idxCpStart, idxCpEnd]));

        % Identify the pulse amplitudes
        ampCp = max(abs(iVecs(:, iSwp)));
    end
    
    % Store in arrays
    idxCpEnds(iSwp) = idxCpEnd;
    idxCpStarts(iSwp) = idxCpStart;
    idxCpMids(iSwp) = idxCpMid;
    ampCps(iSwp) = ampCp;
end

% Identify the pulse response endpoints and midpoints
[idxCpStarts, idxCpEnds] = find_pulse_endpoints(iVecs);

% Identify the pulse response midpoints
idxCpMids = ceil(mean([idxCpStarts, idxCpEnds]), 2);

idxCpEnds = parsedParams.idxBeforeEnd;
idxCpMids = parsedParams.idxMidpoint;

covThreshold = 0.01;            % coefficient of variation threshold

minSweeps = 2;                  % must have at least 2 sweeps

% Count the number of vectors
nVectors = size(vectors, 2);

% If there are no sweeps, this is not an evoked local field potential
%   protocol
if nVectors == 0
    isRepetitivePulse = false;
    return
end

% Not repetitive pulses if there are no vectors recorded
if isempty(vectors)
    isRepetitivePulse = false;
    return;
end

if nanstd(ampCps) / abs(nanmean(ampCps)) < covThreshold && ...
    nanstd(idxCpStarts) / abs(nanmean(idxCpStarts)) < covThreshold

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%