function [idxStartAll, idxEndAll, idxStart2All, idxEnd2All] = ...
            find_pulse_endpoints (vectors)
%% Returns the start and end indices of the first pulse from vector(s)
% Usage: [idxStartAll, idxEndAll, idxStart2All, idxEnd2All] = ...
%           find_pulse_endpoints (vectors)
% Outputs:
%       idxStartAll - indices of pulse start (right before pulse start)
%                   specified as a positive integer column vector
%       idxEndAll   - indices of pulse end (right before pulse end)
%                   specified as a positive integer column vector
%       idxStart2All- indices of pulse start (right after pulse start)
%                   specified as a positive integer column vector
%       idxEnd2All  - indices of pulse end (right after pulse end)
%                   specified as a positive integer column vector
% Arguments:
%       vectors     - vectors containing a pulse
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%
% Requires:
%       
% Used by:
%       cd/count_vectors.m
%       cd/parse_pulse.m
%       cd/find_pulse_response_endpoints.m
%
% File History:
% 2018-07-25 BT - Adapted from find_initial_slopes.m
% 2018-08-10 AL - Change the amplitude to take the value from pulseShifted
%                   rather than from pulse
% 2018-08-10 AL - Now checks number of arguments
% 2018-09-17 AL - Now returns empty indices if there is no pulse
% 2018-10-09 AL - Improved documentation
% 2018-10-09 AL - Now accepts multiple vectors an array or a cell array
% 2018-10-10 AL - Added idxStart2All and idxEnd2All

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
    @(x) isnumeric(x) || iscell(x) && all(cellfun(@isnumeric, x)) );

% Read from the Input Parser
parse(iP, vectors);

%% Preparation
% Count the number of vectors
nVectors = count_vectors(vectors);

%% Do the job
idxStartAll = zeros(nVectors, 1);
idxEndAll = zeros(nVectors, 1);
idxStart2All = zeros(nVectors, 1);
idxEnd2All = zeros(nVectors, 1);
if iscell(vectors)
    %parfor iVec = 1:nVectors
    for iVec = 1:nVectors
        [idxStartAll(iVec), idxEndAll(iVec), ...
            idxStart2All(iVec), idxEnd2All(iVec)] = ...
            find_pulse_endpoints_helper(vectors{iVec});
    end
elseif isnumeric(vectors)
    parfor iVec = 1:nVectors
        [idxStartAll(iVec), idxEndAll(iVec), ...
            idxStart2All(iVec), idxEnd2All(iVec)] = ...
            find_pulse_endpoints_helper(vectors(:, iVec));
    end
else
    error('vectors is not the right type!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idxStart, idxEnd, idxStart2, idxEnd2] = ...
                find_pulse_endpoints_helper(vector)
%% Finds pulse endpoints from a single vector

% Count the number of samples
nSamples = length(vector);

% Subtract the trace by the initial value
vectorShifted = vector - vector(1);

% Find the maximum absolute value and make that the pulse amplitude
[~, idxAbsMax] = max(abs(vectorShifted));
amplitude = vectorShifted(idxAbsMax);

% Find the sign of the amplitude
signAmplitude = sign(amplitude);

% Find the indices of the first half versus the second half
indHalf1 = 1:idxAbsMax;
indHalf2 = idxAbsMax:nSamples;

% Find the indices of the start and end of the pulse
%   Note: Change search direction based on positive or negative pulse
if signAmplitude == 1
    % Find the last point less than 1/4 of the amplitude in the first half
    idxStart = find(vectorShifted(indHalf1) < amplitude * 0.25, 1, 'last');

    % Find the first point greater than 3/4 of the amplitude in the first half
    idxStart2 = find(vectorShifted(indHalf1) > amplitude * 0.75, 1);

    % Find the last point greater than 3/4 of the amplitude in the second half
    idxEndRel = find(vectorShifted(indHalf2) > amplitude * 0.75, 1, 'last');

    % Find the first point less than 1/4 of the amplitude in the second half
    idxEnd2Rel = find(vectorShifted(indHalf2) < amplitude * 0.25, 1);
elseif signAmplitude == -1
    % Find the last point greater than 1/4 of the amplitude in the first half
    idxStart = find(vectorShifted(indHalf1) > amplitude * 0.25, 1, 'last');

    % Find the first point less than 3/4 of the amplitude in the first half
    idxStart2 = find(vectorShifted(indHalf1) < amplitude * 0.75, 1);

    % Find the last point less than 3/4 of the amplitude in the second half
    idxEndRel = find(vectorShifted(indHalf2) < amplitude * 0.75, 1, 'last');

    % Find the first point greater than 1/4 of the amplitude in the second half
    idxEnd2Rel = find(vectorShifted(indHalf2) > amplitude * 0.25, 1);
else
    idxStart = [];
    idxStart2 = [];
    idxEnd = [];
    idxEnd2 = [];
    return;
end

% Shift the indices to correspond to entire vector
idxEnd = idxEndRel + idxAbsMax - 1;
idxEnd2 = idxEnd2Rel + idxAbsMax - 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Find the pulse amplitude at that point
amplitude = pulse(idxAbsMax);

% Force vectors to be a column cell array
if iscell(vectors)
    vectors = vectors(:);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%