function [indStart1, indEnd1, indStart2, indEnd2] = ...
            find_pulse_endpoints (vectors)
%% Returns the start and end indices of the first pulse from vector(s)
% Usage: [indStart1, indEnd1, indStart2, indEnd2] = ...
%           find_pulse_endpoints (vectors)
% Outputs:
%       indStart1    - indices of pulse start (right before pulse start)
%                       specified as a positive integer (or NaN) vector
%       indEnd1      - indices of pulse end (right before pulse end)
%                       specified as a positive integer (or NaN) vector
%       indStart2    - indices of pulse start (right after pulse start)
%                       specified as a positive integer (or NaN) vector
%       indEnd2      - indices of pulse end (right after pulse end)
%                       specified as a positive integer (or NaN) vector
% Arguments:
%       vectors     - vectors containing a pulse
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%
% Requires:
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%       
% Used by:
%       cd/parse_pulse.m
%       cd/find_pulse_response_endpoints.m
%
% File History:
% 2018-07-25 BT - Adapted from find_initial_slopes.m
% 2018-08-10 Change the amplitude to take the value from pulseShifted
%                   rather than from pulse
% 2018-08-10 Now checks number of arguments
% 2018-09-17 Now returns empty indices if there is no pulse
% 2018-10-09 Improved documentation
% 2018-10-09 Now accepts multiple vectors an array or a cell array
% 2018-10-10 Added indStart2 and indEnd2
% 2018-12-15 Now uses force_column_cell.m
% 2018-12-15 Now returns NaN if there is no pulse

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

% Read from the Input Parser
parse(iP, vectors);

%% Preparation
% Force vectors as a cell array of numeric vectors
vectors = force_column_cell(vectors);

%% Do the job
[indStart1, indEnd1, indStart2, indEnd2] = ...
    cellfun(@find_pulse_endpoints_helper, vectors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idxStart1, idxEnd1, idxStart2, idxEnd2] = ...
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
    idxStart1 = find(vectorShifted(indHalf1) < amplitude * 0.25, 1, 'last');

    % Find the first point greater than 3/4 of the amplitude in the first half
    idxStart2 = find(vectorShifted(indHalf1) > amplitude * 0.75, 1);

    % Find the last point greater than 3/4 of the amplitude in the second half
    idxEnd1Rel = find(vectorShifted(indHalf2) > amplitude * 0.75, 1, 'last');

    % Find the first point less than 1/4 of the amplitude in the second half
    idxEnd2Rel = find(vectorShifted(indHalf2) < amplitude * 0.25, 1);
elseif signAmplitude == -1
    % Find the last point greater than 1/4 of the amplitude in the first half
    idxStart1 = find(vectorShifted(indHalf1) > amplitude * 0.25, 1, 'last');

    % Find the first point less than 3/4 of the amplitude in the first half
    idxStart2 = find(vectorShifted(indHalf1) < amplitude * 0.75, 1);

    % Find the last point less than 3/4 of the amplitude in the second half
    idxEnd1Rel = find(vectorShifted(indHalf2) < amplitude * 0.75, 1, 'last');

    % Find the first point greater than 1/4 of the amplitude in the second half
    idxEnd2Rel = find(vectorShifted(indHalf2) > amplitude * 0.25, 1);
else
    idxStart1 = [];
    idxStart2 = [];
    idxEnd1Rel = [];
    idxEnd2Rel = [];
    
end

% Either shift the indices to correspond to entire vector, 
%   or if not found, return as NaN
if isempty(idxStart1)
    idxStart1 = NaN;
end
if isempty(idxStart2)
    idxStart2 = NaN;
end
if isempty(idxEnd1Rel)
    idxEnd1 = NaN;
else
    idxEnd1 = idxEnd1Rel + idxAbsMax - 1;
end
if isempty(idxEnd2Rel)
    idxEnd2 = NaN;
else
    idxEnd2 = idxEnd2Rel + idxAbsMax - 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Find the pulse amplitude at that point
amplitude = pulse(idxAbsMax);

% Force vectors to be a column cell array
if iscell(vectors)
    vectors = vectors(:);
end

% Count the number of vectors
nVectors = count_vectors(vectors);

%% Do the job
indStart1 = zeros(nVectors, 1);
indEnd1 = zeros(nVectors, 1);
indStart2 = zeros(nVectors, 1);
indEnd2 = zeros(nVectors, 1);
if iscell(vectors)
    %parfor iVec = 1:nVectors
    for iVec = 1:nVectors
        [indStart1(iVec), indEnd1(iVec), ...
            indStart2(iVec), indEnd2(iVec)] = ...
            find_pulse_endpoints_helper(vectors{iVec});
    end
elseif isnumeric(vectors)
    for iVec = 1:nVectors
%    parfor iVec = 1:nVectors
        [indStart1(iVec), indEnd1(iVec), ...
            indStart2(iVec), indEnd2(iVec)] = ...
            find_pulse_endpoints_helper(vectors(:, iVec));
    end
else
    error('vectors is not the right type!');
end

idxStart1 = [];
idxStart2 = [];
idxEnd1 = [];
idxEnd2 = [];
return;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%