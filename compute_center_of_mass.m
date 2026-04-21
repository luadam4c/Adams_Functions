function [idxCm, indValid, sigValid] = compute_center_of_mass(signal, varargin)
%% Calculates the center of mass index for a 1D distribution
% Usage: [idxCm, indValid, sigValid] = compute_center_of_mass(signal, indices, varargin)
% Explanation:
%       Determines the weighted average position (center of mass) of a signal.
%
% Example(s):
%       [idxCm, indValid, sigValid] = compute_center_of_mass([0, 1, 2, 1, 0], 'RoundMode', 'integer')
%       [idxCm, indValid, sigValid] = compute_center_of_mass([0, NaN, 2, NaN, 2], 'RoundMode', 'integer')
%       [idxCm, indValid, sigValid] = compute_center_of_mass([0, NaN, 2, NaN, 2], 'RoundMode', 'valid')
%       [idxCm, indValid, sigValid] = compute_center_of_mass([1, 1, 1], [], 'Weights', [0.5, 1, 0.5])
%
% Outputs:
%       idxCm       - The calculated center of mass index
%                   specified as a numeric scalar
%       indValid    - vector of indices corresponding to valid signal values
%                   specified as a numeric vector
%       sigValid    - vector of valid signal values
%                   specified as a numeric vector
%
% Arguments:
%       signal      - vector of signal values corresponding to each index
%                   must be a numeric vector
%       indices     - (opt) vector of indices representing spatial positions
%                   must be a numeric vector
%                   default == (1:length(signal))'
%       varargin    - 'RoundMode': how to round the calculated center of mass
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'none'      - do not round the center of mass
%                       'integer'   - round to the nearest integer
%                       'index'     - round to the nearest index within the indices vector
%                       'valid'     - round to the nearest index that is valid
%                   default == 'none'
%                   - 'Weights': optional vector of weights to dot multiply with signal
%                   must be a numeric vector of the same length as signal
%                   default == []
%
% Requires:
%       None
%
% Used by:
%       scAAV/qupath_get_annotation_data.m

% File History:
% 2026-02-27 Modified from qupath_get_cm.m
% 2026-03-12 Updated to use absolute value of signal values for weighting
% 2026-03-17 Added 'Weights' as an optional argument

%% Hard-coded parameters
validRoundModes = {'none', 'integer', 'index', 'valid'};

%% Default values for optional arguments
indicesDefault  = [];                   % default is set during preparation
roundModeDefault = 'none';              % default round mode
weightsDefault = [];                    % default is no weights

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error('Not enough input arguments, type ''help %s'' for usage', mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'signal', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'indices', indicesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'RoundMode', roundModeDefault, ...
    @(x) any(validatestring(x, validRoundModes)));

addParameter(iP, 'Weights', weightsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Read from the Input Parser
parse(iP, signal, varargin{:});
indices = iP.Results.indices;
roundMode = validatestring(iP.Results.RoundMode, validRoundModes);
weights = iP.Results.Weights;

%% Preparation
% Make sure signal is a column vector
signal = signal(:);

% Apply weights to the signal if provided
if ~isempty(weights)
    weights = weights(:);
    if length(weights) ~= length(signal)
        error('Weights vector must have the same length as the signal vector.');
    end
    signal = signal .* weights;
end

% Set default indices if empty or force as column vector
if isempty(indices)
    indices = (1:length(signal))';
else
    indices = indices(:);
end

% Ensure all values are valid
isValid = ~isnan(signal);
indValid = indices(isValid);
sigValid = signal(isValid);

%% Do the job
% Compute the raw weighted average position (center of mass)
% using the absolute value of the signal as weights
rawCm = sum(indValid .* abs(sigValid)) / sum(abs(sigValid));

% Apply the requested rounding mode
switch roundMode
    case 'none'
        idxCm = rawCm;
    case 'integer'
        idxCm = round(rawCm);
    case 'index'
        [~, minIdx] = min(abs(indices - rawCm));
        idxCm = indices(minIdx);
    case 'valid'
        [~, minIdx] = min(abs(indValid - rawCm));
        idxCm = indValid(minIdx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%