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
%
% Requires:
%       None
%
% Used by:
%       scAAV/qupath_get_annotation_data.m

% File History:
% 2026-02-27 Modified from qupath_get_cm.m

%% Hard-coded parameters
validRoundModes = {'none', 'integer', 'index', 'valid'};

%% Default values for optional arguments
indicesDefault  = [];                   % default is set during preparation
roundModeDefault = 'none';              % default round mode

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

% Read from the Input Parser
parse(iP, signal, varargin{:});
indices = iP.Results.indices;
roundMode = validatestring(iP.Results.RoundMode, validRoundModes);

%% Preparation
% Make sure signal is a column vector
signal = signal(:);

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
rawCm = sum(indValid .* sigValid) / sum(sigValid);

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