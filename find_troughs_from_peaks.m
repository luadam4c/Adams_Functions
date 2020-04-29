function [ampTroughs, indTroughs] = ...
                find_troughs_from_peaks (vec, indPeaks, varargin)
%% Finds troughs of a vector in between given peak indices
% Usage: [ampTroughs, indTroughs] = ...
%               find_troughs_from_peaks (vec, indPeaks, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load sunspot.dat
%       avSpots = sunspot(:, 2);
%       [~, indPeaks] = findpeaks(avSpots);
%       [a, i] = find_troughs_from_peaks(avSpots, indPeaks)
%
% Outputs:
%       ampTroughs  - TODO: Description of ampTroughs
%                   specified as a TODO
%       indTroughs  - TODO: Description of ampTroughs
%                   specified as a TODO
%
% Arguments:
%       vec         - TODO: Description of vec
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/compute_autocorrelogram.m
%       cd/parse_peaks.m

% File History:
% 2020-04-20 Moved from compute_autocorrelogram.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vec');
addRequired(iP, 'indPeaks');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, vec, indPeaks, varargin{:});
% param1 = iP.Results.param1;

%% Preparation
% Count the number of peaks
nPeaks = numel(indPeaks);

%% Do the job
if nPeaks < 2
    % No troughs
    ampTroughs = [];
    indTroughs = [];
else
    % Left peak indices
    indLeftPeak = indPeaks(1:(end-1));

    % Right peak indices
    indRightPeak = indPeaks(2:end);

    % Use the minimums in each interval
    [ampTroughs, indTroughsRel] = ...
        arrayfun(@(x, y) min(vec(x:y)), indLeftPeak, indRightPeak);

    % Compute the original indices
    indTroughs = arrayfun(@(x, y) x + y - 1, indTroughsRel, indLeftPeak);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%