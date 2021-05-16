function [lagBest, corrBest, lagAll, corrAll] = ...
                compute_best_lag (signal1, signal2, varargin)
%% Computes the lag with highest correlation between two signals
% Usage: [lagBest, corrBest, lagAll, corrAll] = ...
%               compute_best_lag (signal1, signal2, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [lagBest, corrBest, lagAll, corrAll] = compute_best_lag([0, 0, 1], [1, 0, 0])
%       [lagBest, corrBest, lagAll, corrAll] = compute_best_lag([1, 0, 0], [0, 0, 1])
%
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
%
% Arguments:
%       signal1     - first signal
%                   must be an array accepted by xcorr
%       signal2     - second signal
%                   must be an array accepted by xcorr
%       varargin    - 'MaxLag': maximum absolute lag in samples
%                   must be empty or NaN or Inf or a positive scalar
%                   default == Inf
%                   - 'ScaleOption': scaling option (see documenttion for xcorr)
%                   must be a string scalar or character vector
%                   default == 'normalized'
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/crosscorr_profile.m

% File History:
% 2021-05-15 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
maxLagDefault = Inf;
scaleOptionDefault = 'normalized';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'signal1');
addRequired(iP, 'signal2');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MaxLag', maxLagDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'ScaleOption', scaleOptionDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, signal1, signal2, varargin{:});
maxLag = iP.Results.MaxLag;
scaleOption = iP.Results.ScaleOption;

%% Do the job
% Compute the cross-correlogram between two signals
if ~isempty(maxLag) && ~isnan(maxLag) && ~isinf(maxLag)
    [corrAll, lagAll] = xcorr(signal1, signal2, maxLag, scaleOption);
else
    [corrAll, lagAll] = xcorr(signal1, signal2, scaleOption);
end

% Find the lag with highest correlation with the largest correlation
[corrBest, I] = max(abs(corrAll));
lagBest = lagAll(I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%