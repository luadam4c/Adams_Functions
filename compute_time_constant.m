function [tauSamples, idxAtTau, valueAtTau] = ...
                compute_time_constant (vectors, varargin)
%% Computes the time constant of vector(s) with a single peak
% Usage: [tauSamples, idxAtTau, valueAtTau] = ...
%               compute_time_constant (vectors, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       tauSamples  - time constant in samples
%                   specified as a TODO
%       idxAtTau    - TODO
%                   specified as a TODO
%       valueAtTau  - TODO
%                   specified as a TODO
%
% Arguments:
%       vectors     - vectors with only a single peak
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for compute_peak_decay()
%
% Requires:
%       cd/compute_peak_decay.m
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/m3ha_plot_gabab_ipsc.m

% File History:
% 2020-01-03 Created by Adam Lu
% TODO: Case for cell array inputs
% TODO: Add 'PeakDirection' as an optional argument

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'vectors', ...
    @(x) assert(isnumeric(x) || iscell(x), ...
                ['vectors must be either a numeric array', ...
                    'or a cell array!']));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, vectors, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the compute_peak_decay() function
otherArguments = iP.Unmatched;

%% Do the job
% Find the maximum values for each vector
[~, indPeaks] = max(vectors, [], 1);

% Compute the decay time constant in samples
[tauSamples, idxAtTau, valueAtTau] = ...
    compute_peak_decay(vectors, indPeaks, otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%