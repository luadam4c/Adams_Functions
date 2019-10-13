function [simData, timeVecs] = m3ha_extract_and_match_time_points(realData, simDataOrig)
%% Interpolate simulated data to match the time points of real data
% Usage: [simData, timeVecs] = m3ha_extract_and_match_time_points(realData, simDataOrig)
% Explanation:
%       This is necessary whenever CVODE (variable time step method) 
%       is applied in NEURON
% Example(s):
%       TODO
% Outputs:
%       simData     - matched simulated data
%                   specified as a cell array of numeric arrays
%       timeVecs    - matched time vectors
%                   specified as a cell array of numeric vectors
% Arguments:    
%       realData    - recorded data
%                   must be a cell array of numeric arrays
%       simDataOrig - unmatched simulated data
%                   must be a cell array of numeric arrays
%
% Requires:
%       cd/match_time_points.m
%
% Used by:    
%       ~/m3ha/optimizer4gabab/run_neuron_once_4compgabab.m

% File History:
% 2018-10-23 Adapted from code in run_neuron_once_4compgabab.m
% 

%% Hard-coded parameters

%% Default values for optional arguments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'realData', ...                 % recorded data
    @(x) assert(iscell(x) && all(cellfun(@isnumeric, x)), ...
        'data must be a cell array of numeric arrays!'));
addRequired(iP, 'simDataOrig', ...              % unmatched simulated data
    @(x) assert(iscell(x) && all(cellfun(@isnumeric, x)), ...
        'data must be a cell array of numeric arrays!'));

% Read from the Input Parser
parse(iP, realData);

%% Do the job
% Get all time vectors from the real data
timeVecs = cellfun(@(x) x(:, 1), realData, 'UniformOutput', false);

% Interpolate to get simulated data
simData = cellfun(@(x, y) match_time_points(x, y), simDataOrig, timeVecs, ...
                    'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

timeVecs = cell(nSweeps, 1);          % time vectors from the real data
simData = cell(nSweeps, 1);           % interpolated simulation data
parfor iSwp = 1:nSweeps
    % Get the time vector from the real data
    timeVecs{iSwp} = realData{iSwp}(:, 1);

    % Interpolate to get new simulated data
    simData{iSwp} = match_time_points(simDataOrig{iSwp}, timeVecs{iSwp});
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
