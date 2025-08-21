function [trialNos, nTrialsEachGroup, nTrialsTotal] = compute_trial_numbers (eventTimeArray, varargin)
%% Computes trial numbers for a cell array of event time arrays
% Usage: [trialNos, nTrialsEachGroup, nTrialsTotal] = compute_trial_numbers (eventTimeArray, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       trialNos     - TODO: Description of trialNos
%                   specified as a TODO
%
% Arguments:
%       eventTimeArray     - TODO: Description of eventTimeArray
%                   must be a numeric array or a cell array of numeric arrays
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/count_vectors.m
%
% Used by:
%       cd/plot_raster.m

% File History:
% 2025-08-21 Extracted from plot_raster.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'eventTimeArray');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, eventTimeArray, varargin{:});
param1 = iP.Results.param1;

%% Do the job
% Get the number of vectors in each event time array
nTrialsEachGroup = cellfun(@count_vectors, eventTimeArray);

% If there is an empty vector, consider it a vector of size 1 for placement
nTrialsEachGroup(nTrialsEachGroup == 0) = 1;

% Get the total number of vectors
nTrialsTotal = sum(nTrialsEachGroup);

% Assign trial numbers grouped by each event time array
trialNos = cell(size(eventTimeArray));
ct = 0; % counts number of trials assigned
for iArray = 1:numel(eventTimeArray)
    % Get the number of vectors in this array
    nVectorThis = nTrialsEachGroup(iArray);

    % Assign the trial numbers in ascending order
    trialNosThis = ct + transpose(1:nVectorThis);

    % Store the trial numbers for this event time array
    trialNos{iArray} = trialNosThis;

    % Update number of trials assigned
    ct = ct + nVectorThis;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%