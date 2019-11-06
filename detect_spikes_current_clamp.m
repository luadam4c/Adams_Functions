function [spikesParams, spikesData] = detect_spikes_current_clamp (vVecs, varargin)
%% Detects spikes from a current clamp recording
% Usage: [spikesParams, spikesData] = detect_spikes_current_clamp (vVecs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       spikesParams- Used and detected parameters, with fields:
%                       minAmpBefore
%                       minAmpAfter
%                   specified as a scalar structure
%       spikesData  - Detected spikes data, with fields:
%                       isSpike
%                       idxSpikes
%                   specified as a scalar structure
%
% Arguments:
%       vVecs       - voltage vectors
%                   must be a numeric array or a cell array of numeric vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
%
% Used by:
%       cd/plot_FI.m

% File History:
% 2019-11-05 Moved from plot_FI.m
% TODO: Make MinAmpBefore & MinAmpAfter optional arguments
% TODO: Make SiMs, tVec, and PlotSpikeDetection optional arguments

%% Hard-coded parameters
% TODO: Make optional parameters
minAmpBefore = 10;     % in mV
minAmpAfter = 5;       % in mV

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

% Add required inputs to the Input Parser
addRequired(iP, 'vVecs', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vVecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, vVecs, varargin{:});
param1 = iP.Results.param1;

%% Preparation
% Determine whether the first argument is a matrix
if isnumeric(vVecs)
    wasNumeric = true;
else
    wasNumeric = false;
end

% Force as a cell array of column vectors
vVecs = force_column_vector(vVecs, 'IgnoreNonVector', false, ...
                            'ForceCellOutput', true);

%% Do the job
% Determine whether each sample point is a spike for each sweep
isSpike = cellfun(@(x) detect_spikes_one_sweep(x, ...
                                minAmpBefore, minAmpAfter), ...
                    vVecs, 'UniformOutput', false);

% Find the spike indices
idxSpikes = cellfun(@(x) find(x), isSpike, 'UniformOutput', false);

% Force as a numeric array if a numeric array is passed in
if wasNumeric
    % Put vectors together into a matrix
    isSpike = force_matrix(isSpike);

    % Extract from the first cell if only one vector is provided
    if numel(idxSpikes) == 1
        idxSpikes = idxSpikes{1};
    end
end

%% Output results
% Output scalars
spikesParams.minAmpBefore = minAmpBefore;
spikesParams.minAmpAfter = minAmpAfter;

% Output vectors
spikesData.isSpike = isSpike;
spikesData.idxSpikes = idxSpikes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isSpike = detect_spikes_one_sweep(vVec, minAmpBefore, minAmpAfter)
%% Determines whether each sample point is a spike
%   Note: Criteria for a spike: 
%   (1) Must be a local maximum 10 mV higher than the previous local minimum
%   (2) Must be 5 mV higher than the minimum value between the spike 
%           and the following spike

%% Preparation
% Count the number of samples
nSamples = numel(vVec);

% Initialize array
isSpike = false(nSamples, 1);

%% Do the job
% No spikes if nSamples is less than 3
if nSamples <= 2
    return
end

% Determine whether each data point is a "local maximum" or a "local minimum" 
%   in the general sense
%   Note: This is probably faster than findpeaks(). Also different
%           in that <= or >= is used on the right side
isLocalMaximum = [false; diff(vVec(1:end-1)) > 0 & diff(vVec(2:end)) <= 0; ...
                    false];
isLocalMinimum = [false; diff(vVec(1:end-1)) < 0 & diff(vVec(2:end)) >= 0; ...
                    false];

% Finds all peaks that satisfy criterion 1
for idx = 2:nSamples-1
    % Only check local maxima
    if isLocalMaximum(idx)
        % Find the index of the previous local minimum
        idxPrevLocalMinimum = idx - 1;
        while ~isLocalMinimum(idxPrevLocalMinimum) && idxPrevLocalMinimum > 1
            idxPrevLocalMinimum = idxPrevLocalMinimum - 1;
        end

        % Make sure a previous local minimum was found
        if idxPrevLocalMinimum > 1
            % Determine whether criterion 1 is satisfied
            isSpike(idx) = vVec(idx) - vVec(idxPrevLocalMinimum) > minAmpBefore;
        end
    end
end

% Finds all peaks that satisfy criterion 2 by filtering out those 
%   that pass criterion 1 in reverse order
for idx = nSamples-1:-1:2
    % Only check peaks that satisfy criterion 1
    if isSpike(idx)
        % Find the index of the following spike
        idxNextSpike = idx + 1;
        while ~isSpike(idxNextSpike) && idxNextSpike < nSamples
            idxNextSpike = idxNextSpike + 1;
        end

        % Determine whether criterion 2 is satisfied
        isSpike(idx) = vVec(idx) - min(vVec(idx:idxNextSpike)) > minAmpAfter;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%