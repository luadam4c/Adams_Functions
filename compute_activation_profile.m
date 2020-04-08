function [percentActivated, timeBins] = compute_activation_profile (cellIds, spikeTimes, varargin)
%% Computes the percent of activated cells in the network over time
% Usage: [percentActivated, timeBins] = compute_activation_profile (cellIds, spikeTimes, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [percentActivated, timeBins] = compute_activation_profile ([0, 1, 0, 2, 3], [3000, 3000, 4500, 6000, 8000]);
%
% Outputs:
%       activationProfile   - TODO: Description of activationProfile
%                           specified as a numeric array
%       timeBins   - TODO: Description of activationProfile
%                           specified as a numeric array
%
% Arguments:
%       cellIds     - TODO: Description of cellIds
%                   must be a TODO
%       spikeTimes  - TODO: Description of cellIds
%                   must be a TODO
%       varargin    - 'TimeBins': TODO: Description of TimeBins
%                   must be a TODO
%                   default == TODO
%                   - 'NCells': number of cells in the network
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/apply_over_cells.m
%       cd/argfun.m
%       cd/compute_grouped_histcounts.m
%       cd/create_error_for_nargin.m
%       cd/create_time_vectors.m
%       cd/find_first_match.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
%       cd/match_format_vector_sets.m
%       cd/unique_custom.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2020-04-08 Created by Adam Lu
% 

%% Hard-coded parameters
binWidthMs = 500;
stimEndMs = 30000;

%% Default values for optional arguments
timeBinsDefault = [];
nCellsDefault = [];

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
addRequired(iP, 'cellIds');
addRequired(iP, 'spikeTimes');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TimeBins', timeBinsDefault);
addParameter(iP, 'NCells', nCellsDefault);

% Read from the Input Parser
parse(iP, cellIds, spikeTimes, varargin{:});
timeBins = iP.Results.TimeBins;
nCells = iP.Results.NCells;

%% Preparation
% Force as column cell arrays of numeric vectors
[cellIds, spikeTimes] = match_format_vector_sets(cellIds, spikeTimes, ...
                                                'ForceCellOutput', true);

% Construct time bins if not provided
if isempty(timeBins)
    nBins = floor(stimEndMs/binWidthMs);
    timeBins = create_time_vectors(nBins, 'SamplingIntervalMs', binWidthMs, ...
                                    'TimeUnits', 'ms');
end

% Decide on number of cells if not provided
if isempty(nCells)
    uniqueIds = apply_over_cells(@(x, y) unique_custom(union(x, y), ...
                                        'IgnoreNan', true), cellIds);
    nCells = numel(uniqueIds);
end

%% Do the job
[percentActivatedCell, timeBinsCell] = ...
    cellfun(@(c, s) compute_one_activation_profile(c, s, timeBins, nCells), ...
            cellIds, spikeTimes, 'UniformOutput', false);

% Force as a matrix
[percentActivated, timeBins] = ...
    argfun(@force_matrix, percentActivatedCell, timeBinsCell);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [percentActivated, timeBins] = ...
                compute_one_activation_profile (cellIds, spikeTimes, ...
                                                timeBins, nCells)

% Force as a column vector
[cellIds, spikeTimes] = argfun(@force_column_vector, cellIds, spikeTimes);

% Find the first spike times for each cell ID
uniqueIds = unique_custom(cellIds, 'IgnoreNan', true);
indFirstMatch = find_first_match(uniqueIds, cellIds);
firstSpikeTimes = spikeTimes(indFirstMatch);

% Add zero to timeBins
if timeBins(1) == 0
    timeBinsWithZero = timeBins;
    timeBins = timeBins(2:end);
else
    timeBinsWithZero = [0; timeBins];
    timeBins = timeBins;
end

% Count the number of cells that had its first spike for each time bin
nCellsFirstSpike = compute_grouped_histcounts(firstSpikeTimes, ...
                                            'Edges', timeBinsWithZero);

% Count the cumulative number of cells activated for each time bin
nCellsActivated = cumsum(nCellsFirstSpike);

% Compute the percent of cells activated over time
percentActivated = nCellsActivated ./ nCells;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%