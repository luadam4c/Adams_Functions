function statsTable = m3ha_compute_ltsburst_statistics (varargin)
%% Computes LTS and burst statistics for some indices (indOfInterest) in ltsDelays, spikesPerPeak, burstDelays & spikesPerBurst
% Usage: statsTable = m3ha_compute_ltsburst_statistics (varargin)
%
% Requires:
%       cd/m3ha_load_sweep_info.m
%       cd/m3ha_select_sweeps_to_fit.m
%
% Used by:
%       cd/m3ha_compute_and_plot_statistics.m
%
% Outputs:
%       statsTable  - a table containing measures as row names
%                       and the following variables:
%                           measureTitle
%                           measureStr
%                           allValues
%                           nValues
%                           meanValue
%                           stdValue
%                           errValue
%                           upper95Value
%                           lower95Value
%                   specified as a table
%
% Arguments:
%       varargin    - 'DataMode': data mode
%                   must be one of:
%                       0 - all data
%                       1 - all of g incr = 100%, 200%, 400%
%                       2 - all of g incr = 100%, 200%, 400% 
%                               but exclude cell-pharm-g_incr sets 
%                               containing problematic sweeps
%                   default == 0
%                   - 'SwpInfo': a table of sweep info, with each row named by 
%                               the matfile base containing the raw data
%                   must a 2D table with row names being file bases
%                       and with the fields:
%                       cellidrow       - cell ID
%                       ltspeaktimes    - LTS peak delays
%                       spikesperpeak   - action potentials per LTS
%                       bursttimes      - burst delays
%                       spikesperburst  - action potentials per burst
%                   default == m3ha_load_sweep_info
%                   - 'PharmCondition': pharm condition
%                   must be a positive integer scalar
%                   default == not provided
%                   - 'GIncrCondition': gIncr condition
%                   must be a scalar
%                   default == not provided
%                   - 'VHoldCondition': vHold condition
%                   must be a scalar
%                   default == not provided

% File History:
% 2016-08-19 Created
% 2016-08-29 Last Modified
% 2017-01-25 - Corrected errors to reflect t-confidence intervals 
%               (from the Gosset's t distribution)
% 2019-11-26 Improved code structure
% 2019-11-27 Now computes means of all measures from the data points
%               for each cell (rather than for each sweep)
% 2019-11-27 Added 'PharmCondition', 'GIncrCondition' & 'VHoldCondition'
%               as optional arguments

%% Hard-coded parameters
% Items to compute
measureTitle = {'LTS onset time (ms)'; 'LTS time jitter (ms)'; ...
                'LTS probability'; 'Spikes per LTS'; ...
                'Burst onset time (ms)'; 'Burst time jitter (ms)'; ...
                'Burst probability'; 'Spikes per burst'};
%   Note: measureStr must match the variables defined in this file!
measureStr = {'ltsOnsetTime'; 'ltsTimeJitter'; ...
                    'ltsProbability'; 'spikesPerLts'; ...
                    'burstOnsetTime'; 'burstTimeJitter'; ...
                    'burstProbability'; 'spikesPerBurst'};
cellIdStr = 'cellidrow';
ltsDelayStr = 'ltspeaktimes';
spikesPerPeakStr = 'spikesperpeak';
burstDelayStr = 'bursttimes';
spikesPerBurstStr = 'spikesperburst';

%% Default values for optional arguments
dataModeDefault = 0;
swpInfoDefault = table.empty;
pharmConditionDefault = [];
gIncrConditionDefault = [];
vHoldConditionDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fixed parameters used in the experiments
cc = 1:1:49;                % Possible cell ID #s

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'DataMode', dataModeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
addParameter(iP, 'SwpInfo', swpInfoDefault, ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addParameter(iP, 'PharmCondition', pharmConditionDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
addParameter(iP, 'GIncrCondition', gIncrConditionDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'VHoldCondition', vHoldConditionDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Read from the Input Parser
parse(iP, varargin{:});
dataMode = iP.Results.DataMode;
swpInfo = iP.Results.SwpInfo;
pharmCondition = iP.Results.PharmCondition;
gIncrCondition = iP.Results.GIncrCondition;
vHoldCondition = iP.Results.VHoldCondition;

% Count the items to compute
nMeasures = numel(measureTitle);

%% Preparation
% Read the sweep info data
if isempty(swpInfo)
    swpInfo = m3ha_load_sweep_info('FileName', dataPath);
end

%% Select sweeps
% Select the sweeps based on data mode
swpInfo = m3ha_select_sweeps_to_fit('SwpInfo', swpInfo, ...
                                    'DataMode', dataMode, ...
                                    'PharmCondition', pharmCondition, ...
                                    'GIncrCondition', gIncrCondition, ...
                                    'VHoldCondition', vHoldCondition);

% Extract whether to use the sweep
toUse = swpInfo.toUse;

% Restrict swpInfo to those sweeps
swpInfoToUse = swpInfo(toUse, :);

%% Determine all possible cells
% Extract the cell IDs
cellIdRow = swpInfoToUse.(cellIdStr);

% Find unique cell IDs
uniqueCellIds = unique(cellIdRow);

% Count the number of cells
nCells = numel(uniqueCellIds);

%% Extract measures
% Extract the measures
ltsDelays = swpInfoToUse.(ltsDelayStr);
spikesPerPeak = swpInfoToUse.(spikesPerPeakStr);
burstDelays = swpInfoToUse.(burstDelayStr);
spikesPerBurst = swpInfoToUse.(spikesPerBurstStr);

% Determine whether each sweep has an LTS or has a burst
hasLts = ltsDelays > 0;
hasBurst = burstDelays > 0;

%% Compute LTS & burst measures for each cell
% Compute the LTS probability for each cell
ltsProbability = arrayfun(@(x) sum(cellIdRow == x & hasLts) / ...
                                sum(cellIdRow == x), uniqueCellIds);

% Compute the burst probability for each cell
burstProbability = arrayfun(@(x) sum(cellIdRow == x & hasBurst) / ...
                                sum(cellIdRow == x), uniqueCellIds);

% Compute means of LTS and burst properties for each cell
[ltsOnsetTime, spikesPerLts, burstOnsetTime, spikesPerBurst] = ...
    argfun(@(x) arrayfun(@(y) nanmean(x(cellIdRow == y)), uniqueCellIds), ...
            ltsDelays, spikesPerPeak, burstDelays, spikesPerBurst);

% Compute standard deviations of LTS and burst properties for each cell
[ltsTimeJitter, burstTimeJitter] = ...
    argfun(@(x) arrayfun(@(y) nanstd(x(cellIdRow == y)), uniqueCellIds), ...
            ltsDelays, burstDelays);

%% Calculate overall LTS & burst statistics
% Collect all values for each measure
%   Note: Each element of the cell array contains a nCells by 1 numeric vector
%           that may contain NaNs
allValues = cellfun(@(x) eval(x), measureStr);

% Compute effective n's (number of values that are not NaNs)
nValues = cellfun(@(x) sum(~isnan(x)), allValues);

% Compute overall statistics
[meanValue, stdValue, errValue, upper95Value, lower95Value] = ...
    argfun(@(x) cellfun(@(y) compute_stats(y, x, 'IgnoreNan', true), ...
                        allValues), ...
            'mean', 'std', 'err', 'upper95', 'lower95');

%% Output results
statsTable = table(measureTitle, measureStr, allValues, nValues, ...
                    meanValue, stdValue, errValue, upper95Value, ...
                    lower95Value, 'RowNames', measureStr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
