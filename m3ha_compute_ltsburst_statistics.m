function statsTable = m3ha_compute_ltsburst_statistics (varargin)
%% Computes LTS and burst statistics for some indices (indOfInterest) in ltsDelays, spikesPerPeak, burstDelays & spikesPerBurst
% Usage: statsTable = m3ha_compute_ltsburst_statistics (varargin)
%
% Explanation:
%       TODO
%
% Example(s):
%       statsTable = m3ha_compute_ltsburst_statistics
%       statsTable = m3ha_compute_ltsburst_statistics('DataMode', 2)
%       statsTable = m3ha_compute_ltsburst_statistics('PharmConditions', 1)
%       statsTable = m3ha_compute_ltsburst_statistics('PharmConditions', 1:4)
%       statsTable = m3ha_compute_ltsburst_statistics('PharmConditions', 1:4, 'GIncrCondition', [100; 200; 400])
%
% Requires:
%       cd/argfun.m
%       cd/compute_stats.m
%       cd/force_column_cell.m
%       cd/force_column_vector.m
%       cd/match_row_count.m
%       cd/m3ha_load_sweep_info.m
%       cd/m3ha_select_sweeps.m
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
%                           stderrValue
%                           errValue
%                           upper95Value
%                           lower95Value
%                   specified as a table
%
% Arguments:
%       varargin    - 'SwpInfo': a table of sweep info, with each row named by 
%                               the matfile base containing the raw data
%                   must a 2D table with row names being file bases
%                       and with the fields:
%                       cellidrow       - cell ID
%                       ltspeaktimes    - LTS peak delays
%                       spikesperpeak   - action potentials per LTS
%                       bursttimes      - burst delays
%                       spikesperburst  - action potentials per burst
%                   default == m3ha_load_sweep_info
%                   - 'DataMode': data mode
%                   must be one of:
%                       0 - all data
%                       1 - all of g incr = 100%, 200%, 400%
%                       2 - all of g incr = 100%, 200%, 400% 
%                               but exclude cell-pharm-g_incr sets 
%                               containing problematic sweeps
%                   default == 0
%                   - 'PharmConditions': pharmacological condition(s)
%                                           to restrict to
%                   must be empty or some of:
%                       1 - control
%                       2 - GAT1 blockade
%                       3 - GAT3 blockade
%                       4 - dual blockade
%                       or a cell array of them (will become 1st dimension)
%                   default == no restrictions
%                   - 'GIncrCondition': conductance amplitude condition(s) (%)
%                                           to restrict to
%                   must be empty or some of: 25, 50, 100, 200, 400, 800
%                       or a cell array of them (will become 2nd dimension)
%                   default == no restrictions
%                   - 'VHoldConditions': holding potential condition(s) (mV)
%                                           to restrict to
%                   must be empty or some of: -60, -65, -70
%                       or a cell array of them (will become 3rd dimension)
%                   default == no restrictions

% File History:
% 2016-08-19 Created
% 2016-08-29 Last Modified
% 2017-01-25 - Corrected errors to reflect t-confidence intervals 
%               (from the Gosset's t distribution)
% 2019-11-26 Improved code structure
% 2019-11-27 Now computes means of all measures from the data points
%               for each cell (rather than for each sweep)
% 2019-11-27 Added 'PharmConditions', 'GIncrConditions' & 'VHoldConditions'
%               as optional arguments
% TODO: Other LTS features that might be of interest
%    ltspeakval = swpInfo.ltspeakval;
%    maxslopeval = swpInfo.maxslopeval;
%    ltspeak2ndder = swpInfo.ltspeak2ndder;
%    ltspeakprom = swpInfo.ltspeakprom;
%    ltspeakwidth = swpInfo.ltspeakwidth;
% TODO: Other burst features that might be of interest
%    maxspikeamp = swpInfo.maxspikeamp;
%    minspikeamp = swpInfo.minspikeamp;
%    spikefrequency = swpInfo.spikefrequency;
%    spikeadaptation = swpInfo.spikeadaptation;

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
ltsDelayStr = 'ltspeaktime';
spikesPerPeakStr = 'spikesperpeak';
burstDelayStr = 'bursttime';
spikesPerBurstStr = 'spikesperburst';

%% Default values for optional arguments
swpInfoDefault = table.empty;
dataModeDefault = 0;
pharmConditionsDefault = [];
gIncrConditionsDefault = [];
vHoldConditionsDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SwpInfo', swpInfoDefault, ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addParameter(iP, 'DataMode', dataModeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
addParameter(iP, 'PharmConditions', pharmConditionsDefault, ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'2d'}));
addParameter(iP, 'GIncrConditions', gIncrConditionsDefault, ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'2d'}));
addParameter(iP, 'VHoldConditions', vHoldConditionsDefault, ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'2d'}));

% Read from the Input Parser
parse(iP, varargin{:});
swpInfo = iP.Results.SwpInfo;
dataMode = iP.Results.DataMode;
pharmConditionsUser = iP.Results.PharmConditions;
gIncrConditionsUser = iP.Results.GIncrConditions;
vHoldConditionsUser = iP.Results.VHoldConditions;

% Count the items to compute
nMeasures = numel(measureTitle);

%% Preparation
% Read the sweep info data
if isempty(swpInfo)
    swpInfo = m3ha_load_sweep_info;
end

%% Create all conditions
if iscell(pharmCondition) || iscell(gIncrCondition) || ...
        iscell(vHoldCondition)
    % Force as cell arrays of column vectors
    [pharmConditionsUser, gIncrConditionsUser, vHoldConditionsUser] = ...
        argfun(@(x) force_column_vector(x, 'ToLinearize', true, ...
                                        'ForceCellOutput', true), ...
                pharmConditionsUser, gIncrConditionsUser, vHoldConditionsUser);

    % Force as column cell arrays
    [pharmConditionsUser, gIncrConditionsUser, vHoldConditionsUser] = ...
        argfun(@force_column_cell, ...
                pharmConditionsUser, gIncrConditionsUser, vHoldConditionsUser);

    % Count the number of values for each type of condition
    nPharm = numel(pharmConditionsUser);
    nGIncr = numel(gIncrConditionsUser);
    nVHold = numel(vHoldConditionsUser);

    % Create condition labels for each condition
    conditionLabel = cell(nPharm, nGIncr, nVHold);
    for iPharm = 1:nPharm

    % Create all set of conditions
    pharmCondition = arrayfun(nPharm, nGIncr, nVHold, ...
                                'UniformOutput', false);
    gIncrCondition = arrayfun(nPharm, nGIncr, nVHold, ...
                                'UniformOutput', false);
    vHoldCondition = arrayfun(nPharm, nGIncr, nVHold, ...
                                'UniformOutput', false);

    % Set flag
    manyConditionsFlag = true;
else
    % There is only one set of conditions
    pharmCondition = pharmConditionsUser;
    gIncrCondition = gIncrConditionsUser;
    vHoldCondition = vHoldConditionsUser;

    % Set flag
    manyConditionsFlag = false;
end

%% Compute statistics
if manyConditionsFlag
    % Compute statistics for each set of conditions
    %   Note: This will return a nPharm x nGIncr x nVhold cell array 
    %           where each element is either a nMeasures x 1 column cell vector
    %            or a nMeasures x 1 column numeric vector
    [allValues, nValues, meanValue, stdValue, ...
            stderrValue, errValue, upper95Value, lower95Value] = ...
        cellfun(@(x, y, z) ...
                m3ha_compute_ltsburst_statistics_helper(swpInfo, dataMode, ...
                x, y, z, cellIdStr, ltsDelayStr, spikesPerPeakStr, ...
                burstDelayStr, spikesPerBurstStr), ...
                pharmCondition, gIncrCondition, vHoldCondition, ...
                'UniformOutput', false);

    % Reorganize cell arrays of cell arrays 
    %   so that measures are grouped together
    %   Note: This will return a nMeasures x 1 cell array,
    %           where each element is a nPharm x nGIncr x nVhold cell array
    allValues = extract_columns(allValues, transpose(1:nMeasures), ...
                            'TreatCnvAsColumns', true, 'OutputMode', 'single');

    % Reorganize cell arrays of numeric vectors 
    %   so that measures are grouped together
    %   Note: This will return a nMeasures x 1 cell array,
    %           where each element is a nPharm x nGIncr x nVhold numeric array
    [nValues, meanValue, stdValue, ...
            stderrValue, errValue, upper95Value, lower95Value] = ...
        argfun(@(x) ...
                arrayfun(@(y) extract_elements(x, 'specific', 'Index', y), ...
                        transpose(1:nMeasures), 'UniformOutput', true), ...
            nValues, meanValue, stdValue, ...
            stderrValue, errValue, upper95Value, lower95Value);
else
    % Compute statistics for this set of conditions
    %   Note: this will return either a nMeasures x 1 column cell vector
    %            or a nMeasures x 1 column numeric vector
    [allValues, nValues, meanValue, stdValue, ...
            stderrValue, errValue, upper95Value, lower95Value] = ...
        m3ha_compute_ltsburst_statistics_helper(swpInfo, dataMode, ...
                pharmCondition, gIncrCondition, vHoldCondition, ...
                cellIdStr, ltsDelayStr, spikesPerPeakStr, ...
                burstDelayStr, spikesPerBurstStr);
end

%% Output results
% Match the number of rows
[dataMode, pharmCondition, gIncrCondition, vHoldCondition] = ...
    argfun(@(x) match_row_count(x, nMeasures), ...
            dataMode, pharmCondition, gIncrCondition, vHoldCondition);

% Create a statistics table
statsTable = table(measureTitle, measureStr, dataMode, ...
                    pharmCondition, gIncrCondition, vHoldCondition, ...
                    allValues, nValues, meanValue, stdValue, ...
                    stderrValue, errValue, upper95Value, lower95Value, ...
                    'RowNames', measureStr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [allValues, nValues, meanValue, stdValue, ...
                    stderrValue, errValue, upper95Value, lower95Value] = ...
                m3ha_compute_ltsburst_statistics_helper(swpInfo, dataMode, ...
                            pharmConditions, gIncrConditions, vHoldConditions, ...
                            cellIdStr, ltsDelayStr, spikesPerPeakStr, ...
                            burstDelayStr, spikesPerBurstStr)
%% Computes LTS and burst statistics for one condition

%% Select sweeps
% Select the sweeps based on data mode
swpInfo = m3ha_select_sweeps('SwpInfo', swpInfo, 'Verbose', false, ...
                                'DataMode', dataMode, ...
                                'PharmConditions', pharmConditions, ...
                                'GIncrConditions', gIncrConditions, ...
                                'VHoldConditions', vHoldConditions);

% Extract whether to use the sweep
toUse = swpInfo.toUse;

% Restrict swpInfo to those sweeps
swpInfoToUse = swpInfo(toUse, :);

%% Determine all possible cells
% Extract the cell IDs
cellIdRow = swpInfoToUse.(cellIdStr);

% Find unique cell IDs
uniqueCellIds = unique(cellIdRow);

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
%   Note: cellfun won't work here because of namespace differences
allValues = cell(nMeasures, 1);
for iMeasure = 1:nMeasures
    allValues{iMeasure} = eval(measureStr{iMeasure});
end

% Compute effective n's (number of values that are not NaNs)
nValues = cellfun(@(x) sum(~isnan(x)), allValues);

% Compute overall statistics
[meanValue, stdValue, stderrValue, errValue, upper95Value, lower95Value] = ...
    argfun(@(x) cellfun(@(y) compute_stats(y, x, 'IgnoreNan', true), ...
                        allValues), ...
            'mean', 'std', 'stderr', 'err', 'upper95', 'lower95');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
