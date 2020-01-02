function statsTable = m3ha_compute_statistics (varargin)
%% Computes LTS and burst statistics for some indices (indOfInterest) in ltsOnsetTimeEachSwp, spikesPerLtsEachSwp, burstOnsetTimeEachSwp & spikesPerBurst
% Usage: statsTable = m3ha_compute_statistics (varargin)
%
% Explanation:
%       TODO
%
% Example(s):
%       statsTable = m3ha_compute_statistics
%       statsTable = m3ha_compute_statistics('DataMode', 2)
%       statsTable = m3ha_compute_statistics('PharmConditions', 1)
%       statsTable = m3ha_compute_statistics('PharmConditions', 1:4)
%       statsTable = m3ha_compute_statistics('PharmConditions', num2cell(1:4))
%       statsTable = m3ha_compute_statistics('PharmConditions', num2cell(1:4), 'GIncrCondition', num2cell([100; 200; 400]))
%
% Requires:
%       cd/argfun.m
%       cd/array_fun.m
%       cd/compute_stats.m
%       cd/first_matching_field.m
%       cd/force_column_cell.m
%       cd/force_column_vector.m
%       cd/match_row_count.m
%       cd/m3ha_load_sweep_info.m
%       cd/m3ha_select_sweeps.m
%
% Used by:
%       cd/m3ha_compute_and_plot_statistics.m
%       cd/m3ha_plot_figure02.m
%       cd/m3ha_plot_figure04.m
%       cd/m3ha_simulate_population.m
%
% Outputs:
%       statsTable  - a table containing measures as row names
%                       and the following variables:
%                           measureTitle
%                           measureStr
%                           dataMode
%                           pharmCondition
%                           gIncrCondition
%                           vHoldCondition
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
% 2019-12-04 Added other LTS features that might be of interest:
% 2019-12-04 Added other burst features that might be of interest
% TODO: Add 'MeasuresToCompute' as an optional argument

%% Hard-coded parameters
% Items to compute
measureTitle = {'LTS amplitude (mV)'; 'LTS maximum slope (V/s)'; ...
                'LTS concavity (V^2/s^2)'; 'LTS prominence (mv)'; ...
                'LTS width (ms)'; 'LTS onset time (ms)'; ...
                'LTS time jitter (ms)'; ...
                'LTS probability'; 'Spikes per LTS'; ...
                'Maximum spike amplitude (mV)'; ...
                'Minimum spike amplitude (mV)'; ...
                'Spike frequency (Hz)'; 'Spike adaptation'; ...
                'Burst onset time (ms)'; 'Burst time jitter (ms)'; ...
                'Burst probability'; 'Spikes per burst'};
%   Note: measureStr must match the variables defined in this file!
measureStr = {'ltsAmplitude'; 'ltsMaxSlope'; ...
                'ltsConcavity'; 'ltsProminence'; ...
                'ltsWidth'; 'ltsOnsetTime'; 'ltsTimeJitter'; ...
                'ltsProbability'; 'spikesPerLts'; ...
                'spikeMaxAmp'; 'spikeMinAmp'; ...
                'spikeFrequency'; 'spikeAdaptation'
                'burstOnsetTime'; 'burstTimeJitter'; ...
                'burstProbability'; 'spikesPerBurst'};

% Columns defined in swpInfo
cellIdStr = 'cellidrow';
ltsAmplitudeStr = {'ltspeakval', 'ltsPeakValue'};
ltsMaxSlopeStr = {'maxslopeval', 'maxSlopeValue'};
ltsConcavityStr = {'peak2ndder', 'peak2ndDer'};
ltsProminenceStr = {'peakprom', 'peakProm'};
ltsWidthStr = {'peakwidth', 'peakWidth'};
ltsOnsetTimeStr = {'ltspeaktime', 'ltsPeakTime'};
spikesPerLtsStr = {'spikesperpeak', 'spikesPerPeak'};
spikeMaxAmpStr = {'maxspikeamp', 'maxSpikeAmp'};
spikeMinAmpStr = {'minspikeamp', 'minSpikeAmp'};
spikeFrequencyStr = {'spikefrequency', 'spikeFrequency'};
spikeAdaptationStr = {'spikeadaptation', 'spikeAdaptation'};
burstOnsetTimeStr = {'bursttime', 'burstTime'};
spikesPerBurstStr = {'spikesperburst', 'spikesPerBurst'};

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
nMeasures = numel(measureStr);

%% Preparation
% Read the sweep info data
if isempty(swpInfo)
    swpInfo = m3ha_load_sweep_info;
end

%% Create all conditions
if iscell(pharmConditionsUser) || iscell(gIncrConditionsUser) || ...
        iscell(vHoldConditionsUser)
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

    % Mesh indices for all set of conditions
    %   Note: meshgrid returns a length(y) by length(x) by length(z) matrix
    %           so one needs to permute the first two dimensions
    %           to give a length(x) by length(y) by length(z) matrix
    % TODO: meshgrid_custom.m
    [indPharmCondition, indGIncrCondition, indVHoldCondition] = ...
        meshgrid(1:nPharm, 1:nGIncr, 1:nVHold);
    [indPharmCondition, indGIncrCondition, indVHoldCondition] = ...
        argfun(@(x) permute(x, [2, 1]), ...
                indPharmCondition, indGIncrCondition, indVHoldCondition);

    % Mesh all set of conditions
    pharmCondition = arrayfun(@(x) pharmConditionsUser{x}, ...
                                indPharmCondition, 'UniformOutput', false);
    gIncrCondition = arrayfun(@(x) gIncrConditionsUser{x}, ...
                                indGIncrCondition, 'UniformOutput', false);
    vHoldCondition = arrayfun(@(x) vHoldConditionsUser{x}, ...
                                indVHoldCondition, 'UniformOutput', false);

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
        array_fun(@(x, y, z) ...
                    m3ha_compute_statistics_helper(swpInfo, dataMode, ...
                        x, y, z, cellIdStr, measureStr, ...
                        ltsAmplitudeStr, ltsMaxSlopeStr, ltsConcavityStr, ...
                        ltsProminenceStr, ltsWidthStr, ltsOnsetTimeStr, ...
                        spikesPerLtsStr, spikeMaxAmpStr, spikeMinAmpStr, ...
                        spikeFrequencyStr, spikeAdaptationStr, ...
                        burstOnsetTimeStr, spikesPerBurstStr), ...
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
                        transpose(1:nMeasures), 'UniformOutput', false), ...
            nValues, meanValue, stdValue, ...
            stderrValue, errValue, upper95Value, lower95Value);
else
    % Compute statistics for this set of conditions
    %   Note: this will return either a nMeasures x 1 column cell vector
    %            or a nMeasures x 1 column numeric vector
    [allValues, nValues, meanValue, stdValue, ...
            stderrValue, errValue, upper95Value, lower95Value] = ...
        m3ha_compute_statistics_helper(swpInfo, dataMode, ...
                pharmCondition, gIncrCondition, vHoldCondition, ...
                cellIdStr, measureStr, ...
                ltsAmplitudeStr, ltsMaxSlopeStr, ltsConcavityStr, ...
                ltsProminenceStr, ltsWidthStr, ltsOnsetTimeStr, ...
                spikesPerLtsStr, spikeMaxAmpStr, spikeMinAmpStr, ...
                spikeFrequencyStr, spikeAdaptationStr, ...
                burstOnsetTimeStr, spikesPerBurstStr);
end

%% Output results
% Put in a cell array
[pharmCondition, gIncrCondition, vHoldCondition] = ...
    argfun(@(x) {x}, pharmCondition, gIncrCondition, vHoldCondition);

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
                m3ha_compute_statistics_helper(swpInfo, dataMode, ...
                        pharmConditions, gIncrConditions, vHoldConditions, ...
                        cellIdStr, measureStr, ...
                        ltsAmplitudeStr, ltsMaxSlopeStr, ltsConcavityStr, ...
                        ltsProminenceStr, ltsWidthStr, ltsOnsetTimeStr, ...
                        spikesPerLtsStr, spikeMaxAmpStr, spikeMinAmpStr, ...
                        spikeFrequencyStr, spikeAdaptationStr, ...
                        burstOnsetTimeStr, spikesPerBurstStr)
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
[ltsAmplitudeEachSwp, ltsMaxSlopeEachSwp, ltsConcavityEachSwp, ...
        ltsProminenceEachSwp, ltsWidthEachSwp, ltsOnsetTimeEachSwp, ...
        spikesPerLtsEachSwp, spikeMaxAmpEachSwp, spikeMinAmpEachSwp, ...
        spikeFrequencyEachSwp, spikeAdaptationEachSwp, ...
        burstOnsetTimeEachSwp, spikesPerBurstEachSwp] = ...
    argfun(@(x) first_matching_field(swpInfoToUse, x), ...
            ltsAmplitudeStr, ltsMaxSlopeStr, ltsConcavityStr, ...
            ltsProminenceStr, ltsWidthStr, ltsOnsetTimeStr, ...
            spikesPerLtsStr, spikeMaxAmpStr, spikeMinAmpStr, ...
            spikeFrequencyStr, spikeAdaptationStr, ...
            burstOnsetTimeStr, spikesPerBurstStr);

% Determine whether each sweep has an LTS or has a burst
hasLts = ltsOnsetTimeEachSwp > 0;
hasBurst = burstOnsetTimeEachSwp > 0;

% Set features to be NaN if there is no LTS
%   Note: This may not be necessary for most features but just in case
%           This is necessary for concavity, prominence and width as of 20191204
%           TODO: Fix and run m3ha_append_lts_properties.m
ltsAmplitudeEachSwp(~hasLts) = NaN;
ltsMaxSlopeEachSwp(~hasLts) = NaN;
ltsConcavityEachSwp(~hasLts) = NaN;
ltsProminenceEachSwp(~hasLts) = NaN;
ltsWidthEachSwp(~hasLts) = NaN;
ltsOnsetTimeEachSwp(~hasLts) = NaN;
spikesPerLtsEachSwp(~hasLts) = NaN;
spikeMaxAmpEachSwp(~hasLts) = NaN;
spikeMinAmpEachSwp(~hasLts) = NaN;
spikeFrequencyEachSwp(~hasLts) = NaN;
spikeAdaptationEachSwp(~hasLts) = NaN;
burstOnsetTimeEachSwp(~hasLts) = NaN;
spikesPerBurstEachSwp(~hasLts) = NaN;

%% Compute LTS & burst measures for each cell
% Compute the LTS probability for each cell
ltsProbability = arrayfun(@(x) sum(cellIdRow == x & hasLts) / ...
                                sum(cellIdRow == x), uniqueCellIds);

% Compute the burst probability for each cell
burstProbability = arrayfun(@(x) sum(cellIdRow == x & hasBurst) / ...
                                sum(cellIdRow == x), uniqueCellIds);

% Compute means of LTS and burst properties for each cell
[ltsAmplitude, ltsMaxSlope, ltsConcavity, ...
        ltsProminence, ltsWidth, ltsOnsetTime, ...
        spikesPerLts, spikeMaxAmp, spikeMinAmp, ...
        spikeFrequency, spikeAdaptation, burstOnsetTime, spikesPerBurst] = ...
    argfun(@(x) arrayfun(@(y) nanmean(x(cellIdRow == y)), uniqueCellIds), ...
            ltsAmplitudeEachSwp, ltsMaxSlopeEachSwp, ltsConcavityEachSwp, ...
            ltsProminenceEachSwp, ltsWidthEachSwp, ltsOnsetTimeEachSwp, ...
            spikesPerLtsEachSwp, spikeMaxAmpEachSwp, spikeMinAmpEachSwp, ...
            spikeFrequencyEachSwp, spikeAdaptationEachSwp, ...
            burstOnsetTimeEachSwp, spikesPerBurstEachSwp);

% Compute standard deviations of LTS and burst properties for each cell
[ltsTimeJitter, burstTimeJitter] = ...
    argfun(@(x) arrayfun(@(y) nanstd(x(cellIdRow == y)), uniqueCellIds), ...
            ltsOnsetTimeEachSwp, burstOnsetTimeEachSwp);

%% Calculate overall LTS & burst statistics
% Collect all values for each measure
%   Note: Each element of the cell array contains a nCells by 1 numeric vector
%           that may contain NaNs
%   Note: cellfun won't work here because of namespace differences
nMeasures = numel(measureStr);
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

swpInfoStrs = {ltsAmplitudeStr; ltsMaxSlopeStr; ltsConcavityStr; ...
                ltsProminenceStr; ltsWidthStr; ltsOnsetTimeStr; ...
                spikesPerLtsStr; spikeMaxAmpStr; spikeMinAmpStr; ...
                spikeFrequencyStr; spikeAdaptationStr; ...
                burstOnsetTimeStr; spikesPerBurstStr};
swpInfoToUse.(x)

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
