function statsTable = m3ha_compute_ltsburst_statistics_new (varargin)
%% Computes LTS and burst statistics for some indices (indOfInterest) in ltsPeakTimes, spikesPerPeak, burstTimes & spikesPerBurst
% Usage: statsTable = m3ha_compute_ltsburst_statistics_new (varargin)
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
%                           allValue
%                           meanValue
%                           stdValue
%                           errValue
%                           upper95Values
%                           lower95Values
%                           nValues
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
%                       cellidrow   - cell ID
%                       prow        - pharmacological condition
%                       grow        - conductance amplitude scaling
%                   default == m3ha_load_sweep_info

% File History:
% 2016-08-19 Created
% 2016-08-29 Last Modified
% 2017-01-25 - Corrected errors to reflect t-confidence intervals 
%               (from the Gosset's t distribution)
% 2019-11-26 Improved code structure

%% Hard-coded parameters
% Items to compute
measureTitle = {'LTS onset time (ms)'; 'LTS time jitter (ms)'; ...
                'LTS probability', 'Spikes per LTS'; ...
                'Burst onset time (ms)'; 'Burst time jitter (ms)'; ...
                'Burst probability'; 'Spikes per burst'};
measureStr = {'ltsOnsetTime', 'ltsTimeJitter', ...
                    'ltsProbability', 'spikesPerLts', ...
                    'burstOnsetTime', 'burstTimeJitter', ...
                    'burstProbability', 'spikesPerBurst'};

%% Default values for optional arguments
dataModeDefault = 0;
swpInfoDefault = table.empty;

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

% Read from the Input Parser
parse(iP, varargin{:});
dataMode = iP.Results.DataMode;
swpInfo = iP.Results.SwpInfo;

% Count the items to compute
nMeasures = numel(measureTitle);

%% Preparation
% Read the sweep info data
if isempty(swpInfo)
    swpInfo = m3ha_load_sweep_info('FileName', dataPath);
end

%% Select sweeps
% Select the sweeps based on data mode
swpInfo = m3ha_select_sweeps_to_fit('DataMode', dataMode);

% Extract whether to use the sweep
toUse = swpInfo.toFit;

% Restrict to those sweeps
swpInfoToUse = swpInfo(toUse, :);

%% Determine all possible cells
% Extract the cell IDs
cellId = swpInfoToUse.cellidrow;

% Find unique cell IDs
uniqueCellIds = unique(cellId);

% Count the number of cells
nCells = numel(uniqueCellIds);

%% Extract measures
% Extract the measures
ltsPeakTimes = swpInfoToUse.ltspeaktimes;
spikesPerPeak = swpInfoToUse.spikesperpeak;
burstTimes = swpInfoToUse.bursttimes;
spikesPerBurst = swpInfoToUse.spikesperburst;

% Determine whether each sweep has an LTS or has a burst
hasLts = ltsPeakTimes > 0;
hasBurst = burstTimes > 0;

%% Compute LTS & burst measures for each cell
% Compute the LTS probability
ltsProbability = arrayfun(@(x) sum(cellId == x & hasLts) / ...
                                sum(cellId == x), uniqueCellIds);

% Compute means of LTS and burst properties
[ltsOnsetTime, spikesPerLts, burstOnsetTime, spikesPerBurst] = ...
    argfun(@(x) arrayfun(@(y) nanmean(x(cellId == y)), ...
                        uniqueCellIds), ...
            ltsPeakTimes, spikesPerPeak, burstTimes, spikesPerBurst);

% 

measureStr = {'ltsOnsetTime', 'ltsTimeJitter', ...
                    'ltsProbability', 'spikesPerLts', ...
                    'burstOnsetTime', 'burstTimeJitter', ...
                    'burstProbability', 'spikesPerBurst'};

%% Initialize stats vectors
allValues = cell(nMeasures, 1);
meanValues = zeros(nMeasures, 1);
stdValues = zeros(nMeasures, 1);
errValues = zeros(nMeasures, 1);
upper95Values = zeros(nMeasures, 1);
lower95Values = zeros(nMeasures, 1);
nValues = zeros(nMeasures, 1);




ct1 = 0;    % counts cells in this condition
ct2 = 0;    % counts cells in this condition with LTSs
ct3 = 0;    % counts cells in this condition with bursts
for ci = 1:length(cc)
    indThisCell = find(uniqueCellIds == cc(ci));
    indThisCellOfInterest = intersect(indOfInterest, indThisCell);
    indThisCellOfInterestHasLts = intersect(indThisCellOfInterest, indHasLts);
    indThisCellOfInterestHasBurst = intersect(indThisCellOfInterest, indHasBurst);
    if ~isempty(indThisCellOfInterest)
        ct1 = ct1 + 1;
        allValues{3}(ct1, 1) = length(indThisCellOfInterestHasLts)/length(indThisCellOfInterest);
        allValues{7}(ct1, 1) = length(indThisCellOfInterestHasBurst)/length(indThisCellOfInterest);
    end
    if ~isempty(indThisCellOfInterestHasLts)
        ct2 = ct2 + 1;
        allValues{1}(ct2, 1) = mean(ltsPeakTimes(indThisCellOfInterestHasLts));
        allValues{2}(ct2, 1) = std(ltsPeakTimes(indThisCellOfInterestHasLts));
        allValues{4}(ct2, 1) = mean(spikesPerPeak(indThisCellOfInterestHasLts));
    end
    if ~isempty(indThisCellOfInterestHasBurst)
        ct3 = ct3 + 1;
        allValues{5}(ct3, 1) = mean(burstTimes(indThisCellOfInterestHasBurst));
        allValues{6}(ct3, 1) = std(burstTimes(indThisCellOfInterestHasBurst));
        allValues{8}(ct3, 1) = mean(spikesPerBurst(indThisCellOfInterestHasBurst));
    end
end

%% Calculate overall LTS & burst statistics
if ~isempty(indOfInterestHasLts)
    meanValues(1) = mean(ltsPeakTimes(indOfInterestHasLts));    % mean of LTS onset time (ms)
    stdValues(1) = std(ltsPeakTimes(indOfInterestHasLts));
                % standard deviation of LTS onset time (ms) over all trials 
%    stdValues(1) = std(allValues{1});
                % standard deviation of LTS onset time (ms) over all cells
    meanValues(2) = mean(allValues{2});    % mean of LTS time jitter (ms)
    stdValues(2) = std(allValues{2});    % standard deviation of LTS time jitter (ms)
    meanValues(3) = mean(allValues{3});    % mean of LTS probability
    stdValues(3) = std(allValues{3});    % standard deviation of LTS probability
    meanValues(4) = mean(spikesPerPeak(indOfInterestHasLts));    % mean of spikes per LTS
    stdValues(4) = std(spikesPerPeak(indOfInterestHasLts));
                % standard deviation of spikes per LTS over all trials 
%    stdValues(4) = std(allValues{4});
                % standard deviation of spikes per LTS over all cells
end
if ~isempty(indOfInterestHasBurst)
    meanValues(5) = mean(burstTimes(indOfInterestHasBurst));    % mean of burst onset time (ms)
    stdValues(5) = std(burstTimes(indOfInterestHasBurst));
                % standard deviation of burst onset time (ms) over all trials 
%    stdValues(5) = std(allValues{1});
                % standard deviation of burst onset time (ms) over all cells
    meanValues(6) = mean(allValues{6});    % mean of burst time jitter (ms)
    stdValues(6) = std(allValues{6});    % standard deviation of burst time jitter (ms)
    meanValues(7) = mean(allValues{7});    % mean of burst probability
    stdValues(7) = std(allValues{7});    % standard deviation of burst probability
    meanValues(8) = mean(spikesPerBurst(indOfInterestHasBurst));    % mean of spikes per burst
    stdValues(8) = std(spikesPerBurst(indOfInterestHasBurst));
                % standard deviation of spikes per burst over all trials 
%    stdValues(8) = std(allValues{4});
                % standard deviation of spikes per burst over all cells
end

%% Calculate 95% confidence intervals of the mean
for bi = 1:nMeasures
    nValues(bi) = length(allValues{bi});
    if nValues(bi) > 0
        % Use 2-sided 95% t-confidence intervals 
        %   (based on Gosset's t distribution with n-1 degrees of freedom)
        errValues(bi) = tinv(0.975, nValues(bi)-1) * stdValues(bi) / sqrt(nValues(bi));
    end
    upper95Values(bi) = meanValues(bi) + errValues(bi);
    lower95Values(bi) = meanValues(bi) - errValues(bi);
end

%% Output results
statsTable = table(measureTitle, measureStr, allValues, nValues, ...
                    meanValues, stdValues, errValues, upper95Values, ...
                    lower95Values, 'RowNames', measureStr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

errValues(bi) = 1.96 * stdValues(bi) / sqrt(nValues(bi));    % 95% confidence intervals

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
