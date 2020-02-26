function [aDrugBurstProbNullValue, aDrugBurstProbAltValue, burstProbStdev] = ...
                metabolismR01_extrapolate_burst_probability
%% Extrapolates burst probability values for A-Drug from GAT1 data
%
% Requires:
%       cd/argfun.m

% File History:
% 2020-02-26 Created by Adam Lu

% Paths
grantDir = '/media/shareX/2020marchR01';
powerDir = fullfile(grantDir, 'Power_Analysis');

% Files
gatStatsFile = 'gat_blockade_pharm_1-4_gincr_200_stats.mat';
aDrugOscFile = 'adrug_oscDurationSec_chevron.csv';
gat1OscFile = 'gat1_oscDurationSec_chevron.csv';

% Strings
controlStr = 'Con';
gat1Str = 'GAT1';
burstProbStr = 'burstProbability';
allValuesStr = 'allValues';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load GAT Blocker stats
m = matfile(fullfile(powerDir, gatStatsFile));
statsTable = m.statsTable;
pharmLabels = m.pharmLabels;

% Determine Control vs. GAT1 index
[idxCon, idxGat1] = ...
    argfun(@(str) find(contains(pharmLabels, str)), controlStr, gat1Str);

% Extract Control vs. GAT1 burst probabilities
burstProbAll = statsTable{burstProbStr, allValuesStr};
burstProbAll = burstProbAll{1};
[conBurstProbValues, gat1BurstProbValues] = ...
    argfun(@(idx) burstProbAll{idx}, idxCon, idxGat1);

% Compute Control vs. GAT1 mean burst probabilities
[conBurstMeanProb, gat1BurstMeanProb] = ...
    argfun(@(allValues) nanmean(allValues), ...
            conBurstProbValues, gat1BurstProbValues);

% Read oscillation table
[aDrugOscTable, gat1OscTable] = ...
    argfun(@(file) readtable(fullfile(powerDir, file)), ...
            aDrugOscFile, gat1OscFile);

% Calculate oscillation duration effect sizes
[aDrugOscEffect, gat1OscEffect] = ...
    argfun(@(table) nanmean(table{:, 2} - table{:, 1}), ...
            aDrugOscTable, gat1OscTable);

% Estimate A-Drug burst probability null value
aDrugBurstProbNullValue = conBurstMeanProb;

% Estimate A-Drug burst probability alt value
aDrugBurstProbAltValue = conBurstMeanProb + ...
                            (gat1BurstMeanProb - conBurstMeanProb) * ...
                            (aDrugOscEffect / gat1OscEffect);

% Estimate burst probability difference standard deviation
burstProbStdev = nanstd(gat1BurstProbValues - conBurstProbValues);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Compute Control vs. GAT1 burst probability standard deviations
[conBurstStdProb, gat1BurstStdProb] = ...
    argfun(@(allValues) nanstd(allValues), ...
            conBurstProbValues, gat1BurstProbValues);
% Estimate burst probability standard deviation
burstProbStdev = sqrt(conBurstStdProb^2 + gat1BurstStdProb^2);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%