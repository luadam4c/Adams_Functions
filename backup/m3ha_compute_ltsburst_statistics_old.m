function [all_stats, mean_stats, std_stats, ct_stats, err_stats, ...
            highbar_stats, lowbar_stats] = ...
                m3ha_compute_statistics (indOfInterest, cellID, ltspeaktime, ...
                                    spikesperpeak, bursttime, spikesperburst)
%% Computes LTS and burst statistics for some indices (indOfInterest) in ltspeaktime, spikesperpeak, bursttime & spikesperburst
%
% Used by:    
%        cd/m3ha_compute_and_plot_statistics.m
%
% Arguments:
%        TODO:
%        swpInfo    - a table of sweep info, with each row named by 
%                       the matfile base containing the raw data
%                   must a 2D table with row names being file bases
%                       and with the fields:
%                       cellidrow   - cell ID
%                       prow        - pharmacological condition
%                       grow        - conductance amplitude scaling
%                   default == m3ha_load_sweep_info

% File History:
% 2016-08-19 Created
% 2016-08-29 Last Modified
% 2017-01-25 - Corrected errors to reflect t-confidence intervals (from the Gosset's t distribution)

%% Fixed parameters used in the experiments
cc = 1:1:49;                % Possible cell ID #s

%% Items to compute
measureTitle = {'LTS onset time (ms)', 'LTS time jitter (ms)', ...
                'LTS probability', 'Spikes per LTS', ...
                'Burst onset time (ms)', 'Burst time jitter (ms)', ...
                'Burst probability', 'Spikes per burst'};

% Count the items to compute
nMeasures = numel(measureTitle);


%% Initialize stats vectors
all_stats = cell(1, nMeasures);
mean_stats = zeros(1, nMeasures);
std_stats = zeros(1, nMeasures);
ct_stats = zeros(1, nMeasures);
err_stats = zeros(1, nMeasures);
highbar_stats = zeros(1, nMeasures);
lowbar_stats = zeros(1, nMeasures);

%% Find the indices with LTS & the indices with bursts
indHasLts = find(ltspeaktime > 0);
indHasBurst = find(bursttime > 0);
indOfInterestHasLts = intersect(indOfInterest, indHasLts);
indOfInterestHasBurst = intersect(indOfInterest, indHasBurst);

%% Calculate LTS & burst statistics for each cell
ct1 = 0;    % counts cells in this condition
ct2 = 0;    % counts cells in this condition with LTSs
ct3 = 0;    % counts cells in this condition with bursts
for ci = 1:length(cc)
    indThisCell = find(cellID == cc(ci));
    indThisCellOfInterest = intersect(indOfInterest, indThisCell);
    indThisCellOfInterestHasLts = intersect(indThisCellOfInterest, indHasLts);
    indThisCellOfInterestHasBurst = intersect(indThisCellOfInterest, indHasBurst);
    if ~isempty(indThisCellOfInterest)
        ct1 = ct1 + 1;
        all_stats{3}(ct1, 1) = length(indThisCellOfInterestHasLts)/length(indThisCellOfInterest);
        all_stats{7}(ct1, 1) = length(indThisCellOfInterestHasBurst)/length(indThisCellOfInterest);
    end
    if ~isempty(indThisCellOfInterestHasLts)
        ct2 = ct2 + 1;
        all_stats{1}(ct2, 1) = mean(ltspeaktime(indThisCellOfInterestHasLts));
        all_stats{2}(ct2, 1) = std(ltspeaktime(indThisCellOfInterestHasLts));
        all_stats{4}(ct2, 1) = mean(spikesperpeak(indThisCellOfInterestHasLts));
    end
    if ~isempty(indThisCellOfInterestHasBurst)
        ct3 = ct3 + 1;
        all_stats{5}(ct3, 1) = mean(bursttime(indThisCellOfInterestHasBurst));
        all_stats{6}(ct3, 1) = std(bursttime(indThisCellOfInterestHasBurst));
        all_stats{8}(ct3, 1) = mean(spikesperburst(indThisCellOfInterestHasBurst));
    end
end

%% Calculate overall LTS & burst statistics
if ~isempty(indOfInterestHasLts)
    mean_stats(1) = mean(ltspeaktime(indOfInterestHasLts));    % mean of LTS onset time (ms)
    std_stats(1) = std(ltspeaktime(indOfInterestHasLts));
                % standard deviation of LTS onset time (ms) over all trials 
%    std_stats(1) = std(all_stats{1});
                % standard deviation of LTS onset time (ms) over all cells
    mean_stats(2) = mean(all_stats{2});    % mean of LTS time jitter (ms)
    std_stats(2) = std(all_stats{2});    % standard deviation of LTS time jitter (ms)
    mean_stats(3) = mean(all_stats{3});    % mean of LTS probability
    std_stats(3) = std(all_stats{3});    % standard deviation of LTS probability
    mean_stats(4) = mean(spikesperpeak(indOfInterestHasLts));    % mean of spikes per LTS
    std_stats(4) = std(spikesperpeak(indOfInterestHasLts));
                % standard deviation of spikes per LTS over all trials 
%    std_stats(4) = std(all_stats{4});
                % standard deviation of spikes per LTS over all cells
end
if ~isempty(indOfInterestHasBurst)
    mean_stats(5) = mean(bursttime(indOfInterestHasBurst));    % mean of burst onset time (ms)
    std_stats(5) = std(bursttime(indOfInterestHasBurst));
                % standard deviation of burst onset time (ms) over all trials 
%    std_stats(5) = std(all_stats{1});
                % standard deviation of burst onset time (ms) over all cells
    mean_stats(6) = mean(all_stats{6});    % mean of burst time jitter (ms)
    std_stats(6) = std(all_stats{6});    % standard deviation of burst time jitter (ms)
    mean_stats(7) = mean(all_stats{7});    % mean of burst probability
    std_stats(7) = std(all_stats{7});    % standard deviation of burst probability
    mean_stats(8) = mean(spikesperburst(indOfInterestHasBurst));    % mean of spikes per burst
    std_stats(8) = std(spikesperburst(indOfInterestHasBurst));
                % standard deviation of spikes per burst over all trials 
%    std_stats(8) = std(all_stats{4});
                % standard deviation of spikes per burst over all cells
end

%% Calculate 95% confidence intervals of the mean
for bi = 1:nMeasures
    ct_stats(bi) = length(all_stats{bi});
    if ct_stats(bi) > 0
        % Use 2-sided 95% t-confidence intervals 
        %   (based on Gosset's t distribution with n-1 degrees of freedom)
        err_stats(bi) = tinv(0.975, ct_stats(bi)-1) * std_stats(bi) / sqrt(ct_stats(bi));
    end
    highbar_stats(bi) = mean_stats(bi) + err_stats(bi);
    lowbar_stats(bi) = mean_stats(bi) - err_stats(bi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

err_stats(bi) = 1.96 * std_stats(bi) / sqrt(ct_stats(bi));    % 95% confidence intervals

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
