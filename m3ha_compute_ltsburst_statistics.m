function [all_stats, mean_stats, std_stats, ct_stats, err_stats, ...
            highbar_stats, lowbar_stats] = ...
                ltsburst_statistics (this_ind, cellID, ltspeaktime, ...
                                    spikesperpeak, bursttime, spikesperburst)
%% Computes LTS and burst statistics for some indices (this_ind) in ltspeaktime, spikesperpeak, bursttime & spikesperburst
%
% Used by:    
%        cd/m3ha_dclampdatalog_analyze.m

% File History:
% 2016-08-19 Created
% 2016-08-29 Last Modified
% 2017-01-25 - Corrected errors to reflect t-confidence intervals (from the Gosset's t distribution)

%% Fixed parameters used in the experiments
cc = 1:1:49;                % Possible cell ID #s

%% Items to compute
statstitle = {'LTS onset time (ms)', 'LTS time jitter (ms)', ...
                'LTS probability', 'Spikes per LTS', ...
                'Burst onset time (ms)', 'Burst time jitter (ms)', ...
                'Burst probability', 'Spikes per burst'};

%% Initialize stats vectors
all_stats = cell(1, length(statstitle));
mean_stats = zeros(1, length(statstitle));
std_stats = zeros(1, length(statstitle));
ct_stats = zeros(1, length(statstitle));
err_stats = zeros(1, length(statstitle));
highbar_stats = zeros(1, length(statstitle));
lowbar_stats = zeros(1, length(statstitle));

%% Find the indices with LTS & the indices with bursts
goodlts_ind = find(ltspeaktime > 0);
goodb_ind = find(bursttime > 0);
this_goodlts_ind = intersect(this_ind, goodlts_ind);
this_goodb_ind = intersect(this_ind, goodb_ind);

%% Calculate LTS & burst statistics for each cell
ct1 = 0;    % counts cells in this condition
ct2 = 0;    % counts cells in this condition with LTSs
ct3 = 0;    % counts cells in this condition with bursts
for ci = 1:length(cc)
    c_ind = find(cellID == cc(ci));
    thisc_ind = intersect(this_ind, c_ind);
    thisc_goodlts_ind = intersect(thisc_ind, goodlts_ind);
    thisc_goodb_ind = intersect(thisc_ind, goodb_ind);
    if ~isempty(thisc_ind)
        ct1 = ct1 + 1;
        all_stats{3}(ct1, 1) = length(thisc_goodlts_ind)/length(thisc_ind);
        all_stats{7}(ct1, 1) = length(thisc_goodb_ind)/length(thisc_ind);
    end
    if ~isempty(thisc_goodlts_ind)
        ct2 = ct2 + 1;
        all_stats{1}(ct2, 1) = mean(ltspeaktime(thisc_goodlts_ind));
        all_stats{2}(ct2, 1) = std(ltspeaktime(thisc_goodlts_ind));
        all_stats{4}(ct2, 1) = mean(spikesperpeak(thisc_goodlts_ind));
    end
    if ~isempty(thisc_goodb_ind)
        ct3 = ct3 + 1;
        all_stats{5}(ct3, 1) = mean(bursttime(thisc_goodb_ind));
        all_stats{6}(ct3, 1) = std(bursttime(thisc_goodb_ind));
        all_stats{8}(ct3, 1) = mean(spikesperburst(thisc_goodb_ind));
    end
end

%% Calculate overall LTS & burst statistics
if ~isempty(this_goodlts_ind)
    mean_stats(1) = mean(ltspeaktime(this_goodlts_ind));    % mean of LTS onset time (ms)
    std_stats(1) = std(ltspeaktime(this_goodlts_ind));
                % standard deviation of LTS onset time (ms) over all trials 
%    std_stats(1) = std(all_stats{1});
                % standard deviation of LTS onset time (ms) over all cells
    mean_stats(2) = mean(all_stats{2});    % mean of LTS time jitter (ms)
    std_stats(2) = std(all_stats{2});    % standard deviation of LTS time jitter (ms)
    mean_stats(3) = mean(all_stats{3});    % mean of LTS probability
    std_stats(3) = std(all_stats{3});    % standard deviation of LTS probability
    mean_stats(4) = mean(spikesperpeak(this_goodlts_ind));    % mean of spikes per LTS
    std_stats(4) = std(spikesperpeak(this_goodlts_ind));
                % standard deviation of spikes per LTS over all trials 
%    std_stats(4) = std(all_stats{4});
                % standard deviation of spikes per LTS over all cells
end
if ~isempty(this_goodb_ind)
    mean_stats(5) = mean(bursttime(this_goodb_ind));    % mean of burst onset time (ms)
    std_stats(5) = std(bursttime(this_goodb_ind));
                % standard deviation of burst onset time (ms) over all trials 
%    std_stats(5) = std(all_stats{1});
                % standard deviation of burst onset time (ms) over all cells
    mean_stats(6) = mean(all_stats{6});    % mean of burst time jitter (ms)
    std_stats(6) = std(all_stats{6});    % standard deviation of burst time jitter (ms)
    mean_stats(7) = mean(all_stats{7});    % mean of burst probability
    std_stats(7) = std(all_stats{7});    % standard deviation of burst probability
    mean_stats(8) = mean(spikesperburst(this_goodb_ind));    % mean of spikes per burst
    std_stats(8) = std(spikesperburst(this_goodb_ind));
                % standard deviation of spikes per burst over all trials 
%    std_stats(8) = std(all_stats{4});
                % standard deviation of spikes per burst over all cells
end

%% Calculate 95% confidence intervals of the mean
for bi = 1:length(statstitle)
    ct_stats(bi) = length(all_stats{bi});
    if ct_stats(bi) > 0
        % Use 2-sided 95% t-confidence intervals 
        %   (based on Gosset's t distribution with n-1 degrees of freedom)
        err_stats(bi) = tinv(0.975, ct_stats(bi)-1) * std_stats(bi) / sqrt(ct_stats(bi));
    end
    highbar_stats(bi) = mean_stats(bi) + err_stats(bi);
    lowbar_stats(bi) = mean_stats(bi) - err_stats(bi);
end

%{
%% OLD CODE

err_stats(bi) = 1.96 * std_stats(bi) / sqrt(ct_stats(bi));    % 95% confidence intervals

%}
