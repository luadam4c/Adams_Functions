function [output1] = m3ha_compute_and_compare_lts_statistics (ltsTablesSim, ltsTablesReal, varargin)
%% Computes statistics for LTS features and compares across (NOT FINISHED)
% Usage: [output1] = m3ha_compute_and_compare_lts_statistics (ltsTablesSim, ltsTablesReal, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
%
% Arguments:
%       reqarg1     - TODO: Description of reqarg1
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       TODO
%       cd/plot_bar.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-11-17 Moved code from m3ha_neuron_run_and_analyze.m but not organized
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'reqarg1');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, reqarg1, varargin{:});
% param1 = iP.Results.param1;

%% Preparation
% TODO

%% Do the job
% TODO

%% Output results
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ltst_all, ltst_mean_real, ltst_mean_sim, ltst_mean_all, ...
    ltst_std_real, ltst_std_sim, ltst_cv_real, ltst_cv_sim, ltst_cv_all, ltst_max_dev_all, ...
    ltsdvdtv_all, ltsdvdtv_mean_real, ltsdvdtv_mean_sim, ltsdvdtv_mean_all, ...
    ltsdvdtv_std_real, ltsdvdtv_std_sim, ltsdvdtv_cv_real, ltsdvdtv_cv_sim, ...
    ltsdvdtv_cv_all, ltsdvdtv_max_dev_all] ...
        = compute_and_compare_statistics(nSweeps, colorMap, ncg, npercg, ...
            cellID, outparams, plotStatisticsFlag, ...
            real_ipeakt, real_ltsv, real_ltst, real_ltsdvdtv, real_ltsdvdtt, ...
            real_pkprom, real_pkwidth, real_pkclass, real_np2der, real_spp, real_btime, real_spb, ...
            sim_ipeakt, sim_ltsv, sim_ltst, sim_ltsdvdtv, sim_ltsdvdtt, ...
            sim_pkprom, sim_pkwidth, sim_pkclass, sim_np2der, sim_spp, sim_btime, sim_spb)
%% Calculate statistics and plot those of simulated data against real data
%% TODO: Make functions out of this

%% Extract from outparams
outFolder = outparams.outFolder;
prefix = outparams.prefix;

%% Compute LTS statistics
ltsp = zeros(ncg, 2);            % LTS probability
ltsv_real = cell(ncg, 1);
ltst_real = cell(ncg, 1);
ltsdvdtv_real = cell(ncg, 1);
ltsdvdtt_real = cell(ncg, 1);
ltsv_sim = cell(ncg, 1);
ltst_sim = cell(ncg, 1);
ltsdvdtv_sim = cell(ncg, 1);
ltsdvdtt_sim = cell(ncg, 1);
both_has_lts_ct = zeros(ncg, 1);
ltst_all = cell(ncg, 1);
ltst_mean_real = zeros(ncg, 1);
ltst_mean_sim = zeros(ncg, 1);
ltst_mean_all = zeros(ncg, 1);
ltst_std_real = zeros(ncg, 1);
ltst_std_sim = zeros(ncg, 1);
ltst_std_all = zeros(ncg, 1);
ltst_cv_real = zeros(ncg, 1);
ltst_cv_sim = zeros(ncg, 1);
ltst_cv_all = zeros(ncg, 1);
ltst_max_dev_all = zeros(ncg, 1);
ltsdvdtv_all = cell(ncg, 1);
ltsdvdtv_mean_real = zeros(ncg, 1);
ltsdvdtv_mean_sim = zeros(ncg, 1);
ltsdvdtv_mean_all = zeros(ncg, 1);
ltsdvdtv_std_real = zeros(ncg, 1);
ltsdvdtv_std_sim = zeros(ncg, 1);
ltsdvdtv_std_all = zeros(ncg, 1);
ltsdvdtv_cv_real = zeros(ncg, 1);
ltsdvdtv_cv_sim = zeros(ncg, 1);
ltsdvdtv_cv_all = zeros(ncg, 1);
ltsdvdtv_max_dev_all = zeros(ncg, 1);
for cgn = 1:ncg            % color group number
    % LTS probability
    ct_real = 0;    % counts sweeps with LTS
    ct_sim = 0;    % counts sweeps with LTS
    for iSwp = 1:nSweeps
        if ceil(iSwp/npercg) == cgn && ~isnan(real_ltst(iSwp))
            ct_real = ct_real + 1;
        end
        if ceil(iSwp/npercg) == cgn && ~isnan(sim_ltst(iSwp))
            ct_sim = ct_sim + 1;
        end
    end
    ltsp(cgn, 1) = ct_real/npercg;
    ltsp(cgn, 2) = ct_sim/npercg;

    % Data for those sweeps that have LTSs both in real and simulated conditions
    ct = 0;         % counts sweeps with LTSs
    for iSwp = 1:nSweeps
        if ceil(iSwp/npercg) == cgn ...
            && ~isnan(real_ltst(iSwp)) && ~isnan(sim_ltst(iSwp))  
            ct = ct + 1;
            ltsv_real{cgn}(1, ct) = real_ltsv(iSwp);
            ltst_real{cgn}(1, ct) = real_ltst(iSwp);
            ltsdvdtv_real{cgn}(1, ct) = real_ltsdvdtv(iSwp);
            ltsdvdtt_real{cgn}(1, ct) = real_ltsdvdtt(iSwp);
            ltsv_sim{cgn}(1, ct) = sim_ltsv(iSwp);
            ltst_sim{cgn}(1, ct) = sim_ltst(iSwp);
            ltsdvdtv_sim{cgn}(1, ct) = sim_ltsdvdtv(iSwp);
            ltsdvdtt_sim{cgn}(1, ct) = sim_ltsdvdtt(iSwp);
        end
    end
    both_has_lts_ct(cgn) = ct;
    if both_has_lts_ct(cgn) > 0
        % LTS peak time data
        ltst_all{cgn} = [ltst_real{cgn} ltst_sim{cgn}];
        ltst_mean_real(cgn) = mean(ltst_real{cgn});
        ltst_mean_sim(cgn) = mean(ltst_sim{cgn});
        ltst_mean_all(cgn) = mean(ltst_all{cgn});
        ltst_std_real(cgn) = std(ltst_real{cgn});
        ltst_std_sim(cgn) = std(ltst_sim{cgn});
        ltst_std_all(cgn) = std([ltst_real{cgn} ltst_sim{cgn}]);
        ltst_cv_real(cgn) = ltst_std_real(cgn)/ltst_mean_real(cgn);
        ltst_cv_sim(cgn) = ltst_std_sim(cgn)/ltst_mean_sim(cgn);
        ltst_cv_all(cgn) = ltst_std_all(cgn)/ltst_mean_all(cgn);
        ltst_max_dev_all(cgn) = max(abs(ltst_all{cgn} - ltst_mean_all(cgn))/ltst_mean_all(cgn));

        % LTS max slope data
        ltsdvdtv_all{cgn} = [ltsdvdtv_real{cgn} ltsdvdtv_sim{cgn}];
        ltsdvdtv_mean_real(cgn) = mean(ltsdvdtv_real{cgn});
        ltsdvdtv_mean_sim(cgn) = mean(ltsdvdtv_sim{cgn});
        ltsdvdtv_mean_all(cgn) = mean(ltsdvdtv_all{cgn});
        ltsdvdtv_std_real(cgn) = std(ltsdvdtv_real{cgn});
        ltsdvdtv_std_sim(cgn) = std(ltsdvdtv_sim{cgn});
        ltsdvdtv_std_all(cgn) = std([ltsdvdtv_real{cgn} ltsdvdtv_sim{cgn}]);
        ltsdvdtv_cv_real(cgn) = ltsdvdtv_std_real(cgn)/ltsdvdtv_mean_real(cgn);
        ltsdvdtv_cv_sim(cgn) = ltsdvdtv_std_sim(cgn)/ltsdvdtv_mean_sim(cgn);
        ltsdvdtv_cv_all(cgn) = ltsdvdtv_std_all(cgn)/ltsdvdtv_mean_all(cgn);
        ltsdvdtv_max_dev_all(cgn) = max(abs(ltsdvdtv_all{cgn} - ltsdvdtv_mean_all(cgn))/ltsdvdtv_mean_all(cgn));
    end
end
if outparams.saveLtsStatsFlag
    save(fullfile(outFolder, [prefix, '_LTS_statistics.mat']), ...
        'ltst_all', 'ltst_mean_real', 'ltst_mean_sim', 'ltst_mean_all', ...
        'ltst_std_real', 'ltst_std_sim', 'ltst_cv_real', 'ltst_cv_sim', 'ltst_cv_all', 'ltst_max_dev_all', ...
        'ltsdvdtv_all', 'ltsdvdtv_mean_real', 'ltsdvdtv_mean_sim', 'ltsdvdtv_mean_all', ...
        'ltsdvdtv_std_real', 'ltsdvdtv_std_sim', 'ltsdvdtv_cv_real', 'ltsdvdtv_cv_sim', ...
        'ltsdvdtv_cv_all', 'ltsdvdtv_max_dev_all', '-v7.3');
end

%% Compute LTS/burst statistics

% Items to compute
statstitle = {'LTS onset time (ms)', 'LTS time jitter (ms)', 'LTS probability', 'Spikes per LTS', ...
        'Burst onset time (ms)', 'Burst time jitter (ms)', 'Burst probability', 'Spikes per burst'};
statsfilename = {'lts_onset_time', 'lts_time_jitter', 'lts_probability', 'spikes_per_lts', ...
                'burst_onset_time', 'burst_time_jitter', 'burst_probability', 'spikes_per_burst'};
pplabel2 = {'Con', 'GAT1', 'GAT3', 'Dual'};

% Initialize stats vectors
all_stats_real = cell(1, length(statstitle));
mean_stats_real = cell(1, length(statstitle));
std_stats_real = cell(1, length(statstitle));
ct_stats_real = cell(1, length(statstitle));
err_stats_real = cell(1, length(statstitle));
highbar_stats_real = cell(1, length(statstitle));
lowbar_stats_real = cell(1, length(statstitle));
all_stats_sim = cell(1, length(statstitle));
mean_stats_sim = cell(1, length(statstitle));
std_stats_sim = cell(1, length(statstitle));
ct_stats_sim = cell(1, length(statstitle));
err_stats_sim = cell(1, length(statstitle));
highbar_stats_sim = cell(1, length(statstitle));
lowbar_stats_sim = cell(1, length(statstitle));
for bi = 1:length(statstitle)
    all_stats_real{bi} = cell(ncg, 1);
    mean_stats_real{bi} = zeros(ncg, 1);
    std_stats_real{bi} = zeros(ncg, 1);
    ct_stats_real{bi} = zeros(ncg, 1);
    err_stats_real{bi} = zeros(ncg, 1);
    highbar_stats_real{bi} = zeros(ncg, 1);
    lowbar_stats_real{bi} = zeros(ncg, 1);

    all_stats_sim{bi} = cell(ncg, 1);
    mean_stats_sim{bi} = zeros(ncg, 1);
    std_stats_sim{bi} = zeros(ncg, 1);
    ct_stats_sim{bi} = zeros(ncg, 1);
    err_stats_sim{bi} = zeros(ncg, 1);
    highbar_stats_sim{bi} = zeros(ncg, 1);
    lowbar_stats_sim{bi} = zeros(ncg, 1);
end

for cgn = 1:ncg            % color group number
    thisp_ind = (cgn - 1) * npercg + (1:npercg);
    [all_stats, mean_stats, std_stats, ct_stats, err_stats, highbar_stats, lowbar_stats] = ...
        ltsburst_statistics(thisp_ind, cellID, real_ltst, real_spp, real_btime, real_spb);
%    [all_stats, mean_stats, std_stats, ct_stats, err_stats, highbar_stats, lowbar_stats] = ...
%        ltsburst_statistics(thisp_ind, cellID, outparams.ltspeaktime, outparams.spikesperpeak, outparams.bursttime, outparams.spikesperburst);
    for bi = 1:length(statstitle)
        all_stats_real{bi}{cgn} = all_stats{bi};
        mean_stats_real{bi}(cgn) = mean_stats(bi);
        std_stats_real{bi}(cgn) = std_stats(bi);
        ct_stats_real{bi}(cgn) = ct_stats(bi);
        err_stats_real{bi}(cgn) = err_stats(bi);
        highbar_stats_real{bi}(cgn) = highbar_stats(bi);
        lowbar_stats_real{bi}(cgn) = lowbar_stats(bi);
    end
    [all_stats, mean_stats, std_stats, ct_stats, err_stats, highbar_stats, lowbar_stats] = ...
        ltsburst_statistics(thisp_ind, cellID, sim_ltst, sim_spp, sim_btime, sim_spb);
    for bi = 1:length(statstitle)
        all_stats_sim{bi}{cgn} = all_stats{bi};
        mean_stats_sim{bi}(cgn) = mean_stats(bi);
        std_stats_sim{bi}(cgn) = std_stats(bi);
        ct_stats_sim{bi}(cgn) = ct_stats(bi);
        err_stats_sim{bi}(cgn) = err_stats(bi);
        highbar_stats_sim{bi}(cgn) = highbar_stats(bi);
        lowbar_stats_sim{bi}(cgn) = lowbar_stats(bi);
    end
end

% Plot statistics if plotStatisticsFlag == 1
if plotStatisticsFlag

    % Plot bar graph comparing LTS probabilities
    fprintf('Plotting bar graph comparing LTS probabilities ...\n');
    if outparams.showStatisticsFlag
        hFig.ltstp_bar = figure(201);
    else
        hFig.ltstp_bar = figure('Visible', 'off');
    end
    set(hFig.ltstp_bar, 'Name', 'Low threshold spike probability');
    clf(hFig.ltstp_bar);
    bar(1:size(colorMap, 1), ltsp);
    legend('Real data', 'Simulated data');
    xlabel('Pharm condition #');
    ylabel('LTS probability');
    title('Low threshold spike probability');
    figName = fullfile(outFolder, [prefix, '_ltstp_bar.png']);
    save_all_figtypes(hFig.ltstp_bar, figName);

    % Plot scatter plot of LTS onset times, don't save yet (axes not fixed)
    fprintf('Plotting scatter plot of LTS onset times ...\n');
    if outparams.showStatisticsFlag
        hFig.ltstcorr = figure(202);
    else
        hFig.ltstcorr = figure('Visible', 'off');
    end
    set(hFig.ltstcorr, 'Name', 'LTS onset times (ms)');
    clf(hFig.ltstcorr);
    for cgn = 1:size(colorMap, 1)            % color group number
        if both_has_lts_ct(cgn) > 0
            if ncg == 4
                subplot(2, 2, cgn); hold on;
                title(['Pharm condition ', num2str(cgn)])
            elseif ncg == 12
                subplot(4, 3, cgn); hold on;
                title(['Pharm condition ', num2str(floor(cgn/3) + 1)])
            end
            xlabel('Real data')
            ylabel('Simulated data')
            plot(ltst_real{cgn}, ltst_sim{cgn}, 'LineStyle', 'none', ...
                'Marker', 'o', 'MarkerEdgeColor', colorMap(cgn, :), 'MarkerFaceColor', colorMap(cgn, :));
        end
    end
    title('Correlation of LTS onset times (ms)')

    % Plot scatter plot of LTS max slopes, don't save yet (axes not fixed)
    fprintf('Plotting scatter plot of LTS max slopes ...\n');
    if outparams.showStatisticsFlag
        hFig.ltsdvdtvcorr = figure(203);
    else
        hFig.ltsdvdtvcorr = figure('Visible', 'off');
    end
    set(hFig.ltsdvdtvcorr, 'Name', 'LTS max slopes (mV/ms)');
    clf(hFig.ltsdvdtvcorr);
    for cgn = 1:size(colorMap, 1)            % color group number
        if both_has_lts_ct(cgn) > 0
            if ncg == 4
                subplot(2, 2, cgn); hold on;
                title(['Pharm condition ', num2str(cgn)])
            elseif ncg == 12
                subplot(4, 3, cgn); hold on;
                title(['Pharm condition ', num2str(floor(cgn/3) + 1)])
            end
            xlabel('Real data')
            ylabel('Simulated data')
            plot(ltsdvdtv_real{cgn}, ltsdvdtv_sim{cgn}, 'LineStyle', 'none', ...
                'Marker', 'o', 'MarkerEdgeColor', colorMap(cgn, :), 'MarkerFaceColor', colorMap(cgn, :));
        end
    end
    title('Correlation of LTS max slopes (mV/ms)')
    

    % Fix axes to make the scales consistent among subplots
    if sum(both_has_lts_ct) ~= 0
        ltst_max_dev_all_max = max(ltst_max_dev_all);
        ltst_std_all_max = max(ltst_std_all);
        if ltst_max_dev_all_max > 3 * ltst_std_all_max
            width = ltst_std_all_max;
        else
            width = ltst_max_dev_all_max;
        end

        % Don't let width be 0
        if width == 0
            width = 0.1;
        end

        set(0, 'CurrentFigure', hFig.ltstcorr);
        for cgn = 1:ncg            % color group number
            if both_has_lts_ct(cgn) > 0
                if ncg == 4
                    subplot(2, 2, cgn);
                elseif ncg == 12
                    subplot(4, 3, cgn);
                end                
                xmin = ltst_mean_all(cgn) * (1 - 1.1 * width);
                xmax = ltst_mean_all(cgn) * (1 + 1.1 * width);
                ymin = xmin;
                ymax = xmax;
                axis([xmin, xmax, ymin, ymax]);
            end
        end
        figName = fullfile(outFolder, [prefix, '_ltstcorr.png']);
        save_all_figtypes(hFig.ltstcorr, figName);

        ltsdvdtv_max_dev_all_max = max(ltsdvdtv_max_dev_all);
        ltsdvdtv_std_all_max = max(ltsdvdtv_std_all);
        if ltsdvdtv_max_dev_all_max > 3 * ltsdvdtv_std_all_max
            width = ltsdvdtv_std_all_max;
        else
            width = ltsdvdtv_max_dev_all_max;
        end

        % Don't let width be 0
        if width == 0
            width = 0.1;
        end

        set(0, 'CurrentFigure' , hFig.ltsdvdtvcorr);
        for cgn = 1:ncg            % color group number
            if both_has_lts_ct(cgn) > 0
                if ncg == 4
                    subplot(2, 2, cgn);
                elseif ncg == 12
                    subplot(4, 3, cgn);
                end                
                xmin = ltsdvdtv_mean_all(cgn) * (1 - 1.1 * width);
                xmax = ltsdvdtv_mean_all(cgn) * (1 + 1.1 * width);
                ymin = xmin;
                ymax = xmax;
                axis([xmin, xmax, ymin, ymax]);
            end
        end
        figName = fullfile(outFolder, [prefix, '_ltsdvdtvcorr.png']);
        save_all_figtypes(hFig.ltsdvdtvcorr, figName);
    end

    %% Create 2D bar graphs for LTS/burst statistics
    ltsburst_stats = cell(1, length(statstitle));    % for parfor
    parfor bi = 1:length(statstitle)
        if outparams.showStatisticsFlag
            ltsburst_stats{bi} = figure(210 + bi);
        else
            ltsburst_stats{bi} = figure('Visible', 'off');
        end
        fprintf('2D bar graph for %s ...\n', statsfilename{bi});
        set(ltsburst_stats{bi}, 'Name', [statsfilename{bi}]);
        clf(ltsburst_stats{bi});

        % Plot means with 95% confidence intervals
        plot_bar([mean_stats_real{bi}, mean_stats_sim{bi}], ...
                    [lowbar_stats_real{bi}, lowbar_stats_sim{bi}], ...
                    [highbar_stats_real{bi}, highbar_stats_sim{bi}], ...
                    'PValues', (1:ncg)', 'PTickLabels', pplabel2, ...
                    'FigHandle', ltsburst_stats{bi});

        if bi == 1
%            ylim([0 2500]);
        elseif bi == 2
%            ylim([0 800]);
        elseif bi == 3
%            ylim([0 1]);
        elseif bi == 4
%            ylim([0 6]);
        end
        legend('Real data', 'Simulated data');
        xlabel('Pharm Condition');
%        ylabel(statstitle{bi});
        title(statstitle{bi});
        figName = fullfile(outFolder, [prefix, '_', statsfilename{bi}]);
%        save_all_figtypes(ltsburst_stats{bi}, figName, {'png', 'fig'});
        save_all_figtypes(ltsburst_stats{bi}, figName, {'png'});
    end

    fprintf('\n');
    hFig.ltsburst_stats = ltsburst_stats;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
