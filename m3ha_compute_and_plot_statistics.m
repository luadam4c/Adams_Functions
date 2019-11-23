function m3ha_compute_and_plot_statistics (fitmode, infolder, outfolder)
%% Plot bar graphs for LTS and burst statistics
% Usage: m3ha_compute_and_plot_statistics (fitmode, infolder, outfolder)
% Arguments: 
%        fitmode        - 0 - all data
%                - 1 - all of g incr = 100%, 200%, 400%
%                - 2 - all of g incr = 100%, 200%, 400% 
%                    but exclude cell-pharm-g_incr sets containing problematic sweeps
%        infolder    - (opt) the directory that contains the matfile to read
%                must be a directory
%                default == //media/adamX/m3ha/data_dclamp/take4/
%        outfolder    - (opt) the directory to output bar graphs
%                    (different subdirectories will be created for each fitmode)
%                must be a directory
%                default == //media/adamX/m3ha/data_dclamp/take4/
% 
%
% Requires:    
%       "infolder"/dclampdatalog_take4.mat
%       cd/m3ha_compute_ltsburst_statistics.m
%       cd/m3ha_find_ind_to_fit.m
%       cd/m3ha_specs_for_fitmode.m
%       cd/m3ha_locate_homedir.m
%
% Used by:    
%       cd/m3ha_parse_dclamp_data.m
%

% File History:
% 2016-08-10 - Created by AL
% 2016-09-05 - Changed method of data import
% 2016-09-06 - Added ANOVA
% 2016-09-08 - Added sum(ct_g) ~= 0 condition 
% 2016-09-13 - Added fitmode
% 2016-09-14 - Added maxnoise & peakclass
% 2016-10-14 - Added suffix; changed the directory name for fitmode == 0 to include suffix ‘_all’
% 2016-10-15 - Made infolder and outfolder optional arguments
% 2016-10-15 - Fixed error when passing fitmode == 0 to m3ha_find_ind_to_fit.m
% 2016-10-31 - Placed suffix into specs_for_fitmode.m
% 2017-01-24 - Plotted data points overlaying boxplots and bargraphs
% 2017-01-24 - Corrected the error bars on the bar graphs (it was half the value previously)
% 2017-01-25 - Corrected error bars on the bar graphs to reflect t-confidence intervals (from the Gosset's t distribution)
% 2018-02-04 - Copy matfile to /take4/ if the suffix is '_all'

%% Flags
debugflag = 0;
bigfontflag = 1; %0;

%% Parameters used in analysis
sig_thr = 0.05;             % p value threshold for determining significance

%% Specify which matfile to use; assumed to be in infolder
filetouse = 'dclampdatalog_take4.mat';

% Items to compute
statstitle = {'LTS onset time (ms)', 'LTS time jitter (ms)', 'LTS probability', 'Spikes per LTS', ...
        'Burst onset time (ms)', 'Burst time jitter (ms)', 'Burst probability', 'Spikes per burst'};
statsfilename = {'lts_onset_time', 'lts_time_jitter', 'lts_probability', 'spikes_per_lts', ...
                'burst_onset_time', 'burst_time_jitter', 'burst_probability', 'spikes_per_burst'};

%% Fixed parameters used in the experiments
vv = [-60; -65; -70];       % Possible Vhold values (mV, LJP-corrected)
vvlabel = {'-60 mV', '-65 mV', '-70 mV'};
gg = [25; 50; 100; 200; 400; 800];  
                            % All possible conductance amplitude scaling %
gglabel_orig = {'25%', '50%', '100%', '200%', '400%', '800%'};
gglabel_tofit = {'100%', '200%', '400%'};
pp = [1; 2; 3; 4];          % Possible pharm conditions (1 - Control; 2 - GAT1 Block; 3 - GAT3 Block; 4 - Dual Block)
pplabel = {'Control', 'GAT1 Block', 'GAT3 Block', 'Dual Block'};
pplabel2 = {'Con', 'GAT1', 'GAT3', 'Dual'};
cc = 1:1:49;                % Possible cell ID #s
ss = 1:1:5;                 % Possible within condition sweep #s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 1
    error('A fitmode is required, type ''help m3ha_compute_and_plot_statistics'' for usage');
elseif isempty(fitmode) || ~isnumeric(fitmode) || ~(fitmode == 0 || fitmode == 1 || fitmode == 2)
    error('fitmode out of range!');
elseif nargin >= 2 && ~isdir(infolder)
    error('infolder must be a directory!');
elseif nargin >= 3 && ~isdir(outfolder)
    error('outfolder must be a directory!');
end

%% Locate home directory
homeDirectory = m3ha_locate_homedir;

%% Set defaults for optional arguments
if nargin < 2
    if debugflag
        infolder = fullfile(homedirectory, '/data_dclamp/take4/');
    else
        infolder = fullfile(homedirectory, '/data_dclamp/take4/');
    end
end
if nargin < 3
    if debugflag
        outfolder = fullfile(homedirectory, '/data_dclamp/take4/debug/');
    else
        outfolder = fullfile(homedirectory, '/data_dclamp/take4/');
    end
end

%% Set font size
if bigfontflag
    axisfontsize = 20;
    textfontsize = 20;
    ylabelfontsize = 20;
    markersize = 20;
else
    axisfontsize = 10;
    textfontsize = 10;
    ylabelfontsize = 11;
    markersize = 6;
end

%% Find path of matfile to use
fullmatfilepath = fullfile(infolder, filetouse);
m = matfile(fullmatfilepath);

%% Print user specifications
fprintf('Using fit mode == %d ... \n', fitmode);
fprintf('Using matfile == %s ... \n', fullmatfilepath);
fprintf('Using sig_thr == %g ... \n', sig_thr);

%% Set suffix according to fitmode
suffix = m3ha_specs_for_fitmode(fitmode);

%% Set conductance amplitude scaling % labels for each fitmode
if fitmode == 0
    gglabel = gglabel_orig;
elseif fitmode == 1 || fitmode == 2
    gglabel = gglabel_tofit;
end

%% Set folders for reading and saving files
outfolder_bar = fullfile(outfolder, ['/bargraphs', suffix, '/']);
if exist(outfolder_bar, 'dir') ~= 7
    mkdir(outfolder_bar);
    fprintf('Made directory %s \n', outfolder_bar);
end

%% Import data
fprintf('Importing data ... \n');
% For m3ha_find_ind_to_fit
fnrow = m.fnrow;
ntraces = numel(fnrow);

% For grouping conditions and for m3ha_find_ind_to_fit
cellID = m.cellidrow;
pharm = m.prow;
gincr = m.grow;

% For grouping conditions
vhold = m.vrow;
swpn = m.swpnrow;

% Statistics to be plotted
ltspeaktime = m.ltspeaktime;
spikesperpeak = m.spikesperpeak;

bursttime = m.bursttime;
spikesperburst = m.spikesperburst;

% Other LTS features that might be of interest
%%% TODO
ltspeakval = m.ltspeakval;
maxslopeval = m.maxslopeval;
ltspeak2ndder = m.ltspeak2ndder;
ltspeakprom = m.ltspeakprom;
ltspeakwidth = m.ltspeakwidth;

% Other burst features that might be of interest
%%% TODO
maxspikeamp = m.maxspikeamp;
minspikeamp = m.minspikeamp;
spikefrequency = m.spikefrequency;
spikeadaptation = m.spikeadaptation;

%% Find/calculate LTS & burst statistics
% Initialize vgp-grouped stats vectors
fprintf('Analyzing data grouped by Vhold-g incr-pharm ... \n');
all_stats_vgp = cell(1, length(statstitle));
mean_stats_vgp = cell(1, length(statstitle));
std_stats_vgp = cell(1, length(statstitle));
ct_stats_vgp = cell(1, length(statstitle));
err_stats_vgp = cell(1, length(statstitle));
highbar_stats_vgp = cell(1, length(statstitle));
lowbar_stats_vgp = cell(1, length(statstitle));
for bi = 1:length(statstitle)
    all_stats_vgp{bi} = cell(length(pp), length(gg), length(vv));
    mean_stats_vgp{bi} = zeros(length(pp), length(gg), length(vv));
    std_stats_vgp{bi} = zeros(length(pp), length(gg), length(vv));
    ct_stats_vgp{bi} = zeros(length(pp), length(gg), length(vv));
    err_stats_vgp{bi} = zeros(length(pp), length(gg), length(vv));
    highbar_stats_vgp{bi} = zeros(length(pp), length(gg), length(vv));
    lowbar_stats_vgp{bi} = zeros(length(pp), length(gg), length(vv));
end

% Group by Vhold, then by g incr, then by pharm condition
scpgv_ind = zeros(length(ss), length(cc), length(pp), length(gg), length(vv));
if fitmode ~= 0
    indtofit = m3ha_find_ind_to_fit(fnrow, cellID, pharm, gincr, fitmode, infolder);
end
for vi = 1:length(vv)
    v_ind = find(vhold == vv(vi));
    for gi = 1:length(gg)
        if (fitmode == 1 || fitmode == 2) && ~(gi == 3 || gi == 4 || gi == 5)
            continue;
        end
        g_ind = find(gincr == gg(gi));
        for hi = 1:length(pp)
            p_ind = find(pharm == pp(hi));
            if fitmode == 0
                vgp_ind = intersect(intersect(v_ind, g_ind), p_ind);
            else
                vgp_ind = intersect(intersect(intersect(v_ind, g_ind), p_ind), indtofit);
            end
            % Find/calculate LTS/burst statistics for each group
            [all_stats, mean_stats, std_stats, ct_stats, err_stats, highbar_stats, lowbar_stats] = ...
                m3ha_compute_ltsburst_statistics(vgp_ind, cellID, ltspeaktime, ...
                            spikesperpeak, bursttime, spikesperburst);
            for bi = 1:length(statstitle)
                all_stats_vgp{bi}{hi, gi, vi} = all_stats{bi};
                mean_stats_vgp{bi}(hi, gi, vi) = mean_stats(bi);
                std_stats_vgp{bi}(hi, gi, vi) = std_stats(bi);
                ct_stats_vgp{bi}(hi, gi, vi) = ct_stats(bi);
                err_stats_vgp{bi}(hi, gi, vi) = err_stats(bi);
                highbar_stats_vgp{bi}(hi, gi, vi) = highbar_stats(bi);
                lowbar_stats_vgp{bi}(hi, gi, vi) = lowbar_stats(bi);
            end

            % Find unique indices for future use
            for ci = 1:length(cc)
                c_ind = find(cellID == cc(ci));
                vgpc_ind = intersect(vgp_ind, c_ind);
                for si = 1:length(ss)
                    s_ind = find(swpn == ss(si));
                    if fitmode == 0
                        unique_ind = intersect(vgpc_ind, s_ind);
                    else
                        unique_ind = intersect(intersect(vgpc_ind, s_ind), indtofit);
                    end
                    if ~isempty(unique_ind)
                        scpgv_ind(si, ci, hi, gi, vi) = unique_ind;
                    end
                end
            end
        end
    end
end

% Initialize gp-grouped stats vectors
fprintf('Analyzing data grouped by g incr-pharm ... \n');
all_stats_gp = cell(1, length(statstitle));
mean_stats_gp = cell(1, length(statstitle));
std_stats_gp = cell(1, length(statstitle));
ct_stats_gp = cell(1, length(statstitle));
err_stats_gp = cell(1, length(statstitle));
highbar_stats_gp = cell(1, length(statstitle));
lowbar_stats_gp = cell(1, length(statstitle));
for bi = 1:length(statstitle)
    all_stats_gp{bi} = cell(length(pp), length(gg));
    mean_stats_gp{bi} = zeros(length(pp), length(gg));
    std_stats_gp{bi} = zeros(length(pp), length(gg));
    ct_stats_gp{bi} = zeros(length(pp), length(gg));
    err_stats_gp{bi} = zeros(length(pp), length(gg));
    highbar_stats_gp{bi} = zeros(length(pp), length(gg));
    lowbar_stats_gp{bi} = zeros(length(pp), length(gg));
end

% Group by g incr, then by pharm condition
for gi = 1:length(gg)
    if (fitmode == 1 || fitmode == 2) && ~(gi == 3 || gi == 4 || gi == 5)
        continue;
    end
    g_ind = find(gincr == gg(gi));
    for hi = 1:length(pp)
        p_ind = find(pharm == pp(hi));
        if fitmode == 0
            gp_ind = intersect(g_ind, p_ind);
        else
            gp_ind = intersect(intersect(g_ind, p_ind), indtofit);
        end
        % Find/calculate LTS/burst statistics for each group
        [all_stats, mean_stats, std_stats, ct_stats, err_stats, highbar_stats, lowbar_stats] = ...
            m3ha_compute_ltsburst_statistics(gp_ind, cellID, ltspeaktime, spikesperpeak, bursttime, spikesperburst);
        for bi = 1:length(statstitle)
            all_stats_gp{bi}{hi, gi} = all_stats{bi};
            mean_stats_gp{bi}(hi, gi) = mean_stats(bi);
            std_stats_gp{bi}(hi, gi) = std_stats(bi);
            ct_stats_gp{bi}(hi, gi) = ct_stats(bi);
            err_stats_gp{bi}(hi, gi) = err_stats(bi);
            highbar_stats_gp{bi}(hi, gi) = highbar_stats(bi);
            lowbar_stats_gp{bi}(hi, gi) = lowbar_stats(bi);
        end
    end
end

% For each statistic, perform ANOVA across pharmacological conditions and compute p values
fprintf('Performing ANOVA ... \n');
pvalue_vg = zeros(length(statstitle), length(gg), length(vv));
                        % p value across pharm conditions for each vg-condition
pvalue_g = zeros(length(statstitle), length(gg));
                        % p value across pharm conditions for each g incr
table_vg = cell(length(statstitle), length(gg), length(vv));
table_g = cell(length(statstitle), length(gg));
%{ 
%For Multiple comparison test (multcompare(stats))
stats_vg = cell(length(statstitle), length(gg), length(vv));
stats_g = cell(length(statstitle), length(gg));
%}
ct_vg = zeros(length(pp), 1);
ct_g = zeros(length(pp), 1);
for bi = 1:length(statstitle)
    for gi = 1:length(gg)
        % If under fitmode 1 & 2, only do this for G incr = 100%, 200%, 400%
        if (fitmode == 1 || fitmode == 2) && ~(gi == 3 || gi == 4 || gi == 5)
            continue;
        end

        % Do some preliminary math
        ct_g = cellfun(@length, all_stats_gp{bi}(:, gi));               % the number of cells in each pharm group
        med_stats_curr = cellfun(@median, all_stats_gp{bi}(:, gi));     % compute median values
        all_stats_curr = cell2mat(all_stats_gp{bi}(:, gi));             % vector of stats

        % If there is no data at all, continue.
        if sum(ct_g) == 0
            continue;
        end

        % Perform One-way ANOVA across pharm conditions for each g incr
        group = [];                        % group #s assigned to each pharm group, a column vector
        for hi = 1:length(pp)
            group = [group; ones(ct_g(hi), 1) * hi];
        end
        [pvalue_g(bi, gi), table_g{bi, gi}] = ...
            anova1(all_stats_curr, group, 'off');               % compute p-value of ANOVA
        h = figure('Visible', 'off');    
        b = boxplot(all_stats_curr, group, 'Notch', 'on');      % create boxplot
        if bigfontflag
            set(b, 'LineWidth', 3);
        end
        hold on;
        ax = gca;
        xtick = get(ax, 'XTick');           % get the location of tick labels on the x axis
        hin0 = 0;                           % counts nonzero pharm conditions
        for hi = 1:length(pp)
            if ct_g(hi) > 0
                hin0 = hin0 + 1;
                line(xtick(hin0)+1/8*[-1 1], med_stats_curr(hi)*[1 1], 'Color', 'r', 'LineWidth', 1);
                for i = 1:ct_g(hi)
                    % Plot a black dot for each cell
                    plot(xtick(hin0)+(i-(ct_g(hi)+1)/2)/(2*ct_g(hi)), all_stats_gp{bi}{hi, gi}(i), ...
                        'k.', 'MarkerSize', markersize);
                end
            end
        end
        set(ax, 'XTickLabel', pplabel2);    % label x axis by pharm conditions
        set(ax, 'FontSize', axisfontsize);
        if pvalue_g(bi, gi) < sig_thr       % label the pvalue, make red if significant
            textcolor = 'r';
        else
            textcolor = 'k';
        end
        text(0.751, ax.YLim(2)*0.95, ['p value = ', num2str(pvalue_g(bi, gi))], ...
            'Color', textcolor, 'FontSize', textfontsize);
        ylabel(statstitle{bi}, 'FontSize', ylabelfontsize);
        if ~bigfontflag
            xlabel('Pharm Condition');
            title(['All traces with G scaled at ', num2str(gg(gi)), '%']);
        end
        figname = fullfile(outfolder_bar, [statsfilename{bi}, '_', num2str(gg(gi)), 'g_boxplot', suffix, '.png']);
        saveas(h, figname);
        close(h);

        % Perform One-way ANOVA across pharm conditions for each vg-condition
        for vi = 1:length(vv)
            % Do some preliminary math
            ct_vg = cellfun(@length, all_stats_vgp{bi}(:, gi, vi));    % the number of cells in each pharm group
            med_stats_curr = cellfun(@median, all_stats_vgp{bi}(:, gi, vi));    % compute median values
            all_stats_curr = cell2mat(all_stats_vgp{bi}(:, gi, vi));    % vector of stats

            % If there is no data at all, continue.
            if sum(ct_vg) == 0
                continue;
            end

            group = [];            % group #s assigned to each pharm group, a column vector
            for hi = 1:length(pp)
                group = [group; ones(ct_vg(hi), 1) * hi];
            end
            [pvalue_vg(bi, gi, vi), table_vg{bi, gi, vi}] = ...
                anova1(all_stats_curr, group, 'off');        
            h = figure('Visible', 'off');
            boxplot(all_stats_curr, group, 'Notch', 'on');            % create boxplot
            ax = gca;
            xtick = get(ax, 'XTick');    % get the location of tick labels on the x axis
            hin0 = 0;            % counts nonzero pharm conditions
            for hi = 1:length(pp)
                if ct_vg(hi) > 0
                    hin0 = hin0 + 1;
                    line(xtick(hin0)+1/8*[-1 1], med_stats_curr(hi)*[1 1], ...
                        'Color', 'r', 'LineWidth', 1);    
                    for i = 1:ct_vg(hi)    % NOTE: this has to be ct_vg(hi), NOT ct_vg(hin0)
                        % Plot a black dot for each cell
                        plot(xtick(hin0)+(i-(ct_g(hi)+1)/2)/(2*ct_vg(hi)), ...
                            all_stats_vgp{bi}{hi, gi, vi}(i), ...
                            'k.', 'MarkerSize', markersize);
                    end
                end
            end
            set(ax, 'XTickLabel', pplabel2);
            set(ax, 'FontSize', axisfontsize);
            if pvalue_g(bi, gi) < sig_thr            % label the pvalue, make red if significant
                textcolor = 'r';
            else
                textcolor = 'k';
            end
            text(0.751, ax.YLim(2)*0.95, ['p value = ', num2str(pvalue_vg(bi, gi, vi))], ...
                'Color', textcolor, 'FontSize', textfontsize);
            ylabel(statstitle{bi}, 'FontSize', ylabelfontsize);
            if ~bigfontflag
                xlabel('Pharm Condition');
                title(['Vhold = ', num2str(vv(vi)), ' mV with G scaled at ', num2str(gg(gi)), '%']);
            end
            figname = fullfile(outfolder_bar, [statsfilename{bi}, ...
                '_', num2str(gg(gi)), 'g_v', num2str(vv(vi)), '_boxplot', suffix,'.png']);
            saveas(h, figname);
        end
    end
end

% Save variables as .mat file
fprintf('Saving variables ... \n');
matFile = fullfile(outfolder_bar, ['m3ha_compute_ltsburst_statistics', suffix, '.mat']);
save(matFile, 'statstitle', 'statsfilename', 'scpgv_ind', ...
    'all_stats_vgp', 'mean_stats_vgp', 'std_stats_vgp', ...
    'ct_stats_vgp', 'err_stats_vgp', 'highbar_stats_vgp', 'lowbar_stats_vgp', ...
    'all_stats_gp', 'mean_stats_gp', 'std_stats_gp', ...
    'ct_stats_gp', 'err_stats_gp', 'highbar_stats_gp', 'lowbar_stats_gp', ... 
    'pvalue_vg', 'pvalue_g', 'table_vg', 'table_g', ...
    '-v7.3');

% Copy .mat file if suffix is '_all'
if strcmp(suffix, '_all')
    copyfile(matFile, outfolder);
end

%% Create 3D bar graph for burst statistics for each Vhold value (Figure 3.4 in Christine's thesis)
fprintf('Plotting 3D bar graphs for burst statistics for each Vhold value ... \n');
parfor bi = 1:length(statstitle)
    for vi = 1:length(vv)
        h = figure('Visible', 'off');
        set(h, 'Name', [statsfilename{bi}, '_', num2str(vv(vi))]);
        clf(h);
        % Plot means
        if fitmode == 0
            bar3(1:length(pp), mean_stats_vgp{bi}(:, :, vi), 0.12, 'detached'); hold on;
        elseif fitmode == 1 || fitmode == 2
            bar3(1:length(pp), mean_stats_vgp{bi}(:, 3:5, vi), 0.12, 'detached'); hold on;
        end
        % Plot 95% confidence intervals
        for gi = 1:length(gglabel)    % x axis
            for hi = 1:length(pp)    % y axis
                barspan = [hi - 0.2; hi + 0.2];
                barspan2 = [mean_stats_vgp{bi}(hi, gi, vi); highbar_stats_vgp{bi}(hi, gi, vi)];
                line([gi; gi], barspan, ...
                    [highbar_stats_vgp{bi}(hi, gi, vi); highbar_stats_vgp{bi}(hi, gi, vi)], ...
                    'Color', 'k');
                line([gi; gi], [hi; hi], barspan2, 'Color', 'k');
            end
        end
        set(gca, 'XTickLabel', gglabel);
        set(gca, 'YTickLabel', pplabel);
         xlabel('IPSC conductance amplitude scaling');
%        ylabel('Pharm Condition');
        zlabel(statstitle{bi});
        title(['Vm = ', vvlabel{vi}]);
        figname = fullfile(outfolder_bar, [statsfilename{bi}, '_', num2str(vv(vi)), suffix, '.png']);
        saveas(h, figname);
        close(h);
    end
end

%% Create 2D bar graph for burst statistics for each Vhold value with IPSC conductance amplitude scaling fixed (Figure 3.5 in Christine's thesis)
fprintf('Plotting 2D bar graphs for burst statistics for each Vhold value with IPSC conductance amplitude scaling fixed ... \n');
for gi = 1:length(gg)
    % If under fitmode 1 & 2, only do this for G incr = 100%, 200%, 400%
    if (fitmode == 1 || fitmode == 2) && ~(gi == 3 || gi == 4 || gi == 5)
        continue;
    end        

    parfor bi = 1:length(statstitle)
        h = figure('Visible', 'off');
        set(h, 'Name', [statsfilename{bi}, '_vsep_', num2str(gg(gi)), 'g']);
        clf(h);
        oldpos = get(h, 'Position');
        newpos = [oldpos(1)-oldpos(3) oldpos(2)-oldpos(4) 3*oldpos(3) 2*oldpos(4)];
        set(h, 'Position', newpos);
        for vi = 1:length(vv)
            % Count the number of cells in each pharm group
            ct_vg = cellfun(@length, all_stats_vgp{bi}(:, gi, vi));

            % Create a 2D bar graph
            subplot(2, length(vv), vi);
            % Plot means
            bar(1:length(pp), mean_stats_vgp{bi}(:, gi, vi), 'c'); hold on;
            % Plot 95% confidence intervals
            errorbar(mean_stats_vgp{bi}(:, gi, vi), err_stats_vgp{bi}(:, gi, vi), 'k.')
            ax = gca;
            set(ax, 'XTickLabel', pplabel2);
            xtick = get(ax, 'XTick');            % get the location of tick labels on the x axis
            for hi = 1:length(pp)
                for i = 1:ct_vg(hi)
                    % Plot a red dot for each cell
                    plot(xtick(hi)+(i-(ct_vg(hi)+1)/2)/(2*ct_vg(hi)), ...
                            all_stats_vgp{bi}{hi, gi, vi}(i), ...
                            'r.', 'MarkerSize', 6);
                end
            end

            if bi == 1
                ylim([0 2500]);
            elseif bi == 2
                ylim([0 800]);
            elseif bi == 3
                ylim([0 1]);
            elseif bi == 4
                ylim([0 6]);
            end
            % Print p value
            ax = gca;
            if pvalue_g(bi, gi) < sig_thr
                text(1, ax.YLim(2)*0.9, ['p value = ', num2str(pvalue_vg(bi, gi, vi))], 'Color', 'r')
            else
                text(1, ax.YLim(2)*0.9, ['p value = ', num2str(pvalue_vg(bi, gi, vi))], 'Color', 'k');
            end            
            xlabel('Pharm Condition');
            if vi == 1
                ylabel(statstitle{bi});
            end
            title(['Vm = ', vvlabel{vi}]);
        end
        figname = fullfile(outfolder_bar, [statsfilename{bi}, '_vsep_', num2str(gg(gi)), 'g', suffix, '.png']);
        saveas(h, figname);
        close(h);
    end
end

%% Create 3D bar graph for burst statistics
fprintf('Plotting 3D bar graphs for all burst statistics ... \n');
parfor bi = 1:length(statstitle)
    h = figure('Visible', 'off');
    set(h, 'Name', statsfilename{bi});
    clf(h);
    % Plot means
    if fitmode == 0
        bar3(1:length(pp), mean_stats_gp{bi}(:, :), 0.12, 'detached'); hold on;
    elseif fitmode == 1 || fitmode == 2
        bar3(1:length(pp), mean_stats_gp{bi}(:, 3:5), 0.12, 'detached'); hold on;
    end
    % Plot 95% confidence intervals
    for gi = 1:length(gglabel)        % x axis
        for hi = 1:length(pp)        % y axis
            barspan = [hi - 0.2; hi + 0.2];
            barspan2 = [mean_stats_gp{bi}(hi, gi); highbar_stats_gp{bi}(hi, gi)];
            line([gi; gi], barspan, [highbar_stats_gp{bi}(hi, gi); highbar_stats_gp{bi}(hi, gi)], 'Color', 'k');
            line([gi; gi], [hi; hi], barspan2, 'Color', 'k');
        end
    end
    set(gca, 'XTickLabel', gglabel);
    set(gca, 'YTickLabel', pplabel);
     xlabel('IPSC conductance amplitude scaling');
%        ylabel('Pharm Condition');
    zlabel(statstitle{bi});
    title('All traces');
    figname = fullfile(outfolder_bar, [statsfilename{bi}, suffix, '.png']);
    saveas(h, figname);
    close(h);
end

%% Create 2D bar graph for burst statistics with IPSC conductance amplitude scaling fixed
fprintf('Plotting 2D bar graphs for all burst statistics with IPSC conductance amplitude scaling fixed ... \n');
for gi = 1:length(gg)
    % If under fitmode 1 & 2, only do this for G incr = 100%, 200%, 400%
    if (fitmode == 1 || fitmode == 2) && ~(gi == 3 || gi == 4 || gi == 5)
        continue;
    end        

    parfor bi = 1:length(statstitle)
        % Count the number of cells in each pharm group
        ct_g = cellfun(@length, all_stats_gp{bi}(:, gi));

        % Create a 2D bar graph
        h = figure('Visible', 'off');
        set(h, 'Name', [statsfilename{bi}, '_', num2str(gg(gi)), 'g']);
        clf(h);

        % Plot means
        bar(1:length(pp), mean_stats_gp{bi}(:, gi), 'c'); hold on;

        % Plot 95% confidence intervals
        errorbar(mean_stats_gp{bi}(:, gi), err_stats_gp{bi}(:, gi), 'k.')
        ax = gca;
        xtick = get(ax, 'XTick');            % get the location of tick labels on the x axis
        for hi = 1:length(pp)
            for i = 1:ct_g(hi)
                % Plot a red dot for each cell
                plot(xtick(hi)+(i-(ct_g(hi)+1)/2)/(2*ct_g(hi)), all_stats_gp{bi}{hi, gi}(i), ...
                    'r.', 'MarkerSize', 6);
            end
        end
        set(ax, 'XTickLabel', pplabel2);
        if bi == 1
            ylim([0 2500]);
        elseif bi == 2
            ylim([0 800]);
        elseif bi == 3
            ylim([0 1]);
        elseif bi == 4
            ylim([0 6]);
        end
        % Print p value
        if pvalue_g(bi, gi) < sig_thr
            text(1, ax.YLim(2)*0.9, ['p value = ', num2str(pvalue_g(bi, gi))], 'Color', 'r');
        else
            text(1, ax.YLim(2)*0.9, ['p value = ', num2str(pvalue_g(bi, gi))], 'Color', 'k');
        end
        xlabel('Pharm Condition');
        ylabel(statstitle{bi});
        title(['All traces with G scaled at ', num2str(gg(gi)), '%']);
        figname = fullfile(outfolder_bar, [statsfilename{bi}, '_', num2str(gg(gi)), 'g', suffix, '.png']);
        saveas(h, figname);
        close(h);
    end
end

%{

        %%% THIS IS WRONG!
        % Calculate the error bar lengths
        err_stats_gp2 = err_stats_gp{bi}(:, gi) / 2;    % the function errorbar plots 2 * err_stats_gp2
        errorbar(mean_stats_gp{bi}(:, gi), err_stats_gp2, 'k.')

if fitmode == 0
    outfolder = '//media/adamX/m3ha/data_dclamp/take4/bargraphs/';
elseif fitmode == 1
    outfolder = '//media/adamX/m3ha/data_dclamp/take4/bargraphs_100-400all/';
elseif fitmode == 2
    outfolder = '//media/adamX/m3ha/data_dclamp/take4/bargraphs_tofit/';
end

        if fitmode == 0
            figname = fullfile(outfolder, [statsfilename{bi}, '_', num2str(vv(vi)), '.png']);
        elseif fitmode == 1 || fitmode == 2
            figname = fullfile(outfolder, [statsfilename{bi}, '_', num2str(vv(vi)), '_tofit.png']);
        end

%% Import data
d = importdata(fullfile(infolder, 'dclampdatalog_take4.csv'));
logheader = d.textdata(1,:);
fn = d.textdata(2:end,1);
dat = d.data;

gabab_amp = m.gabab_amp;
gabab_Trise = m.gabab_Trise;
gabab_TfallFast = m.gabab_TfallFast;
gabab_TfallSlow = m.gabab_TfallSlow;
gabab_w = m.gabab_w;
currpulse = m.currpulse;
Rin = m.Rin;
ioffset_old = m.ioffset_old;
imint = m.imint;
imin = m.imin;
actVhold = m.actVhold;
maxnoise = m.maxnoise;
narrowpeaktime = m.narrowpeaktime;
narrowpeak2ndder = m.narrowpeak2ndder;
peakprom = m.peakprom;
peakclass = m.peakclass;

%% Extract data
% logheader = {'Data filename', 'Cell ID #', ...
%         'Pharm condition', 'Vhold (mV)', 'GABAB IPSC G incr (%)', 'Within condition sweep #', ...
%         'GABAB IPSC G amp (nS)', 'GABAB IPSC G Trise (ms)', 'GABAB IPSC G TfallFast (ms)', ...
%         'GABAB IPSC G TfallSlow (ms)', 'GABAB IPSC G w', ...
%         'Current pulse amplitude (pA)', 'Rinput (MOhm)', 'IPSC offset (ms)', 'IPSC peak time (ms)', 'IPSC peak amplitude (pA)', ...
%         'Actual Vhold (mV)', 'Narrowest peak time (ms)', 'Narrowest peak 2nd derivative (V^2/s^2)', 'Spikes per peak', ...
%         'LTS onset time (ms)', 'LTS amplitude (mV)', 'Burst onset time (ms)', 'Spikes per burst'};
cellID = dat(:, 1);
pharm = dat(:, 2);
vhold = dat(:, 3);
gincr = dat(:, 4);
swpn = dat(:, 5);
gabab_amp = dat(:, 6);
gabab_Trise = dat(:, 7);
gabab_TfallFast = dat(:, 8);
gabab_TfallSlow = dat(:, 9);
gabab_w = dat(:, 10);
currpulse = dat(:, 11);
Rin = dat(:, 12);
ioffset = dat(:, 13);
imint = dat(:, 14);
imin = dat(:, 15);
actVhold = dat(:, 16);
narrowpeaktime = dat(:, 17);
narrowpeak2ndder = dat(:, 18);
spikesperpeak = dat(:, 19);
ltspeaktime = dat(:, 20);
ltspeakval = dat(:, 21);
bursttime = dat(:, 22);
spikesperburst = dat(:, 23);

%% Find burst variance across trials with varying holding potentials
% Group by 49 cells, 4 pharm conditions, 6 G incr -> 1176 conditions
cc = 1:1:49;
pp = [1 2 3 4];
gg = [25 50 100 200 400 800];

ct = 0;      % count of formed groups
goodlts_ind = find(ltst > -1);
for ci = 1:length(cc)
    c_ind = find(cellID == cc(ci));
    for pi = 1:length(pp)
        p_ind = find(pharm == pp(pi));
        for gi = 1:length(gg)
            g_ind = find(gincr == gg(gi));
            cpg_goodlts_ind = intersect(intersect(intersect(c_ind, p_ind), g_ind), goodlts_ind);
            if ~isempty(cpg_goodlts_ind)
                ct = ct + 1;
                samecond_unordered{ct,1} = length(cpg_goodlts_ind);     % number of sweeps in the group
                samecond_unordered{ct,2} = cpg_goodlts_ind;             % original indices of sweeps
                samecond_unordered{ct,3} = std(ltst(cpg_goodlts_ind));  % burst onset time variance
                samecond_unordered{ct,4} = fn(cpg_goodlts_ind);         % file names of sweeps
                samecond_unordered{ct,5} = cc(ci);                      % cell ID
                samecond_unordered{ct,6} = swpn(cpg_goodlts_ind);       % within condition sweep no.
                samecond_unordered{ct,7} = vhold(cpg_goodlts_ind);      % Vhold (mV)
                samecond_unordered{ct,8} = gg(gi);                      % G incr (%)
                samecond_unordered{ct,9} = pp(pi);                      % pharm condition
                samecond_unordered{ct,10} = rinput(cpg_goodlts_ind);    % Rinput (MOhm)
                samecond_unordered{ct,11} = ltst(cpg_goodlts_ind);      % LTS peak time (ms)
                samecond_unordered{ct,12} = ltsv(cpg_goodlts_ind);      % LTS peak val (mV)
                samecond_unordered{ct,13} = ipeak(cpg_goodlts_ind);     % IPSC negative peak (pA)
            end
        end
    end
end

% Sort groups according to burst onset time variance
[trash rearrangement] = sort([samecond_unordered{:,3}], 'descend');
samecond = samecond_unordered(rearrangement,:);

% The following is obsolete
[ltst_std_max ltst_std_max_j] = max([samecond{:,3}]);
most_var_lts_ind = samecond{ltst_std_max_j,2};
most_var_lts_fn = fn(most_var_lts_ind);
most_var_lts_cellID = cellID(most_var_lts_ind);

% Specific cells
samecond_cell25 = samecond(find([samecond{:,5}] == 25),:);
samecond_cell43 = samecond(find([samecond{:,5}] == 43),:);
samecond_cell37 = samecond(find([samecond{:,5}] == 37),:);
samecond_cell29 = samecond(find([samecond{:,5}] == 29),:);
samecond_cell6 = samecond(find([samecond{:,5}] == 6),:);

            if fitmode == 0
                figname = fullfile(outfolder_bar, [statsfilename{bi}, ...
                    '_', num2str(gg(gi)), 'g_v', num2str(vv(vi)), '_boxplot.png']);
            elseif fitmode == 1 || fitmode == 2
            end

if fitmode == 1 || fitmode == 2
    gg = [100; 200; 400];            % Possible conductance amplitude scaling %

    scpgv_indtofit = zeros(length(ss), length(cc), length(pp), length(gg), length(vv));
                vgp_indtofit = intersect(vgp_ind, indtofit);
            if fitmode == 1 || fitmode == 2
                [all_stats, mean_stats, std_stats, ct_stats, err_stats, highbar_stats, lowbar_stats] = ...
                    m3ha_compute_ltsburst_statistics(vgp_indtofit, cellID, ...
                        ltspeaktime, spikesperpeak, bursttime, spikesperburst);
            else
            end
                    unique_indtofit = intersect(unique_ind, indtofit);
                    if ~isempty(unique_indtofit)
                        scpgv_indtofit(si, ci, hi, gi, vi) = unique_indtofit;
                    end

if fitmode == 0
elseif fitmode == 1 || fitmode == 2
    save(fullfile(outfolder_bar, 'm3ha_compute_ltsburst_statistics_tofit.mat'), 'statstitle', 'statsfilename', 'scpgv_indtofit', ...
        'all_stats_vgp', 'mean_stats_vgp', 'std_stats_vgp', ...
        'ct_stats_vgp', 'err_stats_vgp', 'highbar_stats_vgp', 'lowbar_stats_vgp', ...
        'all_stats_gp', 'mean_stats_gp', 'std_stats_gp', ...
        'ct_stats_gp', 'err_stats_gp', 'highbar_stats_gp', 'lowbar_stats_gp', ... 
        'pvalue_vg', 'pvalue_g', 'table_vg', 'table_g', ...
        '-v7.3');
end

            gp_indtofit = intersect(gp_ind, indtofit);
        if fitmode == 0
        elseif fitmode == 1 || fitmode == 2
            [all_stats, mean_stats, std_stats, ct_stats, err_stats, highbar_stats, lowbar_stats] = ...
                m3ha_compute_ltsburst_statistics(gp_indtofit, cellID, ltspeaktime, spikesperpeak, bursttime, spikesperburst);
        end

%}
