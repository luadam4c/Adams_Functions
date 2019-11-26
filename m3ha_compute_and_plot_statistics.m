function m3ha_compute_and_plot_statistics (fitMode, inFolder, outFolder)
%% Plot bar graphs for LTS and burst statistics
% Usage: m3ha_compute_and_plot_statistics (fitMode, inFolder, outFolder)
%
% Arguments: 
%       fitMode     
%       inFolder    - (opt) the directory that contains the matfile to read
%                   must be a directory
%                   default == //media/adamX/m3ha/data_dclamp/take4/
%       outFolder   - (opt) the directory to output bar graphs
%                   (different subdirectories will be created for each fitMode)
%                   must be a directory
%                   default == //media/adamX/m3ha/data_dclamp/take4/
%       varargin    - 'FitMode': fitting mode
%                   must be one of:
%                           - 0 - all data
%                           - 1 - all of g incr = 100%, 200%, 400%
%                           - 2 - all of g incr = 100%, 200%, 400% 
%                                   but exclude cell-pharm-g_incr sets 
%                                   containing problematic sweeps
%                   default == 0
%
% Requires:
%       "inFolder"/dclampdatalog_take4.mat
%       cd/check_dir.m
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
% 2016-09-13 - Added fitMode
% 2016-09-14 - Added maxnoise & peakclass
% 2016-10-14 - Added suffix; changed the directory name for fitMode == 0 to include suffix ‘_all’
% 2016-10-15 - Made inFolder and outFolder optional arguments
% 2016-10-15 - Fixed error when passing fitMode == 0 to m3ha_find_ind_to_fit.m
% 2016-10-31 - Placed suffix into specs_for_fitmode.m
% 2017-01-24 - Plotted data points overlaying boxplots and bargraphs
% 2017-01-24 - Corrected the error bars on the bar graphs (it was half the value previously)
% 2017-01-25 - Corrected error bars on the bar graphs to reflect t-confidence intervals (from the Gosset's t distribution)
% 2018-02-04 - Copy matfile to /take4/ if the suffix is '_all'

%% Flags
debugFlag = false; %true;
bigFontFlag = true; %false;

%% Parameters used in analysis
sigLevel = 0.05;             % p value threshold for determining significance

%% Specify which matfile to use; assumed to be in inFolder
dataFileName = 'dclampdatalog_take4.csv';

figTypes = {'png', 'epsc2'};

% Items to compute
statsTitle = {'LTS onset time (ms)', 'LTS time jitter (ms)', ...
                'LTS probability', 'Spikes per LTS', ...
                'Burst onset time (ms)', 'Burst time jitter (ms)', ...
                'Burst probability', 'Spikes per burst'};
statsFileName = {'lts_onset_time', 'lts_time_jitter', ...
                    'lts_probability', 'spikes_per_lts', ...
                    'burst_onset_time', 'burst_time_jitter', ...
                    'burst_probability', 'spikes_per_burst'};

%% Fixed parameters used in the experiments
vv = [-60; -65; -70];       % Possible Vhold values (mV, LJP-corrected)
vvLabel = {'-60 mV', '-65 mV', '-70 mV'};
gg = [25; 50; 100; 200; 400; 800];  
                            % All possible conductance amplitude scaling %
ggLabelOrig = {'25%', '50%', '100%', '200%', '400%', '800%'};
ggLabelToFit = {'100%', '200%', '400%'};
pp = [1; 2; 3; 4];          % Possible pharm conditions (1 - Control; 2 - GAT1 Block; 3 - GAT3 Block; 4 - Dual Block)
ppLabel = {'Control', 'GAT1 Block', 'GAT3 Block', 'Dual Block'};
ppLabel2 = {'Con', 'GAT1', 'GAT3', 'Dual'};
cc = 1:1:49;                % Possible cell ID #s
ss = 1:1:5;                 % Possible within condition sweep #s

%% Default values for optional arguments
param1Default = [];             % default TODO: Description of fitMode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'fitMode', param1Default, ...
    % TODO: validation function %);

% Read from the Input Parser
parse(iP, varargin{:});
fitMode = iP.Results.fitMode;

%% Check arguments
if nargin < 1
    error('A fitMode is required, type ''help m3ha_compute_and_plot_statistics'' for usage');
elseif isempty(fitMode) || ~isnumeric(fitMode) || ~(fitMode == 0 || fitMode == 1 || fitMode == 2)
    error('fitMode out of range!');
elseif nargin >= 2 && ~isdir(inFolder)
    error('inFolder must be a directory!');
elseif nargin >= 3 && ~isdir(outFolder)
    error('outFolder must be a directory!');
end

%% Locate home directory
homeDirectory = m3ha_locate_homedir;

%% Set defaults for optional arguments
if nargin < 2
    if debugFlag
        inFolder = fullfile(homeDirectory, '/data_dclamp/take4/');
    else
        inFolder = fullfile(homeDirectory, '/data_dclamp/take4/');
    end
end
if nargin < 3
    if debugFlag
        outFolder = fullfile(homeDirectory, '/data_dclamp/take4/debug/');
    else
        outFolder = fullfile(homeDirectory, '/data_dclamp/take4/');
    end
end

%% Set font size
if bigFontFlag
    axisFontSize = 20;
    textFontSize = 20;
    ylabelFontSize = 20;
    markerSize = 20;
else
    axisFontSize = 10;
    textFontSize = 10;
    ylabelFontSize = 11;
    markerSize = 6;
end

%% Find path of data file to use
dataPath = fullfile(inFolder, dataFileName);

%% Print user specifications
fprintf('Using fit mode == %d ... \n', fitMode);
fprintf('Using data file == %s ... \n', dataPath);
fprintf('Using significane level == %g ... \n', sigLevel);

%% Set suffix according to fitMode
suffix = m3ha_specs_for_fitmode(fitMode);

%% Set conductance amplitude scaling % labels for each fitMode
if fitMode == 0
    gglabel = ggLabelOrig;
elseif fitMode == 1 || fitMode == 2
    gglabel = ggLabelToFit;
end

%% Set folders for reading and saving files
outFolderBar = fullfile(outFolder, ['bargraphs', suffix]);
check_dir(outFolderBar);

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
all_stats_vgp = cell(1, length(statsTitle));
mean_stats_vgp = cell(1, length(statsTitle));
std_stats_vgp = cell(1, length(statsTitle));
ct_stats_vgp = cell(1, length(statsTitle));
err_stats_vgp = cell(1, length(statsTitle));
highbar_stats_vgp = cell(1, length(statsTitle));
lowbar_stats_vgp = cell(1, length(statsTitle));
for bi = 1:length(statsTitle)
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
if fitMode ~= 0
    indtofit = m3ha_find_ind_to_fit(fnrow, cellID, pharm, gincr, fitMode, inFolder);
end
for vi = 1:length(vv)
    v_ind = find(vhold == vv(vi));
    for gi = 1:length(gg)
        if (fitMode == 1 || fitMode == 2) && ~(gi == 3 || gi == 4 || gi == 5)
            continue;
        end
        g_ind = find(gincr == gg(gi));
        for hi = 1:length(pp)
            p_ind = find(pharm == pp(hi));
            if fitMode == 0
                vgp_ind = intersect(intersect(v_ind, g_ind), p_ind);
            else
                vgp_ind = intersect(intersect(intersect(v_ind, g_ind), p_ind), indtofit);
            end
            % Find/calculate LTS/burst statistics for each group
            [all_stats, mean_stats, std_stats, ct_stats, err_stats, highbar_stats, lowbar_stats] = ...
                m3ha_compute_ltsburst_statistics(vgp_ind, cellID, ltspeaktime, ...
                            spikesperpeak, bursttime, spikesperburst);
            for bi = 1:length(statsTitle)
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
                    if fitMode == 0
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
all_stats_gp = cell(1, length(statsTitle));
mean_stats_gp = cell(1, length(statsTitle));
std_stats_gp = cell(1, length(statsTitle));
ct_stats_gp = cell(1, length(statsTitle));
err_stats_gp = cell(1, length(statsTitle));
highbar_stats_gp = cell(1, length(statsTitle));
lowbar_stats_gp = cell(1, length(statsTitle));
for bi = 1:length(statsTitle)
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
    if (fitMode == 1 || fitMode == 2) && ~(gi == 3 || gi == 4 || gi == 5)
        continue;
    end
    g_ind = find(gincr == gg(gi));
    for hi = 1:length(pp)
        p_ind = find(pharm == pp(hi));
        if fitMode == 0
            gp_ind = intersect(g_ind, p_ind);
        else
            gp_ind = intersect(intersect(g_ind, p_ind), indtofit);
        end
        % Find/calculate LTS/burst statistics for each group
        [all_stats, mean_stats, std_stats, ct_stats, err_stats, highbar_stats, lowbar_stats] = ...
            m3ha_compute_ltsburst_statistics(gp_ind, cellID, ltspeaktime, spikesperpeak, bursttime, spikesperburst);
        for bi = 1:length(statsTitle)
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
pvalue_vg = zeros(length(statsTitle), length(gg), length(vv));
                        % p value across pharm conditions for each vg-condition
pvalue_g = zeros(length(statsTitle), length(gg));
                        % p value across pharm conditions for each g incr
table_vg = cell(length(statsTitle), length(gg), length(vv));
table_g = cell(length(statsTitle), length(gg));
%{ 
%For Multiple comparison test (multcompare(stats))
stats_vg = cell(length(statsTitle), length(gg), length(vv));
stats_g = cell(length(statsTitle), length(gg));
%}
ct_vg = zeros(length(pp), 1);
ct_g = zeros(length(pp), 1);
for bi = 1:length(statsTitle)
    for gi = 1:length(gg)
        % If under fitMode 1 & 2, only do this for G incr = 100%, 200%, 400%
        if (fitMode == 1 || fitMode == 2) && ~(gi == 3 || gi == 4 || gi == 5)
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
        if bigFontFlag
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
                        'k.', 'MarkerSize', markerSize);
                end
            end
        end
        set(ax, 'XTickLabel', ppLabel2);    % label x axis by pharm conditions
        set(ax, 'FontSize', axisFontSize);
        if pvalue_g(bi, gi) < sigLevel       % label the pvalue, make red if significant
            textcolor = 'r';
        else
            textcolor = 'k';
        end
        text(0.751, ax.YLim(2)*0.95, ['p value = ', num2str(pvalue_g(bi, gi))], ...
            'Color', textcolor, 'FontSize', textFontSize);
        ylabel(statsTitle{bi}, 'FontSize', ylabelFontSize);
        if ~bigFontFlag
            xlabel('Pharm Condition');
            title(['All traces with G scaled at ', num2str(gg(gi)), '%']);
        end
        figname = fullfile(outFolderBar, [statsFileName{bi}, '_', num2str(gg(gi)), 'g_boxplot', suffix, '.png']);
        save_all_figtypes(h, figname, figTypes);
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
                            'k.', 'MarkerSize', markerSize);
                    end
                end
            end
            set(ax, 'XTickLabel', ppLabel2);
            set(ax, 'FontSize', axisFontSize);
            if pvalue_g(bi, gi) < sigLevel            % label the pvalue, make red if significant
                textcolor = 'r';
            else
                textcolor = 'k';
            end
            text(0.751, ax.YLim(2)*0.95, ['p value = ', num2str(pvalue_vg(bi, gi, vi))], ...
                'Color', textcolor, 'FontSize', textFontSize);
            ylabel(statsTitle{bi}, 'FontSize', ylabelFontSize);
            if ~bigFontFlag
                xlabel('Pharm Condition');
                title(['Vhold = ', num2str(vv(vi)), ' mV with G scaled at ', num2str(gg(gi)), '%']);
            end
            figname = fullfile(outFolderBar, [statsFileName{bi}, ...
                '_', num2str(gg(gi)), 'g_v', num2str(vv(vi)), '_boxplot', suffix,'.png']);
            save_all_figtypes(h, figname, figTypes);
        end
    end
end

% Save variables as .mat file
fprintf('Saving variables ... \n');
matFile = fullfile(outFolderBar, ['ltsburst_statistics', suffix, '.mat']);
save(matFile, 'statsTitle', 'statsFileName', 'scpgv_ind', ...
    'all_stats_vgp', 'mean_stats_vgp', 'std_stats_vgp', ...
    'ct_stats_vgp', 'err_stats_vgp', 'highbar_stats_vgp', 'lowbar_stats_vgp', ...
    'all_stats_gp', 'mean_stats_gp', 'std_stats_gp', ...
    'ct_stats_gp', 'err_stats_gp', 'highbar_stats_gp', 'lowbar_stats_gp', ... 
    'pvalue_vg', 'pvalue_g', 'table_vg', 'table_g', ...
    '-v7.3');

% Copy .mat file if suffix is '_all'
if strcmp(suffix, '_all')
    copyfile(matFile, outFolder);
end

%% Create 3D bar graph for burst statistics for each Vhold value (Figure 3.4 in Christine's thesis)
fprintf('Plotting 3D bar graphs for burst statistics for each Vhold value ... \n');
parfor bi = 1:length(statsTitle)
    for vi = 1:length(vv)
        h = figure('Visible', 'off');
        set(h, 'Name', [statsFileName{bi}, '_', num2str(vv(vi))]);
        clf(h);
        % Plot means
        if fitMode == 0
            bar3(1:length(pp), mean_stats_vgp{bi}(:, :, vi), 0.12, 'detached'); hold on;
        elseif fitMode == 1 || fitMode == 2
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
        set(gca, 'YTickLabel', ppLabel);
         xlabel('IPSC conductance amplitude scaling');
%        ylabel('Pharm Condition');
        zlabel(statsTitle{bi});
        title(['Vm = ', vvLabel{vi}]);
        figname = fullfile(outFolderBar, [statsFileName{bi}, '_', num2str(vv(vi)), suffix, '.png']);
        save_all_figtypes(h, figname, figTypes);
        close(h);
    end
end

%% Create 2D bar graph for burst statistics for each Vhold value with IPSC conductance amplitude scaling fixed (Figure 3.5 in Christine's thesis)
fprintf('Plotting 2D bar graphs for burst statistics for each Vhold value with IPSC conductance amplitude scaling fixed ... \n');
for gi = 1:length(gg)
    % If under fitMode 1 & 2, only do this for G incr = 100%, 200%, 400%
    if (fitMode == 1 || fitMode == 2) && ~(gi == 3 || gi == 4 || gi == 5)
        continue;
    end        

    parfor bi = 1:length(statsTitle)
        h = figure('Visible', 'off');
        set(h, 'Name', [statsFileName{bi}, '_vsep_', num2str(gg(gi)), 'g']);
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
            set(ax, 'XTickLabel', ppLabel2);
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
            if pvalue_g(bi, gi) < sigLevel
                text(1, ax.YLim(2)*0.9, ['p value = ', num2str(pvalue_vg(bi, gi, vi))], 'Color', 'r')
            else
                text(1, ax.YLim(2)*0.9, ['p value = ', num2str(pvalue_vg(bi, gi, vi))], 'Color', 'k');
            end            
            xlabel('Pharm Condition');
            if vi == 1
                ylabel(statsTitle{bi});
            end
            title(['Vm = ', vvLabel{vi}]);
        end
        figname = fullfile(outFolderBar, [statsFileName{bi}, '_vsep_', num2str(gg(gi)), 'g', suffix, '.png']);
        save_all_figtypes(h, figname, figTypes);
        close(h);
    end
end

%% Create 3D bar graph for burst statistics
fprintf('Plotting 3D bar graphs for all burst statistics ... \n');
parfor bi = 1:length(statsTitle)
    h = figure('Visible', 'off');
    set(h, 'Name', statsFileName{bi});
    clf(h);
    % Plot means
    if fitMode == 0
        bar3(1:length(pp), mean_stats_gp{bi}(:, :), 0.12, 'detached'); hold on;
    elseif fitMode == 1 || fitMode == 2
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
    set(gca, 'YTickLabel', ppLabel);
     xlabel('IPSC conductance amplitude scaling');
%        ylabel('Pharm Condition');
    zlabel(statsTitle{bi});
    title('All traces');
    figname = fullfile(outFolderBar, [statsFileName{bi}, suffix, '.png']);
    save_all_figtypes(h, figname, figTypes);
    close(h);
end

%% Create 2D bar graph for burst statistics with IPSC conductance amplitude scaling fixed
fprintf('Plotting 2D bar graphs for all burst statistics with IPSC conductance amplitude scaling fixed ... \n');
for gi = 1:length(gg)
    % If under fitMode 1 & 2, only do this for G incr = 100%, 200%, 400%
    if (fitMode == 1 || fitMode == 2) && ~(gi == 3 || gi == 4 || gi == 5)
        continue;
    end        

    parfor bi = 1:length(statsTitle)
        % Count the number of cells in each pharm group
        ct_g = cellfun(@length, all_stats_gp{bi}(:, gi));

        % Create a 2D bar graph
        h = figure('Visible', 'off');
        set(h, 'Name', [statsFileName{bi}, '_', num2str(gg(gi)), 'g']);
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
        set(ax, 'XTickLabel', ppLabel2);
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
        if pvalue_g(bi, gi) < sigLevel
            text(1, ax.YLim(2)*0.9, ['p value = ', num2str(pvalue_g(bi, gi))], 'Color', 'r');
        else
            text(1, ax.YLim(2)*0.9, ['p value = ', num2str(pvalue_g(bi, gi))], 'Color', 'k');
        end
        xlabel('Pharm Condition');
        ylabel(statsTitle{bi});
        title(['All traces with G scaled at ', num2str(gg(gi)), '%']);
        figname = fullfile(outFolderBar, [statsFileName{bi}, '_', num2str(gg(gi)), 'g', suffix, '.png']);
        save_all_figtypes(h, figname, figTypes);
        close(h);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

m = matfile(dataPath);
dataFileName = 'dclampdatalog_take4.mat';

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%