function m3ha_compute_and_plot_statistics (varargin)
%% Plot bar graphs for LTS and burst statistics
% Usage: m3ha_compute_and_plot_statistics (varargin)
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
%                   - 'Directory': working directory
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'InFolder': directory to read files from
%                   must be a string scalar or a character vector
%                   default == same as directory
%                   - 'OutFolder': directory to place output files
%                   must be a string scalar or a character vector
%                   default == same as inFolder
%                   - 'SwpInfo': a table of sweep info, with each row named by 
%                               the matfile base containing the raw data
%                   must a 2D table with row names being file bases
%                       and with the fields:
%                       cellidrow   - cell ID
%                       prow        - pharmacological condition
%                       grow        - conductance amplitude scaling
%                   default == m3ha_load_sweep_info
%                   
% Requires:
%       "inFolder"/dclampdatalog_take4.mat
%       cd/check_dir.m
%       cd/m3ha_compute_ltsburst_statistics.m
%       cd/m3ha_find_ind_to_fit.m TODO: Replace this
%       cd/m3ha_load_sweep_info.m
%       cd/m3ha_locate_homedir.m
%       cd/m3ha_select_sweeps_to_fit.m
%       cd/m3ha_specs_for_fitmode.m
%
% Used by:
%       cd/m3ha_parse_dclamp_data.m
%

% File History:
% 2016-08-10 - Created by AL
% 2016-09-05 - Changed method of data import
% 2016-09-06 - Added ANOVA
% 2016-09-08 - Added sum(ct_g) ~= 0 condition 
% 2016-09-13 - Added dataMode
% 2016-09-14 - Added maxnoise & peakclass
% 2016-10-14 - Added suffix; changed the directory name for dataMode == 0 to include suffix ‘_all’
% 2016-10-15 - Made inFolder and outFolder optional arguments
% 2016-10-15 - Fixed error when passing dataMode == 0 to m3ha_find_ind_to_fit.m
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
measureTitle = {'LTS onset time (ms)', 'LTS time jitter (ms)', ...
                'LTS probability', 'Spikes per LTS', ...
                'Burst onset time (ms)', 'Burst time jitter (ms)', ...
                'Burst probability', 'Spikes per burst'};
measureStr = {'lts_onset_time', 'lts_time_jitter', ...
                    'lts_probability', 'spikes_per_lts', ...
                    'burst_onset_time', 'burst_time_jitter', ...
                    'burst_probability', 'spikes_per_burst'};

%% Fixed parameters used in the experiments
vv = [-60; -65; -70];       % Possible Vhold values (mV, LJP-corrected)
vvLabel = {'-60 mV', '-65 mV', '-70 mV'};
gg = [25; 50; 100; 200; 400; 800];  
                            % All possible conductance amplitude scaling %
ggLabelOrig = {'25%', '50%', '100%', '200%', '400%', '800%'};
ggLabelToUse = {'100%', '200%', '400%'};
pp = [1; 2; 3; 4];          % Possible pharm conditions (1 - Control; 2 - GAT1 Block; 3 - GAT3 Block; 4 - Dual Block)
ppLabel = {'Control', 'GAT1 Block', 'GAT3 Block', 'Dual Block'};
ppLabel2 = {'Con', 'GAT1', 'GAT3', 'Dual'};
cc = 1:1:49;                % Possible cell ID #s
ss = 1:1:5;                 % Possible within condition sweep #s

%% Default values for optional arguments
dataModeDefault = 0;
directoryDefault = '';                  % set later
inFolderDefault = '';                   % set later
outFolderDefault = '';                  % set later
swpInfoDefault = table.empty;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'DataMode', dataModeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'InFolder', inFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SwpInfo', swpInfoDefault, ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Read from the Input Parser
parse(iP, varargin{:});
dataMode = iP.Results.DataMode;
directory = iP.Results.Directory;
inFolder = iP.Results.InFolder;
outFolder = iP.Results.OutFolder;
swpInfo = iP.Results.SwpInfo;

%% Preparation
% Set home directory
if isempty(directory)
    directory = fullfile(m3ha_locate_homedir, 'data_dclamp', 'take4');
end

% Set input directory
if isempty(inFolder)
    inFolder = directory;
end

% Set output directory
if debugFlag
    outFolder = fullfile(directory, 'debug');
else
    outFolder = directory;
end

% Set font size
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

% Find path to data file to use
dataPath = fullfile(inFolder, dataFileName);

% Display user specifications
fprintf('Using fit mode == %d ... \n', dataMode);
fprintf('Using data file == %s ... \n', dataPath);
fprintf('Using significance level == %g ... \n', sigLevel);

% Decide on suffix according to dataMode
suffix = m3ha_specs_for_fitmode(dataMode);

% Decide on conductance amplitude scaling % labels
if dataMode == 0
    gglabel = ggLabelOrig;
elseif dataMode == 1 || dataMode == 2
    gglabel = ggLabelToUse;
end

% Create folder for output figures
figFolder = fullfile(outFolder, strcat('bargraphs', suffix));
check_dir(figFolder);

%% Import data
fprintf('Importing data ... \n');

% Read the sweep info data
if isempty(swpInfo)
    swpInfo = m3ha_load_sweep_info('FileName', dataPath);
end

% Select the sweeps to fit based on data mode
[swpInfo, fileBasesToUse] = m3ha_select_sweeps_to_fit('DataMode', dataMode);

%% TODO: Organize below:

% For m3ha_find_ind_to_fit
fnrow = swpInfo.fnrow;

% For grouping conditions and for m3ha_find_ind_to_fit
cellID = swpInfo.cellidrow;
pharm = swpInfo.prow;
gincr = swpInfo.grow;

% For grouping conditions
vhold = swpInfo.vrow;
swpn = swpInfo.swpnrow;

% Statistics to be plotted
ltspeaktime = swpInfo.ltspeaktime;
spikesperpeak = swpInfo.spikesperpeak;

bursttime = swpInfo.bursttime;
spikesperburst = swpInfo.spikesperburst;

% Other LTS features that might be of interest
%%% TODO
ltspeakval = swpInfo.ltspeakval;
maxslopeval = swpInfo.maxslopeval;
ltspeak2ndder = swpInfo.ltspeak2ndder;
ltspeakprom = swpInfo.ltspeakprom;
ltspeakwidth = swpInfo.ltspeakwidth;

% Other burst features that might be of interest
%%% TODO
maxspikeamp = swpInfo.maxspikeamp;
minspikeamp = swpInfo.minspikeamp;
spikefrequency = swpInfo.spikefrequency;
spikeadaptation = swpInfo.spikeadaptation;

%% Find/calculate LTS & burst statistics
% Initialize vgp-grouped stats vectors
fprintf('Analyzing data grouped by Vhold-g incr-pharm ... \n');
all_stats_vgp = cell(1, length(measureTitle));
mean_stats_vgp = cell(1, length(measureTitle));
std_stats_vgp = cell(1, length(measureTitle));
ct_stats_vgp = cell(1, length(measureTitle));
err_stats_vgp = cell(1, length(measureTitle));
highbar_stats_vgp = cell(1, length(measureTitle));
lowbar_stats_vgp = cell(1, length(measureTitle));
for bi = 1:length(measureTitle)
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
if dataMode ~= 0
    indtofit = m3ha_find_ind_to_fit(fnrow, cellID, pharm, gincr, dataMode, inFolder);
end
for vi = 1:length(vv)
    v_ind = find(vhold == vv(vi));
    for gi = 1:length(gg)
        if (dataMode == 1 || dataMode == 2) && ~(gi == 3 || gi == 4 || gi == 5)
            continue;
        end
        g_ind = find(gincr == gg(gi));
        for hi = 1:length(pp)
            p_ind = find(pharm == pp(hi));
            if dataMode == 0
                vgp_ind = intersect(intersect(v_ind, g_ind), p_ind);
            else
                vgp_ind = intersect(intersect(intersect(v_ind, g_ind), p_ind), indtofit);
            end
            % Find/calculate LTS/burst statistics for each group
            [all_stats, mean_stats, std_stats, ct_stats, err_stats, highbar_stats, lowbar_stats] = ...
                m3ha_compute_ltsburst_statistics(vgp_ind, cellID, ltspeaktime, ...
                            spikesperpeak, bursttime, spikesperburst);
            for bi = 1:length(measureTitle)
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
                    if dataMode == 0
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
all_stats_gp = cell(1, length(measureTitle));
mean_stats_gp = cell(1, length(measureTitle));
std_stats_gp = cell(1, length(measureTitle));
ct_stats_gp = cell(1, length(measureTitle));
err_stats_gp = cell(1, length(measureTitle));
highbar_stats_gp = cell(1, length(measureTitle));
lowbar_stats_gp = cell(1, length(measureTitle));
for bi = 1:length(measureTitle)
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
    if (dataMode == 1 || dataMode == 2) && ~(gi == 3 || gi == 4 || gi == 5)
        continue;
    end
    g_ind = find(gincr == gg(gi));
    for hi = 1:length(pp)
        p_ind = find(pharm == pp(hi));
        if dataMode == 0
            gp_ind = intersect(g_ind, p_ind);
        else
            gp_ind = intersect(intersect(g_ind, p_ind), indtofit);
        end
        % Find/calculate LTS/burst statistics for each group
        [all_stats, mean_stats, std_stats, ct_stats, err_stats, highbar_stats, lowbar_stats] = ...
            m3ha_compute_ltsburst_statistics(gp_ind, cellID, ltspeaktime, spikesperpeak, bursttime, spikesperburst);
        for bi = 1:length(measureTitle)
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
pvalue_vg = zeros(length(measureTitle), length(gg), length(vv));
                        % p value across pharm conditions for each vg-condition
pvalue_g = zeros(length(measureTitle), length(gg));
                        % p value across pharm conditions for each g incr
table_vg = cell(length(measureTitle), length(gg), length(vv));
table_g = cell(length(measureTitle), length(gg));
%{ 
%For Multiple comparison test (multcompare(stats))
stats_vg = cell(length(measureTitle), length(gg), length(vv));
stats_g = cell(length(measureTitle), length(gg));
%}
ct_vg = zeros(length(pp), 1);
ct_g = zeros(length(pp), 1);
for bi = 1:length(measureTitle)
    for gi = 1:length(gg)
        % If under dataMode 1 & 2, only do this for G incr = 100%, 200%, 400%
        if (dataMode == 1 || dataMode == 2) && ~(gi == 3 || gi == 4 || gi == 5)
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
        ylabel(measureTitle{bi}, 'FontSize', ylabelFontSize);
        if ~bigFontFlag
            xlabel('Pharm Condition');
            title(['All traces with G scaled at ', num2str(gg(gi)), '%']);
        end
        figname = fullfile(figFolder, [measureStr{bi}, '_', num2str(gg(gi)), 'g_boxplot', suffix, '.png']);
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
            ylabel(measureTitle{bi}, 'FontSize', ylabelFontSize);
            if ~bigFontFlag
                xlabel('Pharm Condition');
                title(['Vhold = ', num2str(vv(vi)), ' mV with G scaled at ', num2str(gg(gi)), '%']);
            end
            figname = fullfile(figFolder, [measureStr{bi}, ...
                '_', num2str(gg(gi)), 'g_v', num2str(vv(vi)), '_boxplot', suffix,'.png']);
            save_all_figtypes(h, figname, figTypes);
        end
    end
end

% Save variables as .mat file
fprintf('Saving variables ... \n');
matFile = fullfile(figFolder, ['ltsburst_statistics', suffix, '.mat']);
save(matFile, 'measureTitle', 'measureStr', 'scpgv_ind', ...
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
parfor bi = 1:length(measureTitle)
    for vi = 1:length(vv)
        h = figure('Visible', 'off');
        set(h, 'Name', [measureStr{bi}, '_', num2str(vv(vi))]);
        clf(h);
        % Plot means
        if dataMode == 0
            bar3(1:length(pp), mean_stats_vgp{bi}(:, :, vi), 0.12, 'detached'); hold on;
        elseif dataMode == 1 || dataMode == 2
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
        zlabel(measureTitle{bi});
        title(['Vm = ', vvLabel{vi}]);
        figname = fullfile(figFolder, [measureStr{bi}, '_', num2str(vv(vi)), suffix, '.png']);
        save_all_figtypes(h, figname, figTypes);
        close(h);
    end
end

%% Create 2D bar graph for burst statistics for each Vhold value with IPSC conductance amplitude scaling fixed (Figure 3.5 in Christine's thesis)
fprintf('Plotting 2D bar graphs for burst statistics for each Vhold value with IPSC conductance amplitude scaling fixed ... \n');
for gi = 1:length(gg)
    % If under dataMode 1 & 2, only do this for G incr = 100%, 200%, 400%
    if (dataMode == 1 || dataMode == 2) && ~(gi == 3 || gi == 4 || gi == 5)
        continue;
    end        

    parfor bi = 1:length(measureTitle)
        h = figure('Visible', 'off');
        set(h, 'Name', [measureStr{bi}, '_vsep_', num2str(gg(gi)), 'g']);
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
                ylabel(measureTitle{bi});
            end
            title(['Vm = ', vvLabel{vi}]);
        end
        figname = fullfile(figFolder, [measureStr{bi}, '_vsep_', num2str(gg(gi)), 'g', suffix, '.png']);
        save_all_figtypes(h, figname, figTypes);
        close(h);
    end
end

%% Create 3D bar graph for burst statistics
fprintf('Plotting 3D bar graphs for all burst statistics ... \n');
parfor bi = 1:length(measureTitle)
    h = figure('Visible', 'off');
    set(h, 'Name', measureStr{bi});
    clf(h);
    % Plot means
    if dataMode == 0
        bar3(1:length(pp), mean_stats_gp{bi}(:, :), 0.12, 'detached'); hold on;
    elseif dataMode == 1 || dataMode == 2
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
    zlabel(measureTitle{bi});
    title('All traces');
    figname = fullfile(figFolder, [measureStr{bi}, suffix, '.png']);
    save_all_figtypes(h, figname, figTypes);
    close(h);
end

%% Create 2D bar graph for burst statistics with IPSC conductance amplitude scaling fixed
fprintf('Plotting 2D bar graphs for all burst statistics with IPSC conductance amplitude scaling fixed ... \n');
for gi = 1:length(gg)
    % If under dataMode 1 & 2, only do this for G incr = 100%, 200%, 400%
    if (dataMode == 1 || dataMode == 2) && ~(gi == 3 || gi == 4 || gi == 5)
        continue;
    end        

    parfor bi = 1:length(measureTitle)
        % Count the number of cells in each pharm group
        ct_g = cellfun(@length, all_stats_gp{bi}(:, gi));

        % Create a 2D bar graph
        h = figure('Visible', 'off');
        set(h, 'Name', [measureStr{bi}, '_', num2str(gg(gi)), 'g']);
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
        ylabel(measureTitle{bi});
        title(['All traces with G scaled at ', num2str(gg(gi)), '%']);
        figname = fullfile(figFolder, [measureStr{bi}, '_', num2str(gg(gi)), 'g', suffix, '.png']);
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