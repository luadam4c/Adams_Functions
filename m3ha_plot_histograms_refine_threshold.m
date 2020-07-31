function m3ha_plot_histograms_refine_threshold (fitmode, infolder, outfolder, groupmode)
%% Plot histograms for sweep information & passive fit results that will be used for fitting
% Usage: m3ha_plot_histograms_refine_threshold (fitmode, infolder, outfolder, groupmode)
% Arguments: 
%       fitmode     - 0 - all data
%                   - 1 - all of g incr = 100%, 200%, 400%
%                   - 2 - all of g incr = 100%, 200%, 400% 
%                    but exclude cell-pharm-g_incr sets containing problematic sweeps
%       infolder    - (opt) the directory that contains the matfile to read
%                must be a directory
%                default == //media/adamX/m3ha/data_dclamp/take4/
%        outfolder    - (opt) the directory to output histograms 
%                    (different subdirectories will be created for each fitmode)
%                must be a directory
%                default == //media/adamX/m3ha/data_dclamp/take4/
%        groupmode    must be one of the following:
%                'cell_actVhold'
%                'cell'
%                default == 'cell'
% 
% Requires:
%       "infolder"/dclampdatalog_take4.mat
%       "infolder"/trace_comparison.mat
%       cd/m3ha_find_ind_to_fit.m
%       cd/m3ha_locate_homedir.m
%       cd/m3ha_specs_for_datamode.m
%       cd/fit_gaussians_and_refine_threshold.m
%       cd/plot_and_save_histogram.m
%       cd/plot_and_save_boxplot.m
%       cd/find_in_strings.m
%       cd/structs2vecs.m
%
% Used by:
%       cd/m3ha_parse_dclamp_data.m
%

% File History:
% 2016-09-05 - Created
% 2016-09-13 - Added fitmode
% 2016-09-13 - Added max_numComponents
% 2016-09-14 - Added maxnoise & peakclass
% 2016-09-15 - Changed name of function to PlotHistograms.m and added fitmode == 0
% 2016-10-13 - Changed the way vectors are imported so that it takes an arbitrary number of sweep info vectors
%        Now needs find_in_strings.m from /home/Matlab/Adams_Functions
% 2016-10-14 - Added suffix; changed the directory name for fitmode == 0 to include suffix ‘_all’
% 2016-10-15 - Made infolder and outfolder optional arguments
% 2016-10-20 - Removed the two correlation plots now that PlotCorrelations.m is ready
% 2016-10-31 - Placed suffix into specs_for_fitmode.m
% 2016-11-10 - Expanded to use multiple files, added passive params
% 2016-12-01 - Added groupmode, made default to be grouping by cell only
% 2016-12-01 - Changed number of bins of passive parameter histograms to 10
% 2016-12-04 - Changed logvariables to logvariables_params for mpassive
% 2016-12-04 - Added dclampPassiveLog_byswps but plot only logvariables_swpinfo
% 2016-12-05 - Added the bar graphs tau0_tau1_c and tau0_tau1_rising_c
% 2016-12-05 - Moved part of the code to structs2vecs.m
% 2016-12-08 - Changed the grouping for variables in trace_comparison.mat & dclampPassiveLog_byswps_all.mat to be by cell ID#
% 2016-12-08 - Now plots boxplots as well for variables in trace_comparison.mat & dclampPassiveLog_byswps_all.mat
% 2017-02-09 - Plots gaussian fit for RMSE rising and falling phases
% 2017-04-05 - Changed peakclass_labels to include overrules
% 2017-04-05 - Accounted for the fact that vectors in dclampPassiveLog_XXX_tofit.mat are already restricted
% 2017-05-05 - BT - Added additional RMSE analysis plots

%% Flags
debugflag = 1;

%% Specify which matfiles to use; assumed to be in infolder. 
% Files with sweep properties information
% Each must have logheader and logvariables, which must satisfy:
% (1) must be row cell arrays with same length
% (2) the first item is fnrow
swpinfo_files = {'dclampdatalog_take4.mat', 'trace_comparison.mat', 'dclampPassiveLog_byswps_suffix.mat'};
passive_file_bysets = 'dclampPassiveLog_bysets_suffix.mat';
passive_file_bycells = 'dclampPassiveLog_bycells_suffix.mat';

%% Maximum number of components in the Gaussian mixture fit
max_numComponents_1 = 3;    % for narrowest 2nd derivative
max_numComponents_2 = 15;   % for RMSE rising and falling phases

%% Parameters used for data reorganization
lts_thr = -0.0023;        % 2nd derivative in kV/s^2 below which defines an LTS peak
lts_thr_alt = -0.0081823; % Alternative threshold that determines "gray area"

%% Peak Classification (corresponds to peakclass == 1 ~ 9):
peakclass_labels = {'Not LTS by prominence', ...
            'Not LTS by narrowness', ...
            'Not LTS by shape', ...
            'Not LTS by overrule', ...
            'LTS by overrule', ...
            'LTS with no burst; contentious', ...
            'LTS with burst; contentious', ...
            'LTS with no burst; definite', ...
            'LTS with burst; definite'};

%% Passive fit grouping
VholdBC = [-62.5; -67.5; -72.5];            % Possible holding potential bin centers (LJP-corrected)
VholdBC_labels = {'-62.5 mV', '-67.5 mV', '-72.5 mV'};    % for legend

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 1
    error('A fitmode is required, type ''help m3ha_plot_histograms_refine_threshold'' for usage');
elseif isempty(fitmode) || ~isnumeric(fitmode) || ~(fitmode == 0 || fitmode == 1 || fitmode == 2)
    error('fitmode out of range!, type ''help m3ha_plot_histograms_refine_threshold'' for usage');
elseif nargin >= 2 && ~isdir(infolder)
    error('infolder must be a directory!');
elseif nargin >= 3 && ~isdir(outfolder)
    error('outfolder must be a directory!');
elseif nargin >= 4 && ~ischar(groupmode)
    error('groupmode must be either ''cell_actVhold'' or ''cell''!')
elseif nargin >= 4 && ischar(groupmode) && (~strcmp(groupmode, 'cell_actVhold') || ~strcmp(groupmode, 'cell'))
    error('groupmode must be either ''cell_actVhold'' or ''cell''!')
end

%% Locate home directory
homeDirectory = m3ha_locate_homedir;

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
    functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsdirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsdirectory does not exist!');
end
addpath_custom(fullfile(functionsdirectory, '/Adams_Functions/'));        % for %%%
addpath_custom(fullfile(functionsdirectory, '/Downloaded_Functions/'));    % for %%%
addpath_custom(fullfile(functionsdirectory, '/Brians_Functions/'));        % for reorder.m

%% Set defaults for optional arguments
if nargin < 2
    if debugflag
        infolder = fullfile(homedirectory, '/data_dclamp/take4/debug/');
    else
        infolder = fullfile(homedirectory, '/data_dclamp/take4/');
    end
end
if nargin < 3
    if debugflag
        outfolder = fullfile(homedirectory, '/data_dclamp/take4/debug/');
%        outfolder = '/media/shareX/share/Brian/Analyze_fitting_errors/debug/';
    else
        outfolder = fullfile(homedirectory, '/data_dclamp/take4/');
%        outfolder = '/media/shareX/share/Brian/Analyze_fitting_errors/';
    end
end
if nargin < 4
    groupmode = 'cell';
end

% Set suffix and title modification according to fitmode
[suffix, titleMod] = m3ha_specs_for_datamode(fitmode);

%% Choose passive file
if strcmp(groupmode, 'cell_actVhold')
    passive_file = passive_file_bysets;
elseif strcmp(groupmode, 'cell')
    passive_file = passive_file_bycells;
end

%% Add suffixes to files
swpinfo_files{3} = strrep(swpinfo_files{3}, '_suffix', suffix);
passive_file = strrep(passive_file, '_suffix', suffix);

%% Find path of matfile to use
nfiles = numel(swpinfo_files);
fullmatfilepaths = cell(1, nfiles);
m = cell(1, nfiles);
for f = 1:nfiles
    fullmatfilepaths{f} = fullfile(infolder, swpinfo_files{f});
    if exist(fullmatfilepaths{f}, 'file') ~= 2
        error('The file %s doesn''t exist!', fullmatfilepaths{f});
    else
        m{f} = matfile(fullmatfilepaths{f});
    end
end
passive_path = fullfile(infolder, passive_file);
if exist(passive_path, 'file') ~= 2
    error('The file %s doesn''t exist!', passive_path);
else
    mpassive = matfile(passive_path);
end

%% Print user specifications
fprintf('Using fit mode == %d ... \n', fitmode);
for f = 1:nfiles
    fprintf('Using matfile == %s ... \n', fullmatfilepaths{f});    
end
fprintf('Using max_numComponents_1 == %d ... \n', max_numComponents_1);
fprintf('Using max_numComponents_2 == %d ... \n', max_numComponents_2);
fprintf('Using lts_thr == %g ... \n', lts_thr);
fprintf('Using lts_thr_alt == %g ... \n', lts_thr_alt);

%% Create folders for saving files
outfolder_hist = fullfile(outfolder, ['/histograms', suffix, '/']);
outfolder_box = fullfile(outfolder, ['/boxplots', suffix, '/']);
outfolder_bar = fullfile(outfolder, ['/bargraphs', suffix, '/']);
outfolders = {outfolder_hist, outfolder_box, outfolder_bar};
for f = 1:numel(outfolders)
    if exist(outfolders{f}, 'dir') ~= 7
        mkdir(outfolders{f});
        fprintf('Made directory %s \n', outfolders{f});
    end
end

%% Import data
% All sweep info with be plotted except data filename (fnrow)
hist_labels = [];
hist_var_names = [];
nvars = zeros(1, nfiles);                % number of variables in each file
for f = 1:nfiles
    if f == 3                % dclampPassiveLog_byswps
        logheader = m{f}.logheader_swpinfo;
        logvariables = m{f}.logvariables_swpinfo;
    else
        logheader = m{f}.logheader;
        logvariables = m{f}.logvariables;        
    end
    nvars(f) = numel(logheader) - 1;
    hist_labels = [hist_labels, logheader(1, 2:end)];             % labels are for all sweep info except data filename
    hist_var_names = [hist_var_names, logvariables(1, 2:end)];    % variable names corresponding to each label
end
nhists = numel(hist_var_names);            % number of histograms to plot
nvars_orig = nvars(1);                    % number of original variables

% Read in the sweep information from matfile
hist_vectors = cell(nhists, 1);            % stores vectors for the histograms
ct = 0;                        % count variables
for f = 1:nfiles
    for r = 1:nvars(f)
        ct = ct + 1;
        hist_vectors{ct} = m{f}.(hist_var_names{ct});
    end
end

%% Import passive fit data
passive_labels = mpassive.logheader_params;
passive_var_names = mpassive.logvariables_params;
num_passive = numel(passive_var_names);    % number of passive parameters
if strcmp(groupmode, 'cell_actVhold')
    celln_set = mpassive.celln_set;
    VholdBC_set = mpassive.VholdBC_set;
    params_set = mpassive.params_set;

    % Extract vectors of passive parameters from params_set
    passive_vectors = structs2vecs(params_set);

    % Convert VholdBC_set values to 1, 2 or 3
    VholdBCn_set = zeros(size(VholdBC_set));
    for t = 1:numel(VholdBC_labels)
        VholdBCn_set(VholdBC_set == VholdBC(t)) = t;
    end

    % For parfor to work
    celln_cell = [];

elseif strcmp(groupmode, 'cell')
    celln_cell = mpassive.celln_cell;
    params_cell = mpassive.params_cell;
    params_S_R2_cell = mpassive.params_S_R2_cell;

    % Extract vectors of passive parameters from params_cell
    passive_vectors = structs2vecs(params_cell);
    passive_S_R2_vectors = structs2vecs(params_S_R2_cell);

    % For parfor to work
    celln_set = [];
    VholdBC_set = [];
    VholdBCn_set = [];
end

%% Restrict data to fit if fitmode > 0
if fitmode > 0
    %% Find indices of fnrow in dclampdatalog_take4.mat that will be used for fitting
    fnrow_old = m{1}.fnrow;                % file names for each sweep, used for indtofit
    cellidrow_old = m{1}.cellidrow;        % Cell ID # for each sweep, used for indtofit
    prow_old = m{1}.prow;                % Pharm condition for each sweep, used for indtofit
    grow_old = m{1}.grow;                % G incr for each sweep, used for indtofit
    indtofit = m3ha_find_ind_to_fit(fnrow_old, cellidrow_old, prow_old, ...
                                    grow_old, fitmode, infolder);

    %% Update variables for histograms
    for k = 1:sum(nvars(1:2))            % vectors in dclampPassiveLog_XXX_tofit.mat are already restricted
        hist_vectors{k} = hist_vectors{k}(indtofit);    % Restrict sweep information to those of interest for fitting
    end
end

%% Create file names
hist_filenames = cell(1, nhists);        % stores file names for the histograms
for k = 1:nhists
    hist_filenames{k} = [hist_var_names{k}, suffix, '.png'];
end
box_filenames = cell(1, nhists);        % stores file names for the boxplots
for k = 1:nhists
    box_filenames{k} = [hist_var_names{k}, '_c', suffix, '.png'];
end
passive_hist_filenames = cell(1, nhists);    % stores file names for the histograms
for k = 1:num_passive
    passive_hist_filenames{k} = [passive_var_names{k}, suffix, '.png'];
    if strcmp(groupmode, 'cell_actVhold')
        passive_bar_filenames{k} = [passive_var_names{k}, '_cv', suffix, '.png'];
    elseif strcmp(groupmode, 'cell')
        passive_bar_filenames{k} = [passive_var_names{k}, '_c', suffix, '.png'];
    end
end

%% Plot and save histograms
% Extract vectors needed
cellidrow = hist_vectors{find_in_strings('cellidrow', hist_var_names)};
peakclass = hist_vectors{find_in_strings('peakclass', hist_var_names)};

% Create labels for cell ID
total_ncells = max(cellidrow);
cellidrow_labels = cell(1, total_ncells);
for c = 1:total_ncells
    cellidrow_labels{c} = ['Cell #', num2str(c)];
end

% Plot and save histograms
parfor k = 1:nhists
    h1 = figure(5000 + k);
    set(h1, 'Visible', 'off');
    clf(h1);
    h2 = figure(6000 + k);
    set(h2, 'Visible', 'off');
    clf(h2);
    if k <= nvars(1)        % variables in dclampdatalog_take4.mat
        plot_and_save_histogram(h1, hist_vectors{k}, hist_labels{k}, '# of sweeps', ...
                    outfolder_hist, hist_filenames{k}, 50, ...
                    titleMod, peakclass, peakclass_labels);
    else                    % variables in trace_comparison.mat & dclampPassiveLog_byswps_all.mat
        plot_and_save_histogram(h1, hist_vectors{k}, hist_labels{k}, '# of sweeps', ...
                    outfolder_hist, hist_filenames{k}, 50, ...
                    titleMod, cellidrow, cellidrow_labels);
        plot_and_save_boxplot(h2, hist_vectors{k}, hist_labels{k}, ...
                    outfolder_box, box_filenames{k}, ...
                    titleMod, cellidrow, 'Cell ID #');
    end
    close all
end
parfor k = 1:num_passive
    h = figure(10000 + k);
    set(h, 'Visible', 'off');
    clf(h);
    if strcmp(groupmode, 'cell_actVhold')
        plot_and_save_histogram(h, passive_vectors{k}, passive_labels{k}, '# of sets', ...
                    outfolder_hist, passive_hist_filenames{k}, 50, ...
                    titleMod, VholdBCn_set, VholdBC_labels);
    elseif strcmp(groupmode, 'cell')
        plot_and_save_histogram(h, passive_vectors{k}, passive_labels{k}, '# of cells', ...
                    outfolder_hist, passive_hist_filenames{k}, 10, ...
                    titleMod);
    end

    h = figure(15000 + k);
    set(h, 'Visible', 'off');
    clf(h);
    if strcmp(groupmode, 'cell_actVhold')
        gscatter(celln_set, passive_vectors{k}, VholdBC_set, [], 'o');
    elseif strcmp(groupmode, 'cell')
        plot(celln_cell, passive_vectors{k}, 'o');
    end
    xlabel('Cell number')
    ylabel(passive_labels{k})
    title(sprintf('All values for %s %s', passive_labels{k}, titleMod));
    figname = fullfile(outfolder_bar, passive_bar_filenames{k});
    saveas(h, figname);
    close all
end

%% tau0 & tau1
if strcmp(groupmode, 'cell')
    tau0 = passive_vectors{find_in_strings('tau0', passive_var_names)};
    tau1 = passive_vectors{find_in_strings('tau1', passive_var_names)};
    tau0_S_R2 = passive_S_R2_vectors{find_in_strings('tau0', passive_var_names)};
    tau1_S_R2 = passive_S_R2_vectors{find_in_strings('tau1', passive_var_names)};
    figname = fullfile(outfolder_bar, ['tau0_tau1_c', suffix, '.png']);
    plot_tau0_tau1(celln_cell, tau0, tau1, figname, titleMod, suffix);
    figname = fullfile(outfolder_bar, ['tau0_tau1_rising_c', suffix, '.png']);
    plot_tau0_tau1(celln_cell, tau0_S_R2, tau1_S_R2, figname, titleMod, suffix);
end

%% Refine threshold to use to define an LTS
% The function find_in_strings.m is in /home/Matlab/Adams_Functions/
fprintf('Finding new LTS threshold ...\n');

% Extract vectors needed
bursttime = hist_vectors{find_in_strings('bursttime', hist_var_names)};
peak2ndder = hist_vectors{find_in_strings('peak2ndder', hist_var_names)};
rmse_R_row = hist_vectors{find_in_strings('rmse_R_row', hist_var_names)};
rmse_F_row = hist_vectors{find_in_strings('rmse_F_row', hist_var_names)};

% Find sweeps with bursts and overall burst threshold
burst_thr = max(peak2ndder(bursttime > 0));
fprintf('Possible 2nd derivative threshold for bursts is %g\n', burst_thr);

% Create file name for figures
figname = ['Narrowest_peak_2ndder_nospont', suffix, '.png'];
peakclass_labels{3} = 'Not LTS by shape (taken out)';
[BestModel, numComponents, Mu_best, Std_best, Prop_best, minAIC, lts_thr_new, lts_thr_alt_new, new_peakclass, min_threshold, ~] = ...
    fit_gaussians_and_refine_threshold(peak2ndder', figname, ...
       'Narrowest Peak 2nd Derivatives', 'MaxNumComponents', max_numComponents_1, ...
       'OldThr', lts_thr, 'ThrMin', burst_thr, 'OutFolder', outfolder_hist, ...
       'FitMode', fitmode, 'PeakClass', peakclass', 'PeakClassLabels', peakclass_labels, ...
       'MinThreshold', 0, 'ThresMode', 'minFirstTwoComponents');
fprintf('New LTS threshold is %g kV/s^2\n', lts_thr_new);
fprintf('Alternate LTS threshold is %g kV/s^2\n', lts_thr_alt_new);

% Change legend position in fit_gaussians_and_refine_threshold 
% Following two plots are RMSE with fit
% Create file name for figures
figname = ['rmse_R_row_Fit', suffix, '.png'];
[BestModelR, numComponentsR, Mu_bestR, Std_bestR, Prop_bestR, minAICR, lts_thr_newR, lts_thr_alt_newR, new_peakclassR, min_thresholdR, ~] = ...
    fit_gaussians_and_refine_threshold(rmse_R_row', figname, ...
       'RMSE (mV) in the rising phase', 'MaxNumComponents', max_numComponents_2, ...
       'OldThr', 0, 'ThrMin', 0, 'OutFolder', outfolder_hist, ...
       'FitMode', fitmode, 'PeakClass', cellidrow', 'PeakClassLabels', cellidrow_labels, ...
       'MinThreshold', 0, 'ThresMode', 'minFirstTwoComponents');
fprintf('New LTS threshold is %g kV/s^2\n', lts_thr_newR);
fprintf('Alternate LTS threshold is %g kV/s^2\n', lts_thr_alt_newR);

figname = ['rmse_F_row_Fit', suffix, '.png'];
[BestModelF, numComponentsF, Mu_bestF, Std_bestF, Prop_bestF, minAICF, lts_thr_newF, lts_thr_alt_newF, new_peakclassF, min_thresholdF, ~] = ...
    fit_gaussians_and_refine_threshold(rmse_F_row', figname, ...
       'RMSE (mV) in the falling phase', 'MaxNumComponents', max_numComponents_2, ...
       'OldThr', 0, 'ThrMin', 0, 'OutFolder', outfolder_hist, ...
       'FitMode', fitmode, 'PeakClass', new_peakclassR', 'PeakClassLabels', cellidrow_labels, ...
       'MinThreshold', min_thresholdR, 'ThresMode', 'minFirstTwoComponents');
fprintf('New LTS threshold is %g kV/s^2\n', lts_thr_newF);
fprintf('Alternate LTS threshold is %g kV/s^2\n', lts_thr_alt_newF);

% Following two plots are RMSE examining above and below threshold
% Create file name for figures
figname = ['rmse_R_row_Fit_threshold', suffix, '.png'];
[BestModelR, numComponentsR, Mu_bestR, Std_bestR, Prop_bestR, minAICR, lts_thr_newR, lts_thr_alt_newR, new_peakclassR, min_thresholdR, ~] = ...
    fit_gaussians_and_refine_threshold(rmse_R_row', figname, ...
       'RMSE (mV) in the rising phase', 'MaxNumComponents', max_numComponents_2, ...
       'OldThr', 0, 'ThrMin', 0, 'OutFolder', outfolder_hist, ...
       'FitMode', fitmode, 'PeakClass', cellidrow', 'PeakClassLabels', cellidrow_labels, ...
       'MinThreshold', 0, 'ThresMode', 'minFirstTwoComponents');
fprintf('New LTS threshold is %g kV/s^2\n', lts_thr_newR);
fprintf('Alternate LTS threshold is %g kV/s^2\n', lts_thr_alt_newR);

% Create file name for figures
figname = ['rmse_F_row_Fit_threshold', suffix, '.png'];
[BestModelF, numComponentsF, Mu_bestF, Std_bestF, Prop_bestF, minAICF, lts_thr_newF, lts_thr_alt_newF, new_peakclassF, min_thresholdF, ~] = ...
    fit_gaussians_and_refine_threshold(rmse_F_row', figname, ...
       'RMSE (mV) in the falling phase', 'MaxNumComponents', max_numComponents_2, ...
       'OldThr', 0, 'ThrMin', 0, 'OutFolder', outfolder_hist, ...
       'FitMode', fitmode, 'PeakClass', new_peakclassR', 'PeakClassLabels', cellidrow_labels, ...
       'MinThreshold', 0, 'ThresMode', 'minFirstTwoComponents');
fprintf('New LTS threshold is %g kV/s^2\n', lts_thr_newF);
fprintf('Alternate LTS threshold is %g kV/s^2\n', lts_thr_alt_newF);

% Following two plots are RMSE with cells 5, 14, and 19 being examined (change cells in fit_gaussians_and_refine_threshold.m)
% Create file name for figures
figname = ['rmse_R_row_Fit_traces', suffix, '.png'];
[BestModelR, numComponentsR, Mu_bestR, Std_bestR, Prop_bestR, minAICR, lts_thr_newR, lts_thr_alt_newR, new_peakclassR, min_thresholdR, ~] = ...
    fit_gaussians_and_refine_threshold(rmse_R_row', figname, ...
       'RMSE (mV) in the rising phase', 'MaxNumComponents', max_numComponents_2, ...
       'OldThr', 0, 'ThrMin', 0, 'OutFolder', outfolder_hist, ...
       'FitMode', fitmode, 'PeakClass', cellidrow', 'PeakClassLabels', cellidrow_labels, ...
       'MinThreshold', 0, 'ThresMode', 'minFirstTwoComponents');
fprintf('New LTS threshold is %g kV/s^2\n', lts_thr_newR);
fprintf('Alternate LTS threshold is %g kV/s^2\n', lts_thr_alt_newR);

% Create file name for figures
figname = ['rmse_F_row_Fit_traces', suffix, '.png'];
[BestModelF, numComponentsF, Mu_bestF, Std_bestF, Prop_bestF, minAICF, lts_thr_newF, lts_thr_alt_newF, new_peakclassF, min_thresholdF, ~] = ...
    fit_gaussians_and_refine_threshold(rmse_F_row', figname, ...
       'RMSE (mV) in the falling phase', 'MaxNumComponents', max_numComponents_2, ...
       'OldThr', 0, 'ThrMin', 0, 'OutFolder', outfolder_hist, ...
       'FitMode', fitmode, 'PeakClass', cellidrow', 'PeakClassLabels', cellidrow_labels, ...
       'MinThreshold', 0, 'ThresMode', 'minFirstTwoComponents');
%        fitmode, new_peakclassR', cellidrow_labels, 0);
fprintf('New LTS threshold is %g kV/s^2\n', lts_thr_newF);
fprintf('Alternate LTS threshold is %g kV/s^2\n', lts_thr_alt_newF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_tau0_tau1(celln_cell, tau0, tau1, figname, titleMod, suffix)

tau0_min = min(tau0);

h = figure(floor((1+rand())*10^6));
%set(h, 'Visible', 'off');
clf(h);
plot(celln_cell, tau0, 'ro', 'Displayname', 'tau0'); hold on;
plot(celln_cell, tau1, 'bd', 'Displayname', 'tau1');
line([0 50], [tau0_min tau0_min], 'Color', 'g', 'LineStyle', '--', ...
    'Displayname', ['Minimum of tau0 == ', num2str(tau0_min), ' (ms)']);
xlabel('Cell number')
ylabel('Time constant (ms)')
legend('location', 'northwest');
title(sprintf('All values for Time constant (ms) %s', titleMod));
saveas(h, figname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

% THIS DOESN'T WORK!
% Construct file names for the histograms
hist_filenames = {};
for f = 1:numel(hist_var)
    varname = @(x) inputname(1);        % Creates an anonymous function that prints the input argument as a string
    hist_filenames{f} = varname(hist_var{f});    % Returns an empty string
end

file_info = whos(m);        % information about the matfile
ct = 0;
for k = 1:numel(file_info)
    if strcmp(file_info(k).name, 'logheader')
        hist_label = m.logheader(1, 2:end);    % Labels are for all sweep info except data filename
    elseif strcmp(file_info(k).name, 'fnrow')
        fnrow = m.fnrow;            % File name for each sweep
    elseif isnumeric(m.(file_info(k).name))        % must be an array of doubles
        ct = ct + 1;
        swpinfo{ct, 1} = m.(file_info(k).name);    % Read in the sweep information from matfile
        swpinfo{ct, 2} = file_info(k).name;    % Store the variable name for that sweep information
    end
end
nswpinfo = ct;                        % number of sweep information vectors to plot

for k = 1:nswpinfo
    if fitmode == 0        % all data
        hist_filenames{k, 1} = [swpinfo{k, 2}, '_all.png'];
    elseif fitmode == 1    % all of g incr = 100%, 200%, 400%
        hist_filenames{k, 1} = [swpinfo{k, 2}, '_100-400all.png.png'];
    elseif fitmode == 2    % all of g incr = 100%, 200%, 400% 
                % but exclude cell-pharm-g_incr sets containing problematic sweeps
        hist_filenames{k, 1} = [swpinfo{k, 2}, '_tofit.png'];
    end
end

fnrow = m.fnrow;
cellidrow = m.cellidrow;
prow = m.prow;
vrow = m.vrow;
grow = m.grow;
swpnrow = m.swpnrow;
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
peak2ndder = m.peak2ndder;
peakprom = m.peakprom;
peakclass = m.peakclass;
spikesperpeak = m.spikesperpeak;
ltspeaktime = m.ltspeaktime;
ltspeakval = m.ltspeakval;
bursttime = m.bursttime;
spikesperburst = m.spikesperburst;

numswps = numel(swpinfo{1});    % number of sweeps

%load(fullfile(data_dir, 'burst_statistics.mat'), 'vgpcs_ind');

    fnrow = fnrow(indtofit);
    cellidrow = cellidrow(indtofit);
    prow = prow(indtofit);
    vrow = vrow(indtofit);
    grow = grow(indtofit);
    swpnrow = swpnrow(indtofit);
    gabab_amp = gabab_amp(indtofit);
    gabab_Trise = gabab_Trise(indtofit);
    gabab_TfallFast = gabab_TfallFast(indtofit);
    gabab_TfallSlow = gabab_TfallSlow(indtofit);
    gabab_w = gabab_w(indtofit);
    currpulse = currpulse(indtofit);
    Rin = Rin(indtofit);
    ioffset_old = ioffset_old(indtofit);
    imint = imint(indtofit);
    imin = imin(indtofit);
    actVhold = actVhold(indtofit);
    maxnoise = maxnoise(indtofit);
    narrowpeaktime = narrowpeaktime(indtofit);
    peak2ndder = peak2ndder(indtofit);
    peakprom = peakprom(indtofit);
    peakclass = peakclass(indtofit);
    spikesperpeak = spikesperpeak(indtofit);
    ltspeaktime = ltspeaktime(indtofit);
    ltspeakval = ltspeakval(indtofit);
    bursttime = bursttime(indtofit);
    spikesperburst = spikesperburst(indtofit);

hist_label = {'Data filename', 'Cell ID #', ...
        'Pharm condition', 'Vhold (mV)', 'GABAB IPSC G incr (%)', 'Within condition sweep #', ...
        'GABAB IPSC G amp (nS)', 'GABAB IPSC G Trise (ms)', 'GABAB IPSC G TfallFast (ms)', ...
        'GABAB IPSC G TfallSlow (ms)', 'GABAB IPSC G w', ...
        'Current pulse amplitude (pA)', 'Rinput (MOhm)', ...
        'IPSC offset (obsolete) (ms)', 'IPSC peak time (ms)', 'IPSC peak amplitude (pA)', ...
        'Actual Vhold (mV)', 'Maximum noise (mV)', ...
        'Narrowest peak time (ms)', 'Narrowest peak 2nd derivative (kV/s^2)', ...
        'Peak prominence (mV)', 'Peak class #', 'Spikes per peak', ...
        'LTS onset time (ms)', 'LTS amplitude (mV)', ...
        'Burst onset time (ms)', 'Spikes per burst'};
hist_var = {fnrow, cellidrow, ...
        prow, vrow, grow, swpnrow, ...
        gabab_amp, gabab_Trise, gabab_TfallFast, ...
        gabab_TfallSlow, gabab_w, ...
        currpulse, Rin, ioffset_old, imint, imin, ...
        actVhold, maxnoise, narrowpeaktime, peak2ndder, ...
        peakprom, peakclass, spikesperpeak, ...
        ltspeaktime, ltspeakval, bursttime, spikesperburst};

if fitmode == 0        % all data
    hist_filenames = {'', 'cellidrow_all.png', ...
            'prow_all.png', 'vrow_all.png', 'grow_all.png', 'swpnrow_all.png', ...
            'gabab_amp_all.png', 'gabab_Trise_all.png', 'gabab_TfallFast_all.png', ...
            'gabab_TfallSlow_all.png', 'gabab_w_all.png', ...
            'currpulse_all.png', 'Rinput_all.png', ...
            'IPSCoffset_old_all.png', 'IPSCpeaktime_all.png', 'IPSCpeakamp_all.png', ...
            'actVhold_all.png', 'maxnoise_all.png', ...
            'Narrowest_peaktime_all.png', 'Narrowest_peak_2ndder_all.png', ...
            'peakprom_all.png', 'peakclass_all.png', 'spikesperpeak_all.png'...
            'LTSpeaktime_all.png', 'LTSpeakvalue_all.png', ...
            'bursttime_all.png', 'spikesperburst_all.png'};
elseif fitmode == 1    % all of g incr = 100%, 200%, 400%
    hist_filenames = {'', 'cellidrow_100-400all.png', ...
            'prow_100-400all.png', 'vrow_100-400all.png', 'grow_100-400all.png', 'swpnrow_100-400all.png', ...
            'gabab_amp_100-400all.png', 'gabab_Trise_100-400all.png', 'gabab_TfallFast_100-400all.png', ...
            'gabab_TfallSlow_100-400all.png', 'gabab_w_100-400all.png', ...
            'currpulse_100-400all.png', 'Rinput_100-400all.png', ...
            'IPSCoffset_old_100-400all.png', 'IPSCpeaktime_100-400all.png', 'IPSCpeakamp_100-400all.png', ...
            'actVhold_100-400all.png', 'maxnoise_100-400all.png', ...
            'Narrowest_peaktime_100-400all.png', 'Narrowest_peak_2ndder_100-400all.png', ...
            'peakprom_100-400all.png', 'peakclass_100-400all.png', 'spikesperpeak_100-400all.png'...
            'LTSpeaktime_100-400all.png', 'LTSpeakvalue_100-400all.png', ...
            'bursttime_100-400all.png', 'spikesperburst_100-400all.png'};
elseif fitmode == 2    % all of g incr = 100%, 200%, 400% but exclude cell-pharm-g_incr sets containing problematic sweeps
    hist_filenames = {'', 'cellidrow_tofit.png', ...
            'prow_tofit.png', 'vrow_tofit.png', 'grow_tofit.png', 'swpnrow_tofit.png', ...
            'gabab_amp_tofit.png', 'gabab_Trise_tofit.png', 'gabab_TfallFast_tofit.png', ...
            'gabab_TfallSlow_tofit.png', 'gabab_w_tofit.png', ...
            'currpulse_tofit.png', 'Rinput_tofit.png', ...
            'IPSCoffset_old_tofit.png', 'IPSCpeaktime_tofit.png', 'IPSCpeakamp_tofit.png', ...
            'actVhold_tofit.png', 'maxnoise_tofit.png', ...
            'Narrowest_peaktime_tofit.png', 'Narrowest_peak_2ndder_tofit.png', ...
            'peakprom_tofit.png', 'peakclass_tofit.png', 'spikesperpeak_tofit.png'...
            'LTSpeaktime_tofit.png', 'LTSpeakvalue_tofit.png', ...
            'bursttime_tofit.png', 'spikesperburst_tofit.png'};
end

parfor f = 2:numel(hist_var)
    h = figure('Visible', 'off');
    plot_and_save_histogram(h, hist_var{f}, hist_label{f}, outfolder_hist, hist_filenames{f}, 50, ...
                fitmode, peakclass, peakclass_labels);
    close(h);
end

goodb_ind = find(bursttime > 0);
burst_thr = max(peak2ndder(goodb_ind));

%% Plot input resistance versus current pulse amplitude
% Extract vectors needed
currpulse = hist_vectors{find_in_strings('currpulse', hist_var_names)};
Rin = hist_vectors{find_in_strings('Rin', hist_var_names)};
h = figure('Visible', 'off');
plot(currpulse, Rin, 'o');
xlim([-300 0]);
xlabel('Current pulse amplitude (pA)');
ylabel('Rinput (MOhm)');
if fitmode == 0
    title('Input resistance variation (all)')
elseif fitmode == 1
    title('Input resistance variation (100%, 200%, 400% g incr)')
elseif fitmode == 2
    title('Input resistance variation (for fitting)')
end
figname = fullfile(outfolder_hist, ['Rin_vs_cpa', suffix, '.png']);
saveas(h, figname);
close(h);

%% Plot input resistance versus holding membrane potential
% Extract vectors needed
actVhold = hist_vectors{find_in_strings('actVhold', hist_var_names)};
h = figure('Visible', 'off');
plot(actVhold, Rin, 'o');
xlim([-80 -55]);
xlabel('Holding membrane potential (mV)');
ylabel('Rinput (MOhm)');
if fitmode == 0
    title('Input resistance variation (all)')
elseif fitmode == 1
    title('Input resistance variation (100%, 200%, 400% g incr)')
elseif fitmode == 2
    title('Input resistance variation (for fitting)')
end
figname = fullfile(outfolder_hist, ['Rin_vs_actVhold', suffix, '.png']);
saveas(h, figname);
close(h);

    VholdBC(VholdBC == str2double(VholdBC_labels{t})) = t;

        empty_pos = find(cellfun('isempty', hist_vectors));
        hist_vectors{empty_pos(1)} = m{f}.(hist_var_names{r});    % Read in the sweep information from matfile

passive_labels = mpassive.logheader;
passive_var_names = mpassive.logvariables;

    numsets = numel(params_set);            % number of sets
    passive_mat = zeros(num_passive, numsets);
    for h = 1:numsets
        if isstruct(params_set{h})
            % Convert the params structure into a column vector
            passive_mat(:, h) = cell2mat(struct2cell(params_set{h}));    
        end
    end
    passive_vectors = mat2cell(passive_mat, ones(1, num_passive));    

    numcells = numel(params_cell);            % number of sets
    passive_mat = zeros(num_passive, numcells);
    for h = 1:numcells
        if isstruct(params_cell{h})
            % Convert the params structure into a column vector
            passive_mat(:, h) = cell2mat(struct2cell(params_cell{h}));    
        end
        if isstruct(params_S_R2_cell{h})
            % Convert the params structure into a column vector
            passive_S_R2_mat(:, h) = cell2mat(struct2cell(params_S_R2_cell{h}));    
        end
    end
    passive_vectors = mat2cell(passive_mat, ones(1, num_passive));
    passive_S_R2_vectors = mat2cell(passive_S_R2_mat, ones(1, num_passive));
    
%%% XXX: Is this OLD CODE? I actually think this is a better place for all this BT - addressed in reorder.m, falling reliant on computations done in rising first
%%%         You should make different cellidrow & cellidrow_labels here so that you can generate both types of figures at once
% 
%cell4 = find(cellidrow == 4);
%cell5 = find(cellidrow == 5);
%cell19 = find(cellidrow == 19);
%takeout_indices = [cell4 cell5 cell19];
%cellidrow(takeout_indices) = [];
%rmse_R_row(takeout_indices) = [];
%cellidrow_labels([4, 5, 19]) = [];
%max_peakclass = max(peakclass);
%y = 1;
%for x = 1:max_peakclass-3
%    while peakclass(y) == x
%        y = y + 1;
%    end
%    while peakclass(y) ~= x + 1
%        peakclass(y) = x + 1;
%        y = y + 1;
%    end
%end
    fit_gaussians_and_refine_threshold(peak2ndder', max_numComponents_1, ...
        lts_thr, burst_thr, 'Narrowest Peak 2nd Derivatives', outfolder_hist, figname, ...
        fitmode, peakclass', peakclass_labels, 0);
    fit_gaussians_and_refine_threshold(rmse_R_row', max_numComponents_2, ...
        0, 0, 'RMSE (mV) in the rising phase', outfolder_hist, figname, ...
        fitmode, cellidrow', cellidrow_labels, 0);
    fit_gaussians_and_refine_threshold(rmse_F_row', max_numComponents_2, ...
        0, 0, 'RMSE (mV) in the falling phase', outfolder_hist, figname, ...
        fitmode, new_peakclassR', cellidrow_labels, min_thresholdR);
    fit_gaussians_and_refine_threshold(rmse_R_row', max_numComponents_2, ...
        0, 0, 'RMSE (mV) in the rising phase', outfolder_hist, figname, ...
        fitmode, cellidrow', cellidrow_labels, 0);
    fit_gaussians_and_refine_threshold(rmse_F_row', max_numComponents_2, ...
        0, 0, 'RMSE (mV) in the falling phase', outfolder_hist, figname, ...
        fitmode, new_peakclassR', cellidrow_labels, 0);
    fit_gaussians_and_refine_threshold(rmse_R_row', max_numComponents_2, ...
        0, 0, 'RMSE (mV) in the rising phase', outfolder_hist, figname, ...
        fitmode, cellidrow', cellidrow_labels, 0);
    fit_gaussians_and_refine_threshold(rmse_F_row', max_numComponents_2, ...
        0, 0, 'RMSE (mV) in the falling phase', outfolder_hist, figname, ...
        fitmode, cellidrow', cellidrow_labels, 0);

%}
