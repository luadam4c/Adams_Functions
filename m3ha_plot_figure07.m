% m3ha_plot_figure07.m
%% Plots Figure 07 for the GAT Blocker paper
%
% Requires:
%       cd/addvars_custom.m
%       cd/all_files.m
%       cd/all_subdirs.m
%       cd/apply_iteratively.m
%       cd/apply_over_cells.m
%       cd/archive_dependent_scripts.m
%       cd/argfun.m
%       cd/array_fun.m
%       cd/compute_combined_trace.m
%       cd/convert_to_char.m
%       cd/create_labels_from_numbers.m
%       cd/create_label_from_sequence.m
%       cd/create_subplots.m
%       cd/decide_on_colormap.m
%       cd/extract_common_directory.m
%       cd/extract_common_prefix.m
%       cd/extract_fileparts.m
%       cd/extract_substrings.m
%       cd/extractFrom.m
%       cd/find_matching_files.m
%       cd/force_column_cell.m
%       cd/force_matrix.m
%       cd/ismatch.m
%       cd/lower_first_char.m
%       cd/m3ha_extract_candidate_label.m
%       cd/m3ha_network_analyze_spikes.m
%       cd/m3ha_network_plot_gabab.m
%       cd/m3ha_network_plot_essential.m
%       cd/m3ha_plot_violin.m
%       cd/plot_grouped_jitter.m
%       cd/plot_scale_bar.m
%       cd/plot_tuning_curve.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/sscanf_full.m
%       cd/unique_custom.m
%       cd/update_figure_for_corel.m

% File History:
% 2020-01-30 Modified from m3ha_plot_figure05.m
% 2020-02-06 Added plot200CellExamples and plot2CellM2h
% 2020-03-10 Updated pharm labels
% 2020-04-09 Added combineActivationProfiles
% 2019-04-28 Changed timeToStabilize from 2000 to 3000

%% Hard-coded parameters
% Flags
plotIpscComparison = false; %true;
plot2CellEssential = false; %true;
plot2CellM2h = false; %true;

analyze2CellSpikes = false; %true;
plotAnalysis2Cell = false; %true;
backupPrevious2Cell = false; %true;
combine2CellPopulation = false; %true;
plot2CellViolins = false; %true;

plot200CellExamples = false; %true;

analyze200CellSpikes = false; %true;
plotAnalysis200Cell = false;
backupPrevious200Cell = false;
combineActivationProfiles = false; %true;
combine200CellPopulation = false; %true;
plot200CellViolins = false; %true;
plot200CellGroupByCellJitters = false; %true;
combineEach200CellNetwork = false; %true;
plot200CellGroupByEpasJitters = false; %true;
plot200CellCumDist = false; %true;

archiveScriptsFlag = true;

% Directories
parentDirectory = fullfile('/media', 'adamX', 'm3ha');
figure07Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure07');
figure08Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure08');
networkDirectory = fullfile(parentDirectory, 'network_model');

% exampleIterName2Cell = '20200131T1345_using_bestparams_20200126_singleneuronfitting101';  % 20200131
% exampleIterName2Cell = '20200205T1353_using_bestparams_20200203_manual_singleneuronfitting0-102_2cell_examples';
% popIterName2Cell = '20200204T1042_using_bestparams_20200203_manual_singleneuronfitting0-102_vtraub_-65_2cell_spikes';
% exampleIterName200Cell = '20200204T1239_using_bestparams_20200203_manual_singleneuronfitting0-102_200cell_spikes';
% popIterName200Cell = exampleIterName200Cell;
% rankNumsToUse = [2, 4, 5, 7, 9, 10, 12, 13, 16, 20, 21, 23, 25, 29];
% popIterName2Cell = '20200208T1230_using_bestparams_20200203_manual_singleneuronfitting0-102_2cell_spikes';
% popIterName2Cell = '20200305T2334_using_bestparams_20200203_manual_singleneuronfitting0-102_2cell_REgpas_varied';
% popIterName2Cell = '20200306T1724_using_bestparams_20200203_manual_singleneuronfitting0-102_2cell_gpas_varied';
% popIterName2Cell = '20200308T2306_using_bestparams_20200203_manual_singleneuronfitting0-102_2cell_TCepas_varied';
% popIterName200Cell = exampleIterName200Cell;
% popIterName200Cell = '20200309T1346_using_bestparams_20200203_manual_singleneuronfitting0-102_200cell_TCepas_varied';
% popIterName2Cell = '20200309T0013_using_bestparams_20200203_manual_singleneuronfitting0-102_2cell_TCepas_varied';
% popIterName2Cell = '20200311T2144_using_bestparams_20200203_manual_singleneuronfitting0-102_2cell_TCepas_varied';
% popIterName200Cell = '20200312T0130_using_bestparams_20200203_manual_singleneuronfitting0-102_200cell_TCepas_varied';
% exampleIterName200Cell = '20200208T1429_using_bestparams_20200203_manual_singleneuronfitting0-102_200cell_spikes';
% exampleSeedNumberDirName = 'seedNumber_5';      % Use seed number 5 (TCepas = -70)

exampleIterName2Cell = '20200207T1554_using_bestparams_20200203_manual_singleneuronfitting0-102_REena88_TCena88_2cell_examples';
popIterName2Cell = '20200418_using_bestparams_20200203_manual_singleneuronfitting0-102';
exampleIterName200Cell = '20200408_using_bestparams_20200203_manual_singleneuronfitting0-102';
exampleSeedNumberDirName = 'seedNumber_21';      % Use seed number 21 (TCepas = -70)
popIterName200Cell = '20200408_using_bestparams_20200203_manual_singleneuronfitting0-102';
candCellSheetName = 'candidate_cells.csv';
oscParamsSuffix = 'oscillation_params';

% Well-fitted, good 2-cell network response
rankNumsToUse = [2:4, 6, 8:11, 14, 18, 19, 21, 23];

% Files

% Analysis settings
% Should be consistent with m3ha_plot_figure03.m & m3ha_plot_figure07.m
exampleCellNames2Cell = {'D101310'; 'G101310'};
% exampleCellNames200Cell = {'D101310'; 'G101310'; 'hetero4'; 'hetero8'; 'hetero12'};
exampleCellNames200Cell = {'D101310'; 'hetero12'};

gIncr = 200;                % Original dynamic clamp gIncr value
pharmConditions = (1:4)';   % Pharmacological conditions
                            %   1 - Control
                            %   2 - GAT 1 Block
                            %   3 - GAT 3 Block
                            %   4 - Dual Block
measuresOfInterest = {'oscillationProbability'; 'meanOscPeriod2Ms'; ...
                        'meanOscIndex4'; 'meanPercentActive'; ...
                        'meanOscDurationSec'; 'meanPercentActiveTC'; ...
                        'meanHalfActiveLatencyMsTC'};
measureTitles = {'Oscillation Probability'; 'Oscillation Period (ms)'; ...
                    'Oscillatory Index'; 'Active Cells (%)'; ...
                    'Oscillation Duration (sec)'; 'Active TC Cells (%)'; ...
                    'Half Activation Time (ms)'};
measuresOfInterestJitter = {'oscPeriod2Ms'; ...
                        'oscIndex4'; 'percentActive'; ...
                        'oscDurationSec'; 'percentActiveTC'; ...
                        'halfActiveLatencyMsTC'};
measureTitlesJitter = measureTitles(2:end);

% The following must be consistent with m3ha_net.hoc
timeToStabilize = 2000;         % padded time (ms) to make sure initial value 
                                %   of simulations are stabilized

% Plot settings
ipscFigWidth = 8.5;
ipscFigHeight = 6;
xLimits2CellGabab = timeToStabilize + [800, 2800];
xLimits2CellEssential = timeToStabilize + [800, 3100];
yLimitsGabab = [-1, 15];
% yLimitsEssential = {[], [], [], [], [], [], []};
yLimitsEssential = {[-100, 100], [-100, 100], [-1, 12], [-15, 5], ...
                    [1e-10, 1], [1e-10, 1], [1e-8, 1e0]};
yTicksEssential = {[-75, 75], [-75, 75], [0, 5, 10], [-10, -5, 0], ...
                    [1e-8, 1e-2], [1e-8, 1e-2], [1e-6, 1e-2]};
yLimitsM2h = [1e-10, 1];
yTicksM2h = [1e-8, 1e-2];
essential2CellFigWidth = 9.35;
essential2CellFigHeight = 1 * 7;
m2h2CellFigWidth = 8.5;
m2h2CellFigHeight = 1;
example200CellFigWidth = 8.5;
example200CellFigHeight = 3;
pharmLabelsShort = {'{\it s}Con', '{\it s}GAT1', ...
                    '{\it s}GAT3', '{\it s}Dual'};

% epasToPlot = [];
epasToPlot = [-74; -70; -66; -62];

% candidateLabels200Cell = {'candidateIDs_32'; 'candidateIDs_2,14,32,35'; ...
%                                   'candidateIDs_2,14,20,29-30,32,35-36'};
candidateLabels200Cell = {};

% candidateLabels200CellEcdfs = {'candidateIDs_2,14,32,35'; 'candidateIDs_2'; ...
%                 'candidateIDs_14'; 'candidateIDs_32'; 'candidateIDs_35'};
% candidateLabels200CellEcdfs = {'candidateIDs_2,14,20,29-30,32,35-36'; ...
%                     'candidateIDs_2'; 'candidateIDs_14'; 'candidateIDs_20'; ...
%                     'candidateIDs_29'; 'candidateIDs_30'; 'candidateIDs_32'; ...
%                     'candidateIDs_35'; 'candidateIDs_36'};
candidateLabels200CellEcdfs = {'candidateIDs_2,7,11,13-14,20,27,29-30,32,35-36'; ...
                    'candidateIDs_2'; 'candidateIDs_7'; 'candidateIDs_11'; ...
                    'candidateIDs_13'; 'candidateIDs_14'; 'candidateIDs_20'; ...
                    'candidateIDs_27'; 'candidateIDs_29'; 'candidateIDs_30'; ...
                    'candidateIDs_32'; 'candidateIDs_35'; 'candidateIDs_36'};

% candidateLabels200CellActivationProfiles = {};
candidateLabels200CellActivationProfiles = {'candidateIDs_32'; ...
                    'candidateIDs_2,7,11,13-14,20,27,29-30,32,35-36'};

cellNameStr = 'cellName';
epasStr = 'TCepas';

figTypes = {'png', 'epsc'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Find the directory for this iteration
exampleIterDir2Cell = fullfile(networkDirectory, exampleIterName2Cell);
exampleIterDir200Cell = fullfile(networkDirectory, exampleIterName200Cell);
popIterDir2Cell = fullfile(networkDirectory, popIterName2Cell);
popIterDir200Cell = fullfile(networkDirectory, popIterName200Cell);

% Find all possible candidate labels
if isempty(candidateLabels200Cell)
    candidateLabels200Cell = find_candidate_labels(popIterDir200Cell);

    % TEMP: remove
    toRemove = contains(candidateLabels200Cell, ...
                            'candidateIDs_7,13-14,22,32,36');
    candidateLabels200Cell = candidateLabels200Cell(~toRemove);
end

% Construct the full path to the candidate cell spreadsheet
candCellSheetPath = fullfile(networkDirectory, candCellSheetName);

% Create a rank string
rankStr = ['rank', create_label_from_sequence(rankNumsToUse)];

% Create a condition label
[conditionLabel2Cell, conditionLabel200Cell] = ...
    argfun(@(x) [x, '_', rankStr, '_gIncr', num2str(gIncr)], ...
            popIterName2Cell, popIterName200Cell);

% Create a population data spreadsheet name
popDataSheetName2Cell = [popIterName2Cell, '_', rankStr, '_', ...
                            oscParamsSuffix, '.csv'];
popDataSheetName200Cell = [popIterName200Cell, '_', rankStr, '_', ...
                            oscParamsSuffix, '.csv'];

% Create a network data spreadsheet names
networkSheetNames = strcat(popIterName200Cell, '_', ...
                            candidateLabels200Cell, '_', ...
                            oscParamsSuffix, '.csv');

networkSheetNamesEcdfs = strcat(popIterName200Cell, '_', ...
                                candidateLabels200CellEcdfs, '_', ...
                                oscParamsSuffix, '.csv');

% Contruct the full path to the population data spreadsheet
popDataPath2Cell = fullfile(figure07Dir, popDataSheetName2Cell);
popDataPath200Cell = fullfile(figure08Dir, popDataSheetName200Cell);
networkDataPaths = fullfile(figure08Dir, networkSheetNames);
networkDataPathsEcdfs = fullfile(figure08Dir, networkSheetNamesEcdfs);

% Construct stats table paths
networkStatLabels = strcat(popIterName200Cell, '_', candidateLabels200Cell, ...
                            '_gIncr', num2str(gIncr));
statsGroupByEpasPaths = ...
    fullfile(figure08Dir, strcat(networkStatLabels, '_groupByEpas_stats.mat'));

% Create color maps
colorMapPharm = decide_on_colormap([], 4);
colorMapPharmCell = arrayfun(@(x) colorMapPharm(x, :), ...
                            transpose(1:4), 'UniformOutput', false);


%% Find example files and directories
if plotIpscComparison || plot2CellEssential || plot2CellM2h
    % Find example network directories
    [~, exampleDirs2Cell] = ...
        cellfun(@(x) all_subdirs('Directory', exampleIterDir2Cell, ...
                                'Keyword', x), ...
                exampleCellNames2Cell, 'UniformOutput', false);
end
if plot200CellExamples
    % Select seed number directory
    seedNumberDir = fullfile(exampleIterDir200Cell, exampleSeedNumberDirName);

    % Find example network directories
    [~, exampleDirs200Cell] = ...
        cellfun(@(x) all_subdirs('Directory', seedNumberDir, ...
                                'Keyword', x, 'Recursive', true), ...
                exampleCellNames200Cell, 'UniformOutput', false);
end

%% Plots figures for comparing dynamic clamp ipsc
if plotIpscComparison
    cellfun(@(x, y) plot_ipsc_comparison(x, exampleIterName2Cell, gIncr, y, ...
                                        figure07Dir, figTypes, ...
                                        ipscFigWidth, ipscFigHeight, ...
                                        xLimits2CellGabab, yLimitsGabab), ...
            exampleCellNames2Cell, exampleDirs2Cell);
end

%% Plots example 2-cell networks
if plot2CellEssential
    cellfun(@(a, b) ...
        cellfun(@(x, y) plot_2cell_examples(x, exampleIterName2Cell, ...
                            gIncr, a, y, figure07Dir, figTypes, ...
                            essential2CellFigWidth, essential2CellFigHeight, ...
                            xLimits2CellEssential, yLimitsEssential, ...
                            'essential', b, yTicksEssential), ...
                exampleCellNames2Cell, exampleDirs2Cell), ...
        num2cell(pharmConditions), colorMapPharmCell);
end

%% Plots m2h of example 2-cell networks
if plot2CellM2h
    cellfun(@(a, b) ...
        cellfun(@(x, y) plot_2cell_examples(x, exampleIterName2Cell, ...
                            gIncr, a, y, figure07Dir, figTypes, ...
                            m2h2CellFigWidth, m2h2CellFigHeight, ...
                            xLimits2CellEssential, yLimitsM2h, ...
                            'm2h', b, yTicksM2h), ...
                exampleCellNames2Cell, exampleDirs2Cell), ...
        num2cell(pharmConditions), colorMapPharmCell);
end

%% Plots example 200-cell networks
if plot200CellExamples
    arrayfun(@(z) ...
        cellfun(@(x, y) plot_200cell_examples(x, exampleIterName200Cell, ...
                            gIncr, z, y, figure08Dir, figTypes, ...
                        example200CellFigWidth, example200CellFigHeight), ...
                exampleCellNames200Cell, exampleDirs200Cell), ...
        pharmConditions);
end

%% Analyzes spikes for all 2-cell networks
if analyze2CellSpikes
    reanalyze_network_spikes(popIterDir2Cell, backupPrevious2Cell, ...
                                plotAnalysis2Cell);
end

%% Combines quantification over all 2-cell networks
if combine2CellPopulation
    combine_osc_params(popIterDir2Cell, candCellSheetPath, ...
                            rankNumsToUse, popDataPath2Cell, {});
end

%% Analyzes spikes for all 200-cell networks
if analyze200CellSpikes
    reanalyze_network_spikes(popIterDir200Cell, backupPrevious200Cell, ...
                                plotAnalysis200Cell);
end

%% Combines quantification over all homogeneous 200-cell networks
if combine200CellPopulation
    combine_osc_params(popIterDir200Cell, candCellSheetPath, ...
                            rankNumsToUse, popDataPath200Cell, {});
end

%% Combines activation profiles over seed numbers for each 200-cell network
if combineActivationProfiles
    combine_activation_profiles(popIterDir200Cell, figure08Dir, epasToPlot, ...
                                candidateLabels200CellActivationProfiles);
end

%% Plots oscillation measures over pharm condition 
%       across all 2-cell networks
if plot2CellViolins
    % Construct stats table path
    stats2dPath2Cell = ...
        fullfile(figure07Dir, strcat(conditionLabel2Cell, '_stats.mat'));

    % Compute statistics if not done already
    if ~isfile(stats2dPath2Cell)
        % Compute statistics for all features
        disp('Computing statistics for violin plots ...');
        statsTable = m3ha_network_compute_statistics(popDataPath2Cell, ...
                                gIncr, measuresOfInterest, measureTitles, ...
                                'mean', cellNameStr);

        % Generate labels
        conditionLabel = conditionLabel2Cell;
        pharmLabels = pharmLabelsShort;

        % Save stats table
        save(stats2dPath2Cell, 'statsTable', 'pharmLabels', ...
                            'conditionLabel', '-v7.3');
    end

    % Plot violin plots
    m3ha_plot_violin(stats2dPath2Cell, 'RowsToPlot', measuresOfInterest, ...
                    'OutFolder', figure07Dir);
end

%% Plots mean oscillation measures over pharm condition 
%       across all 200-cell networks
if plot200CellViolins
    % Construct stats table path
    stats2dPath200Cell = ...
        fullfile(figure08Dir, strcat(conditionLabel200Cell, '_stats.mat'));

    % Compute statistics if not done already
    if ~isfile(stats2dPath200Cell)
        % Compute statistics for all features
        disp('Computing statistics for violin plots ...');
        statsTable = m3ha_network_compute_statistics(popDataPath200Cell, ...
                                gIncr, measuresOfInterest, measureTitles, ...
                                'mean', cellNameStr);

        % Generate labels
        conditionLabel = conditionLabel200Cell;
        pharmLabels = pharmLabelsShort;

        % Save stats table
        save(stats2dPath200Cell, 'statsTable', 'pharmLabels', ...
                            'conditionLabel', '-v7.3');
    end

    % Plot violin plots
    m3ha_plot_violin(stats2dPath200Cell, 'RowsToPlot', measuresOfInterest, ...
                    'OutFolder', figure08Dir);
end

%% Plots oscillation measures grouped by cells
if plot200CellGroupByCellJitters
    % Construct stats table path
    statsGroupByCellPath200Cell = ...
        fullfile(figure08Dir, strcat(conditionLabel200Cell, ...
                    '_groupByCell_stats.mat'));

    % Compute statistics if not done already
    if ~isfile(statsGroupByCellPath200Cell)
        % Compute statistics for all features
        disp('Computing statistics for grouped jitter plots ...');
        statsTable = m3ha_network_compute_statistics(popDataPath200Cell, ...
                            gIncr, measuresOfInterestJitter, ...
                            measureTitlesJitter, 'grouped', cellNameStr);

        % Generate labels
        conditionLabel = conditionLabel200Cell;
        pharmLabels = pharmLabelsShort;

        % Save stats table
        save(statsGroupByCellPath200Cell, 'statsTable', 'pharmLabels', ...
                            'conditionLabel', '-v7.3');
    end

    % Plot jitter plots
    m3ha_plot_all_jitters(statsGroupByCellPath200Cell, figure08Dir, ...
                                measuresOfInterestJitter);
end

%% Combines quantification over each 200-cell networks
if combineEach200CellNetwork
    cellfun(@(dataPath, candidateLabel) ...
                combine_osc_params(popIterDir200Cell, candCellSheetPath, ...
                            rankNumsToUse, dataPath, candidateLabel), ...
            networkDataPaths, candidateLabels200Cell, 'UniformOutput', false);
end

%% Plots oscillation measures grouped by TCepas
if plot200CellGroupByEpasJitters
    % Compute statistics if not done already
    cellfun(@(a, b, c) ...
            m3ha_network_compute_and_save_statistics(a, b, gIncr, ...
                            measuresOfInterestJitter, measureTitlesJitter, ...
                            'grouped', epasStr, c, pharmLabelsShort), ...
            statsGroupByEpasPaths, networkDataPaths, networkStatLabels);

    % Plot jitter plots
    cellfun(@(a) m3ha_plot_all_jitters(a, figure08Dir, ...
                            measuresOfInterestJitter), ...
            statsGroupByEpasPaths);
end

%% Plots cumulative distribution plots
if plot200CellCumDist
   m3ha_plot_all_cumulative_distributions(networkDataPathsEcdfs, ...
                                        measuresOfInterestJitter, ...
                                        measureTitlesJitter, candCellSheetPath); 
end

%% Archive all scripts for this run
if archiveScriptsFlag
    if plot200CellExamples || analyze200CellSpikes || plotAnalysis200Cell || ...
            backupPrevious200Cell || combineActivationProfiles || ...
            combine200CellPopulation || plot200CellViolins || ...
            plot200CellGroupByCellJitters || combineEach200CellNetwork || ...
            plot200CellGroupByEpasJitters || plot200CellCumDist
        archive_dependent_scripts(mfilename, 'OutFolder', figure08Dir);
    else
        archive_dependent_scripts(mfilename, 'OutFolder', figure07Dir);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_ipsc_comparison (cellName, popIterName2Cell, gIncr, ...
                                inFolder, outFolder, ...
                                figTypes, figWidth, figHeight, ...
                                xLimits, yLimits)
% Plot an IPSC comparison plot

% Create a gIncr string
gIncrStr = ['gIncr', num2str(gIncr)];

% Create figure names
figPathBase = fullfile(outFolder, [cellName, '_', popIterName2Cell, ...
                        '_', gIncrStr, '_gabab_ipsc_comparison']);
figPathBaseOrig = [figPathBase, '_orig'];

% Create the figure
fig = set_figure_properties('AlwaysNew', true);

% Plot comparison
m3ha_network_plot_gabab('SaveNewFlag', false, 'InFolder', inFolder, ...
                        'XLimits', xLimits, 'YLimits', yLimits, ...
                        'FigTitle', 'suppress', ...
                        'AmpScaleFactor', gIncr);

% Save original figure
drawnow;
save_all_figtypes(fig, figPathBaseOrig, 'png');

% Plot a scale bar
plot_scale_bar('x', 'XBarUnits', 'ms', 'XBarLength', 200, ...
                'XPosNormalized', 0.9, 'YPosNormalized', 0.9);

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Width', figWidth, 'Height', figHeight, ...
                        'RemoveXRulers', true, 'AlignSubplots', true);

% Save the figure
drawnow;
save_all_figtypes(fig, figPathBase, figTypes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_2cell_examples (cellName, iterName, gIncr, pharm, ...
                            inFolder, outFolder, figTypes, ...
                            figWidth, figHeight, xLimits, yLimits, ...
                            plotType, colorMap, yTickLocs)
% Plot 2-cell network examples

% Create a gIncr string
gIncrStr = ['gIncr', num2str(gIncr)];
pharmStr = ['pharm', num2str(pharm)];

% Create figure names
figPathBase = fullfile(outFolder, [cellName, '_', iterName, ...
                        '_', gIncrStr, '_', pharmStr, '_2cell_', plotType]);
figPathBaseOrig = [figPathBase, '_orig'];

% Create the figure
fig = set_figure_properties('AlwaysNew', true);

% Plot example
handles = ...
    m3ha_network_plot_essential('SaveNewFlag', false, 'InFolder', inFolder, ...
                        'XLimits', xLimits, 'YLimits', yLimits, ...
                        'FigTitle', 'suppress', ...
                        'AmpScaleFactor', gIncr, 'PharmCondition', pharm, ...
                        'PlotType', plotType, 'Color', colorMap);

% Save original figure
drawnow;
save_all_figtypes(fig, figPathBaseOrig, 'png');

% Fine tune
switch plotType
case 'essential'
    % Get all subplots
    subPlots = handles.subPlots;

    % Plot a scale bar in the first subplot
    subplot(subPlots(1));
    plot_scale_bar('x', 'XBarUnits', 'ms', 'XBarLength', 200, ...
                    'XPosNormalized', 0.9, 'YPosNormalized', 0.9);
end

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Width', figWidth, 'Height', figHeight, ...
                        'YTickLocs', yTickLocs, ...
                        'RemoveXRulers', true, 'AlignSubplots', true);

% Save the figure
drawnow;
save_all_figtypes(fig, figPathBase, figTypes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_200cell_examples (cellName, iterName, gIncrDclamp, pharm, ...
                            inFolder, outFolder, figTypes, figWidth, figHeight)
% Plot 200-cell network examples

% Get the gIncr value for the network
gIncr = gIncrDclamp / 12;

% Create strings
gIncrStr = ['gIncr', num2str(gIncrDclamp)];
pharmStr = ['pharm', num2str(pharm)];

% Find the appropriate simulation number
simNumber = m3ha_network_find_sim_number(inFolder, pharm, gIncr);

% Create figure names
figPathBase = fullfile(outFolder, [cellName, '_', iterName, ...
                        '_', gIncrStr, '_', pharmStr, '_200cell_example']);
figPathBaseOrig = [figPathBase, '_orig'];

%% Full figure
% Create the figure
fig = set_figure_properties('AlwaysNew', true);

% Plot spike raster plot
m3ha_network_raster_plot(inFolder, 'OutFolder', outFolder, ...
                        'SingleTrialNum', simNumber, ...
                        'PlotSpikes', true, 'PlotTuning', false, ...
                        'PlotOnly', true);

% Save original figure
drawnow;
save_all_figtypes(fig, figPathBaseOrig, 'png');

% Update figure for CorelDraw
update_figure_for_corel(fig, 'RemoveXLabels', true, 'RemoveYLabels', true, ...
                        'RemoveTitles', true, 'RemoveXRulers', true);
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                            'Width', figWidth, 'Height', figHeight);

% Plot a scale bar only for the Dual Blockade condition
if pharm == 4
    plot_scale_bar('x', 'XBarUnits', 'sec', 'XBarLength', 2, ...
                    'XPosNormalized', 0.6, 'YPosNormalized', 0.2);
end

% Save the figure
drawnow;
save_all_figtypes(fig, figPathBase, figTypes);

% Close all figures
close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function simNumber = m3ha_network_find_sim_number (inFolder, pharm, gIncr)
%% TODO: Move this to m3ha_network_raster_plot.m

% Create strings
paramsPrefix = 'sim_params';
pharmStr = ['pCond_', num2str(pharm)];
gIncrStr = ['gIncr_', num2str(gIncr)];

% Create the keyword
keyword = [pharmStr, '_', gIncrStr];

% Find the sim params file
[~, paramPath] = all_files('Directory', inFolder, 'Prefix', paramsPrefix, ...
                            'Keyword', keyword, 'MaxNum', 1);

% Find the corresponding simNumber
paramTable = readtable(paramPath, 'ReadRowNames', true);
simNumber = paramTable{'simNumber', 'Value'};

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function reanalyze_network_spikes(popIterDir, backupPrevious, plotAnalysis)
%% Re-analyzes spikes for all networks

%% Hard-coded parameters
oscParamsSuffix = 'oscillation_params';

if backupPrevious
    % Create a backup suffix
    oscParamsBackupSuffix = ['oscillation_params_backup_', create_time_stamp];

    % Locate all oscillation parameter paths
    [~, oscParamPaths] = ...
        all_files('Directory', popIterDir, ...
                    'Suffix', oscParamsSuffix, 'Extension', 'csv', ...
                    'Recursive', true, 'ForceCellOutput', true);
                
    % Create backup paths
    oscParamBackupPaths = ...
        replace(oscParamPaths, oscParamsSuffix, oscParamsBackupSuffix);

    % Backup parameters files
    cellfun(@(x, y) movefile(x, y), oscParamPaths, oscParamBackupPaths);
end

% Find all network subdirectories
[~, netSimDirs] = all_subdirs('Directory', popIterDir, 'Level', 2);

% Analyze spikes for all network subdirectories
array_fun(@(x) m3ha_network_analyze_spikes('Infolder', x, ...
                'PlotFlag', plotAnalysis), ...
            netSimDirs, 'UniformOutput', false);
% cellfun(@(x) m3ha_network_analyze_spikes('Infolder', x, ...
%                 'PlotFlag', plotAnalysis), ...
%         netSimDirs, 'UniformOutput', false);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combinedTable = combine_osc_params (popIterDir, candCellSheetPath, ...
                                rankNumsToUse, combinedPath, candidateLabels)
% TODO: Add candidate IDs

%% Hard-coded parameters
rankNumStr = 'rankNum';
seedNumStr = 'seedNumber';
cellNameStr = 'cellName';
oscParamsSuffix = 'oscillation_params';
candIdStr = 'candidateId';
candLabelStr = 'candidateLabel';

%% Do the job
% Read the candidate cell table
candCellTable = readtable(candCellSheetPath, 'ReadRowNames', true);

% Find the cell names to use from the table
if isempty(candidateLabels)
    % Extract cell names and candidate IDs
    rowConditions = {rankNumStr, rankNumsToUse};
    [cellNamesToUse, candIdsToUse] = ...
        argfun(@(x) extract_vars(candCellTable, x, ...
                                'RowConditions', rowConditions), ...
                cellNameStr, candIdStr);

    % Create candidate labels
    candidateLabels = create_labels_from_numbers(candIdsToUse, ...
                                                'Prefix', 'candidateIDs_');
else
    % Force as a cell array
    candidateLabels = force_column_cell(candidateLabels);

    % Extract cell names to use
    cellNamesToUse = candidateLabel2cellName(candidateLabels, candCellTable);
end

% Find all seed number subdirectories
[~, seedNumDirs] = all_subdirs('Directory', popIterDir, 'Recursive', false, ...
                                'Prefix', 'seedNumber');

% Extract all tables for each seed number
oscParamTablesCell = ...
    cellfun(@(seedNumDir) retrieve_osc_param_tables(seedNumDir, ...
                        candidateLabels, cellNamesToUse, oscParamsSuffix, ...
                        seedNumStr, cellNameStr, candLabelStr), ...
            seedNumDirs, 'UniformOutput', false);

% Vertically concatenate the cell arrays
oscParamTables = apply_over_cells(@vertcat, oscParamTablesCell);

% Vertically concatenate the tables
combinedTable = apply_over_cells(@vertcat, oscParamTables);

% Join the candidate cell info to the table
if any(ismatch(cellNamesToUse, '[A-Z][0-9]{6}', 'MatchMode', 'regexp'))
    combinedTable = join(combinedTable, candCellTable, 'Keys', cellNameStr);
end

% Save the table
writetable(combinedTable, combinedPath);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cellNameToUse = candidateLabel2cellName (candidateLabel, candCellTable)

% Hard-coded parameters
candIdStr = 'candidateId';
cellNameStr = 'cellName';

if iscell(candidateLabel)
    cellNameToUse = ...
        cellfun(@(x) candidateLabel2cellName(x, candCellTable), ...
                candidateLabel, 'UniformOutput', false);
    return
end

candIdToUse = sscanf_full(candidateLabel, '%d');

if numel(candIdToUse) == 1
    cellNameToUse = extract_vars(candCellTable, cellNameStr, ...
                                'RowConditions', {candIdStr, candIdToUse});
    if iscell(cellNameToUse)
        cellNameToUse = cellNameToUse{1};
    end
else
    indNeg = find(candIdToUse < 0);

    if ~isempty(indNeg)
        rangeStart = candIdToUse(indNeg - 1);
        rangeEnd = -candIdToUse(indNeg);
        candIdsToRemove = candIdToUse(indNeg);
        candIdsToAdd = arrayfun(@(x, y) x:y, rangeStart, rangeEnd, ...
                                'UniformOutput', false);
        candIdToUse = setdiff(union(candIdToUse, ...
                                union_over_cells(candIdsToAdd)), ...
                                candIdsToRemove);
    end

    cellNameToUse = ['hetero', num2str(numel(candIdToUse))];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function candidateLabels = find_candidate_labels (popIterDir)

%% Hard-coded parameters
oscDataSuffix = 'oscillation_data';

% Find all oscillation data matfiles with this candidate label
[~, oscDataPaths] = ...
    all_files('Directory', popIterDir, 'Recursive', true, ...
                'Suffix', oscDataSuffix, 'Extension', 'mat');

% Extract all candidate label strings
candidateStrs = ...
    m3ha_extract_candidate_label(oscDataPaths, 'FromBaseName', true);

% Find unique candidate labels
candidateLabels = unique(candidateStrs);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combine_activation_profiles (popIterDir, outFolder, epasToPlot, ...
                                        candidateLabels)

cellfun(@(c) combine_activation_profiles_helper(c, popIterDir, ...
                                                outFolder, epasToPlot), ...
        candidateLabels);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combine_activation_profiles_helper (candidateLabel, popIterDir, ...
                                                outFolder, epasToPlot)

%% Hard-coded parameters
oscDataSuffix = 'oscillation_data';
seedNumStr = 'seedNumber';
seedNumLabelRegExp = [seedNumStr, '_[0-9]*'];
oscDataStr = 'oscData';
oscParamsStr = 'oscParams';

figHeight = 4.3;
figWidth = 4.3;

%% Do the job
% Find all oscillation data matfiles with this candidate label
[~, oscDataPaths] = ...
    all_files('Directory', popIterDir, 'Recursive', true, ...
                'Keyword', candidateLabel, 'Suffix', oscDataSuffix, ...
                'Extension', 'mat');

% Extract the seed number labels
seedNumLabels = extract_substrings(oscDataPaths, 'RegExp', seedNumLabelRegExp);

% Extract the base name
popIterDirName = extract_fileparts(popIterDir, 'base');

% Keep only oscillation data matfiles under a seed number directory
toKeep = ~isemptycell(seedNumLabels);
oscDataPaths = oscDataPaths(toKeep);
seedNumLabels = seedNumLabels(toKeep);

% Extract the seed numbers
seedNums = cellfun(@(x) sscanf_full(x, '%d'), seedNumLabels);

% Extract all data tables
oscDataMatFiles = cellfun(@matfile, oscDataPaths, 'UniformOutput', false);
oscDataTables = cellfun(@(m) m.(oscDataStr), ...
                        oscDataMatFiles, 'UniformOutput', false);

% Extract all condition strings from the first params table
oscParamsTable = oscDataMatFiles{1}.(oscParamsStr);
condStrs = oscParamsTable.Properties.RowNames;
nCells = oscParamsTable{:, 'nCells'};

% Create a figure path for each condition string
figPathBases = fullfile(outFolder, strcat(popIterDirName, '_', ...
                    candidateLabel, '_', condStrs, '_activation_profile'));
figTitleBases = replace(strcat(candidateLabel, '_', condStrs), '_', '\_');

% TEMP: TODO
condStrs = {1; 2; 3; 4};

% Plot mean activation profiles
cellfun(@(a, b, c, d) plot_mean_activation_profiles(a, b, c, d, ...
                        seedNums, oscDataTables, epasToPlot, ...
                        figHeight, figWidth), ...
        figPathBases, figTitleBases, condStrs, num2cell(nCells));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_mean_activation_profiles (figPathBase, figTitleBase, ...
                                        condStr, nCells, seedNums, ...
                                        oscDataTables, epasToPlot, ...
                                        figHeight, figWidth)

% Compute the corresponding TCepas value
TCepasValues = -75 + mod(seedNums, 16);

% Unique epas values
uniqueEpas = unique(TCepasValues);

% Restrict unique epas values
if ~isempty(epasToPlot)
    uniqueEpas = intersect(uniqueEpas, epasToPlot);
end

% Create epas labels
uniqueEpasLabels = create_labels_from_numbers(uniqueEpas, 'Prefix', 'epas = ');

% Count unique epas values
nEpas = numel(uniqueEpas);

% Extract the activation profiles for TC
[timeBinsSeconds, percentActivatedTC] = ...
    argfun(@(colStr) cellfun(@(x) x{condStr, colStr}{1}, ...
                        oscDataTables, 'UniformOutput', false), ...
            'timeBinsSeconds', 'percentActivatedTC');

% Group activation profiles by TCepas value
percentActivatedEachEpas = ...
    arrayfun(@(epas) percentActivatedTC(TCepasValues == epas), ...
            uniqueEpas, 'UniformOutput', false);

% Compute the combined trace for each TCepas value
[meanAct, lowerAct, upperAct] = ...
    argfun(@(method) ...
            cellfun(@(traces) compute_combined_trace(traces, method), ...
                    percentActivatedEachEpas, 'UniformOutput', false), ...
            'mean', 'lower95', 'upper95');

% Force as matrices
[meanAct, lowerAct, upperAct] = ...
    argfun(@force_matrix, meanAct, lowerAct, upperAct);

% Decide on colors
% TODO: 'ForceCellOutput' for decide_on_colormap.m
% colors = decide_on_colormap([], nEpas);
% colorsCell = arrayfun(@(i) colors(i, :), transpose(1:nEpas), ...
%                         'UniformOutput', false);

% Create a figure
fig = set_figure_properties('AlwaysNew', true);

% Hold on
hold on;

% Plot the mean activation profiles
handles = plot_tuning_curve(timeBinsSeconds{1}, meanAct, ...
                            'LowerCI', lowerAct, 'UpperCI', upperAct, ...
                            'ColumnLabels', uniqueEpasLabels, ...                            
                            'LineWidth', 1);

ylim([0, nCells]);
xlabel('Time (s)');
ylabel('Percent Activated (%)');
title(['Activation profile for ', figTitleBase]);
legend(handles.curves, 'location', 'southeast');

% Save the figure
save_all_figtypes(fig, [figPathBase, '_orig.png'], 'png');

% Update for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Height', figHeight, 'Width', figWidth);

% Save the figure
save_all_figtypes(fig, [figPathBase, '.png'], {'png', 'epsc'});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function oscParamTables = ...
                retrieve_osc_param_tables (seedNumDir, candidateLabelsToUse, ...
                                        cellNamesToUse, oscParamsSuffix, ...
                                        seedNumStr, cellNameStr, candLabelStr)

% Hard-coded parameters
epasStr = 'TCepas';

% Force as a column cell array
cellNamesToUse = force_column_cell(cellNamesToUse);

% Extract the seed number
seedNumber = sscanf_full(extract_fileparts(seedNumDir, 'base'), '%d');

% Compute the TC epas
TCepas = -75 + mod(seedNumber, 16);

% Locate corresponding oscillation parameter paths
[~, oscParamPaths] = ...
    find_matching_files(candidateLabelsToUse, 'Directory', seedNumDir, ...
                        'Suffix', oscParamsSuffix, 'Extension', 'csv', ...
                        'Recursive', true, 'ForceCellOutput', true);

% Read the oscillation parameter tables
oscParamTables = cellfun(@readtable, oscParamPaths, 'UniformOutput', false);

% Add the TC epas to the tables
oscParamTables = ...
    cellfun(@(x, y) addvars_custom(x, TCepas, ...
                            'NewVariableNames', epasStr, 'Before', 1), ...
            oscParamTables, 'UniformOutput', false);

% Add the seed number to the tables
oscParamTables = ...
    cellfun(@(x, y) addvars_custom(x, seedNumber, ...
                            'NewVariableNames', seedNumStr, 'Before', 1), ...
            oscParamTables, 'UniformOutput', false);

% Add the candidate label to the tables
oscParamTables = ...
    cellfun(@(x, y) addvars_custom(x, {y}, 'NewVariableNames', candLabelStr, ...
                                    'Before', 1), ...
            oscParamTables, candidateLabelsToUse, 'UniformOutput', false);

% Add the cell name to the tables
oscParamTables = ...
    cellfun(@(x, y) addvars_custom(x, {y}, 'NewVariableNames', cellNameStr, ...
                                    'Before', 1), ...
            oscParamTables, cellNamesToUse, 'UniformOutput', false);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m3ha_network_compute_and_save_statistics(statsPath, dataPath, gIncr, ...
                            measuresOfInterestJitter, measureTitlesJitter, ...
                            method, groupNameStr, conditionLabel, pharmLabels)

if ~isfile(statsPath)
    % Compute statistics for all features
    disp('Computing statistics for grouped jitter plots ...');
    statsTable = m3ha_network_compute_statistics(dataPath, ...
                        gIncr, measuresOfInterestJitter, ...
                        measureTitlesJitter, method, groupNameStr);

    % Save stats table
    save(statsPath, 'statsTable', 'pharmLabels', ...
                        'conditionLabel', '-v7.3');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function statsTable = m3ha_network_compute_statistics (popDataPath, gIncr, ...
                                measureStr, measureTitle, method, groupNameStr)
%% Computes all statistics for the 2-cell network

%% Hard-coded parameters
if isempty(groupNameStr)
    groupNameStr = 'cellName';
end
seedNumStr = 'seedNumber';
gIncrStr = 'gIncr';
pharmStr = 'pCond';
dclamp2NetworkAmpRatio = 12;

%% Do the job
% Read the data table
popDataTable = readtable(popDataPath);

% Restrict to the rows with given gIncr
rowsToUse = round(popDataTable.(gIncrStr) * dclamp2NetworkAmpRatio) == gIncr;

% Change the measure strings the original non-averaged strings
switch method
    case 'mean'
        measureStrNoMean = ...
            replace(measureStr, {'mean', 'oscillationProbability'}, ...
                                            {'', 'hasOscillation'});
        measureStrOrig = lower_first_char(measureStrNoMean);
    case 'grouped'
        measureStrOrig = measureStr;
    otherwise
        error('Not implemented yet!');
end

% Locate the columns of interest
colsOfInterest = [{groupNameStr}; {seedNumStr}; {pharmStr}; measureStrOrig];

% Extract the table of interest
popTableOfInterest = popDataTable(rowsToUse, colsOfInterest);

% Compute statistics for each measure of interest
[allValues, pharmCondition, uniqueGroupValues] = ...
    cellfun(@(x) m3ha_network_stats_helper(popTableOfInterest, seedNumStr, ...
                                        pharmStr, groupNameStr, x, method), ...
                    measureStrOrig, 'UniformOutput', false);

% Create the statistics table
statsTable = table(measureTitle, measureStr, pharmCondition, ...
                    uniqueGroupValues, allValues, 'RowNames', measureStr);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [allValuesEachPharm, pharmCondition, uniqueGroupValues] = ...
                m3ha_network_stats_helper (popDataTable, seedNumStr, ...
                                    pharmStr, groupNameStr, measureStr, method)
%% Computes the statistics for one measure

%% Do the job
% Extract from table
pharmAll = popDataTable.(pharmStr);
groupValueAll = popDataTable.(groupNameStr);

% Get all possible pharmacological conditions
pharmCondition = force_column_cell(num2cell(unique(pharmAll, 'sorted')));

% Get all group values
if iscell(groupValueAll)
    uniqueGroupValues = ...
        force_column_cell(unique_custom(groupValueAll, 'IgnoreNaN', true));
else
    uniqueGroupValues = ...
        force_column_vector(unique_custom(groupValueAll, 'IgnoreNaN', true));
end

% Find corresponding row numbers
if iscell(uniqueGroupValues)
    rowsEachGroupEachPharm = ...
        cellfun(@(p) cellfun(@(c) pharmAll == p & strcmp(groupValueAll, c), ...
                            uniqueGroupValues, 'UniformOutput', false), ...
                pharmCondition, 'UniformOutput', false);
else
    rowsEachGroupEachPharm = ...
        cellfun(@(p) arrayfun(@(c) pharmAll == p & ismatch(groupValueAll, c), ...
                                uniqueGroupValues, 'UniformOutput', false), ...
                pharmCondition, 'UniformOutput', false);
end

% Get mean values across iterations for all cells
%    for each pharm condition, for this measure
switch method
    case 'mean'
        allValuesEachPharm = ...
            cellfun(@(a) ...
                    cellfun(@(b) nanmean(popDataTable{b, measureStr}), a), ... 
                rowsEachGroupEachPharm, 'UniformOutput', false);
    case 'grouped'
        allValuesEachPharm = ...
            cellfun(@(a) ...
                    cellfun(@(b) popDataTable{b, measureStr}, ...
                            a, 'UniformOutput', false), ... 
                rowsEachGroupEachPharm, 'UniformOutput', false);
    otherwise
        error('Not implemented yet!');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m3ha_plot_all_jitters (statsPath, outFolder, rowsToPlot)
% TODO: Pull out as its own function

%% Hard-coded parameters
figWidth = 6;
figHeight = 3;
figTypes = {'png', 'epsc'};
otherArguments = struct;

%% Preparation
% Set default output directory
if isempty(outFolder)
    outFolder = extract_fileparts(statsPath, 'directory');
end

% Load stats table
disp('Loading statistics for grouped jitter plots ...');
if isfile(statsPath)
    load(statsPath, 'statsTable', 'pharmLabels', 'conditionLabel');
else
    fprintf('%s does not exist!\n', statsPath);
    return;
end

% Restrict to measures to plot
if ~(ischar(rowsToPlot) && strcmp(rowsToPlot, 'all'))
    statsTable = statsTable(rowsToPlot, :);
end

% Extract variables
allMeasureTitles = statsTable.measureTitle;
allMeasureStrs = statsTable.measureStr;
allValues = statsTable.allValues;
uniqueGroupValues = statsTable.uniqueGroupValues;

% Create figure bases
allFigBases = combine_strings({allMeasureStrs, conditionLabel});

% Create full path bases
allFigPathBases = fullfile(outFolder, allFigBases);

% Create figure title
if contains(conditionLabel, 'candidate')
    figTitle = extractFrom(conditionLabel, 'candidate');
elseif contains(conditionLabel, 'rank')
    figTitle = extractFrom(conditionLabel, 'rank');
else
    figTitle = conditionLabel
end
figTitle = replace(figTitle, '_', '\_');

%% Do the job
% Plot all grouped jitter plots
disp('Plotting grouped jitter plots ...');
handles = ...
    cellfun(@(a, b, c, d) m3ha_plot_grouped_jitter(...
                            a, b, c, pharmLabels, d, figTitle, ...
                            figHeight, figWidth, ...
                            figTypes, otherArguments), ...
            allValues, uniqueGroupValues, allMeasureTitles, allFigPathBases);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = m3ha_plot_grouped_jitter (allValues, uniqueGroupValues, ...
                                                measureTitle, pharmLabels, ...
                                                figPathBase, figTitle, ...
                                                figHeight, figWidth, ...
                                                figTypes, otherArguments)

% Hard-coded parameters
xTickAngle = 320;

% Create figure
fig = set_figure_properties('AlwaysNew', true);

% Convert to character array or a cell array of character arrays
uniqueGroupLabels = convert_to_char(uniqueGroupValues);

% Plot groups as a grouped jitter plot
jitters = plot_grouped_jitter(allValues, 'XTickLabels', pharmLabels, ...
                        'XTickAngle', xTickAngle, 'YLabel', measureTitle, ...
                        'GroupingLabels', uniqueGroupLabels, otherArguments);

% Create title
title(figTitle);

% Save the figure
save_all_figtypes(fig, [figPathBase, '_orig'], 'png');

% Set y axis limits based on measureTitle
switch measureTitle
    case 'LTS probability'
        ylim([0, 1]);
    case 'LTS onset time (ms)'
        ylim([0, 2000]);
    case 'Spikes Per LTS'
        ylim([0, 6.5]);
    case 'LTS maximum slope (V/s)'
        ylim([0, 5]);
    case 'LTS amplitude (mV)'
        ylim([-75, -45]);
        yticks(-75:10:-45);
    otherwise
        % Do nothing
end

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Height', figHeight, 'Width', figWidth, ...
                        'RemoveTitle', true, 'RemoveLegend', true, ...
                        'ScatterMarkerSize', 3);

% Save the figure
save_all_figtypes(fig, figPathBase, figTypes);

% Save in handles
handles.fig = fig;
handles.jitters = jitters;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m3ha_plot_all_cumulative_distributions (networkDataPaths, ...
                        measuresOfInterest, measureLabels, candCellSheetPath)

% Hard-coded parameters
groupVarName = 'pCond';
figWidth = 4.3 * 2;
figHeight = 4.3 * 2;
figTypes = {'png', 'epsc'};
ecdfsSuffix = 'ecdfs';
yLabel = '% of simulations';
% linkAxesOption = 'xy';
linkAxesOption = 'y';

% Read all tables
dataTables = cellfun(@readtable, networkDataPaths, 'UniformOutput', false);

% Extract all candidate label strings
candidateLabels = ...
    m3ha_extract_candidate_label(networkDataPaths, 'FromBaseName', true);

% Read the candidate cell table
candCellTable = readtable(candCellSheetPath, 'ReadRowNames', true);

% Find corresponding cell names
cellNamesToUse = candidateLabel2cellName(candidateLabels, candCellTable);

% Extract common directory
commonDir = extract_common_directory(networkDataPaths);

% Extract iteration string
beforeCandStr = extractBefore(networkDataPaths, '_candidateIDs');
beforeCandStr = extract_fileparts(beforeCandStr, 'base');
iterStr = extract_common_prefix(beforeCandStr);

% Create figure path bases
figPathBases = fullfile(commonDir, strcat(measuresOfInterest, '_', ...
                                            iterStr, '_', ecdfsSuffix));

% Plot for each measure of interest
cellfun(@(a, b, c) plot_cumulative_distributions (dataTables, a, ...
                                            groupVarName, cellNamesToUse, ...
                                            b, yLabel, ...
                                            linkAxesOption, ...
                                            c, figTypes, ...
                                            figWidth, figHeight), ...
        measuresOfInterest, measureLabels, figPathBases, ...
        'UniformOutput', false);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_cumulative_distributions (dataTables, measureOfInterest, ...
                                            groupVarName, tableNames, ...
                                            measureLabel, yLabel, ...
                                            linkAxesOption, ...
                                            figPathBase, figTypes, ...
                                            figWidth, figHeight)
% TODO: Pull out as its own function

% Force as a column cell array
dataTables = force_column_cell(dataTables);

% Extract the grouping vectors
groupingVecs = cellfun(@(T) T.(groupVarName), dataTables, ...
                        'UniformOutput', false);

% Find the unique group values and sort them
uniqueGroupValues = unique_custom(union_over_cells(groupingVecs));

% Count the number of groups
nGroups = numel(uniqueGroupValues);

% Count the number of tables
nTables = numel(dataTables);

% Decide on a color map
colorMap = decide_on_colormap([], nTables, 'ForceCellOutput', true);

% Extract the values from each table for each group
dataEachGroup = ...
    arrayfun(@(g) extract_vars(dataTables, measureOfInterest, ...
                                    'RowConditions', {groupVarName, g}), ...
            uniqueGroupValues, 'UniformOutput', false);

% Compute ecdfs for each table and each group
[ecdfValuesEachGroup, xValuesEachGroup] = ...
    cellfun(@(x) compute_ecdfs(x), ...
            dataEachGroup, 'UniformOutput', false);

% Set default y axis label
if isempty(yLabel)
    yLabel = 'Cumulative %';
end

% Create figure
[fig, ax] = create_subplots(nGroups, 'AlwaysNew', true);

% Plot all cdfs for each group
ecdfs = cell(nGroups, 1);
for iGroup = 1:nGroups
    % Make current axes
    subplot(ax(iGroup));

    % Extract all cdfs and x values for this gruop
    ecdfValuesEachTable = ecdfValuesEachGroup{iGroup};
    xValuesEachTable = xValuesEachGroup{iGroup};
    groupValueThis = uniqueGroupValues(iGroup);
    groupValueThisStr = convert_to_char(groupValueThis);

    % Hold on
    wasHold = hold_on;

    % Plot stairs
    ecdfsThis = cellfun(@(a, b, c) stairs(a, 100 * b, 'Color', c, ...
                                            'LineWidth', 0.5), ...
                        xValuesEachTable, ecdfValuesEachTable, colorMap);

    % Hold off
    hold_off(wasHold);

    % X axis limits
    % TODO

    % Y axis limits
    ylim([0, 100]);

    % X axis Label
    xlabel(measureLabel);

    % Y axis Label
    ylabel(yLabel);

    % Title
    title(sprintf('%s = %s', groupVarName, groupValueThisStr));

    % Modify line widths
    set(ecdfsThis(contains(tableNames, 'hetero')), 'LineWidth', 2);

    % Legend
    legend(tableNames, 'location', 'northeast');

    ecdfs{iGroup} = ecdfsThis;
end

% Link the x and y axes
linkaxes(ax, linkAxesOption);

% Save the figure
save_all_figtypes(fig, [figPathBase, '_orig'], 'png');

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Height', figHeight, 'Width', figWidth);

% Save the figure
save_all_figtypes(fig, figPathBase, figTypes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ecdfValues, xValues] = compute_ecdfs(data)
%% Computes ecdfs for a group of vectors and aligns the x limits
% TODO: Pull out as its own function

%% Preparation
% Force as volumn vectors
data = force_column_vector(data, 'IgnoreNonVectors', true);

%% Do the job
if iscell(data)
    % Compute ecdf values for each vector
    [ecdfValues, xValues] = ...
        cellfun(@(y) ecdf(y), data, 'UniformOutput', false);

    % Compute the minimum and maximum x values
    minXValue = apply_iteratively(@min, xValues);
    maxXValue = apply_iteratively(@max, xValues);

    % Make sure the vectors include the boundary values
    [ecdfValues, xValues] = ...
        cellfun(@(x, y) add_boundary_values(x, y, minXValue, maxXValue), ...
                ecdfValues, xValues, 'UniformOutput', false);
else
    [ecdfValues, xValues] = ecdf(data);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ecdfValues, xValues] = add_boundary_values (ecdfValues, xValues, ...
                                                        minXValue, maxXValue)
% Add left boundary if necessary
if ecdfValues(1) ~= 0 || xValues(1) ~= minXValue
    xValues = [minXValue; xValues];
    ecdfValues = [0; ecdfValues];
end

% Add right boundary if necessary
if ecdfValues(end) ~= 1 || xValues(end) ~= maxXValue
    xValues = [xValues; maxXValue];
    ecdfValues = [ecdfValues; 1];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
