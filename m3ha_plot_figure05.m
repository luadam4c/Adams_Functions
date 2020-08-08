% m3ha_plot_figure05.m
%% Plots Figure 05, Figure 06 & Figure 07 for the GAT Blocker paper
%
% Requires:
%       cd/archive_dependent_scripts.m
%       cd/check_dir.m
%       cd/create_labels_from_numbers.m
%       cd/create_plot_movie.m
%       cd/create_time_vectors.m
%       cd/decide_on_colormap.m
%       cd/convert_to_char.m
%       cd/extract_fileparts.m
%       cd/find_matching_files.m
%       cd/ismember_custom.m
%       cd/m3ha_load_sweep_info.m
%       cd/m3ha_compute_gabab_ipsc.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/plot_scale_bar.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/update_neuron_scripts.m
%       cd/update_figure_for_corel.m

% File History:
% 2019-12-29 Modified from m3ha_plot_figure03.m
% 2020-01-29 Separated outputs to figure05Dir and figure06Dir
% 2020-04-09 Changed the y axis limits of m2h discrepancy plot to [1e-1, 1e3]
% 2020-04-13 Added plotVoltageVsOpdFig5
% 2020-05-04 Added createPlotMovie
% 2020-08-04 Splitted Figure06 to Figure06 and Figure07
% 2020-08-06 Updated pharmLabelsLong
% 2020-08-06 Updated legend placement in movies

%% Hard-coded parameters
% Flags
updateScripts = false; %true;
simulateIpscr = false; %true;
plotAllVoltages = false; %true;
plotAllTotalCurrents = false; %true;
plotAllComponentCurrents = false; %true;
plotDend2ITproperties = false; %true;

simulateTauhModes = false; %true;
plotSomaVoltage = false; %true;

computeIpscVariation = false; %true;
simulateIpscVariation = false; %true;
plotEssential = false; %true;

plotM2hFig5 = false; %true;
plotM2hTauh = false; %true;
plotM2hGabab = false; %true;

plotVoltageVsOpdFig5 = false; %true;
plotVoltageVsOpdTauh = false; %true;
plotVoltageVsOpdGabab = false; %true;

createPlotMovieFig5 = false; %true;
createPlotMovieTauh = false;
createPlotMovieGabab = false; %true;

createPhasePlotOnlyMovieFig5 = false; %true;
createPhasePlotOnlyMovieTauh = false;
createPhasePlotOnlyMovieGabab = false; %true;

simulateNoITSoma = false; %true;

archiveScriptsFlag = true;

% Directories
parentDirectory = fullfile('/media', 'adamX', 'm3ha');
figure02Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure02');
figure05Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure05');
figure06Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure06');
figure07Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure07');
gababIpscDir = figure07Dir;
fitDirectory = fullfile(parentDirectory, 'optimizer4gabab');

% Files
sweepInfoFile = 'dclampdatalog_take4.csv';
datalogPath = fullfile(figure02Dir, sweepInfoFile);
paramFileSuffix = 'params';

% Analysis settings
% exampleCellNames = {'D101310'; 'C101210'};
% exampleCellNames = {'D101310'; 'M101210'};
% exampleCellNames = {'D101310'; 'G101310'};
exampleCellNames = {'D101310'};

% Must be consistent with m3ha_compute_gabab_ipsc.m
% gababIpscSheetBases = {'gababipsc_gat3_vary_amp', ...
%                         'gababipsc_dual_vary_amp', ...
%                         'gababipsc_gat3_vary_tau', ...
%                         'gababipsc_dual_vary_tau', ...
%                         'gababipsc_vary_dual_to_gat3_to_gat1', ...
%                         'gababipsc_original', ...
%                         'gababipsc_gat3_vary_amp2'};
gababIpscSheetBases = {'gababipsc_dual_vary_tau'};

% Simulation settings
dataModeIpscr = 2;                  % data mode for IPSC response
                                    %   0 - all data
                                    %   1 - all of g incr = 100%, 200%, 400% 
                                    %   2 - same g incr but exclude 
                                    %       cell-pharm-g_incr sets 
                                    %       containing problematic sweeps
rowmodeIpscr = 1;                   % row mode for IPSC response
                                    %   1 - each row is a pharm condition
                                    %   2 - each row is a pharm, g incr pair
attemptNumberIpscr = 7;             % attempt number for IPSC response
                                    %   1 - Use 4 traces @ 200% gIncr 
                                    %           for this data mode
                                    %   2 - Use all traces @ 200% gIncr 
                                    %           for this data mode
                                    %   3 - Use all traces for this data mode
                                    %   4 - Use 1 trace for each pharm x gIncr 
                                    %           pair for this data mode
                                    %   5 - Use 4 traces @ 400% gIncr 
                                    %       for this data mode
                                    %   6 - Same as 4 but prioritize least vHold
                                    %   7 - Same as 1 but prioritize least vHold
                                    %   8 - Same as 5 but prioritize least vHold
                                    %   9 - Use 4 traces @ 100% gIncr 
                                    %       for this data mode
                                    %   10 - Same as 7 but prioritize least actVHold
% tauhModesAll = 4:5;
% tauhModesAll = 1:5;
% tauhModesAll = 6:7;
tauhModesAll = 1:7;

gababIpscSetNumber = 2; %1;

newParamsNoITSoma = {'pcabarITSoma', 0};

% The following must be consistent with singleneuron4compgabab.hoc
timeToStabilize = 3000;         % padded time (ms) to make sure initial value 
                                %   of simulations are stabilized

% Plot settings
colorMapPharm = [];                 % use m3ha default
colorMapVary = @jet;                % rainbow colors

overlappedFig5Width = 10.5; %4.7; %5.7;  % For Figure 5
% overlappedFig5Width = 8;       % For Figure 5 Supplement
overlappedTauhWidth = 8;
overlappedGababWidth = 4;
overlappedFigHeightPerRow = 1.5;
overlappedXLimits = timeToStabilize + [800, 2800]; %[800, 2000];
allVoltagesYLimits = {[-95, -20], [-95, -20], [-95, -20], [-95, -20], ...
                        [-4, 2], [-0.3, 0.3], [0, 10], [-4, 2]};
allTotalCurrentsYLimits = {[-4, 2], [-4, 2], [-0.3, 0.3], [-15, 5], ...
                            [-0.3, 0.3], [-5, 15], [-0.3, 0.3], [-0.3, 0.3]};
allComponentCurrentsYLimits = {[-15, 5], [-10, 5], [-10, 5], [-10, 5], ...
                            [-5, 15], [-5, 10], [-5, 10], [-5, 10]};
dend2ITpropertiesYLimits = {[-10, 5], [0, 1], [0, 1], [0, 1], ...
                            [0, 1], [1e-7, 1e0], [1e-7, 1e0], ...
                            [1e-8, 1e0], [1e-8, 1e0], [1e-1, 1e2], ...
                            [0, 0.03], [-0.001, 0.001]};
somaVoltageYLimits = {[-95, -25], [1e-8, 1e0]};
essentialYLimits = {[-110, -40], [0, 10], [-0.5, 0.1], ...
                            [-20, 5], [1e-8, 1e0]};

allVoltagesYTickLocs = {-80:20:-40, -80:20:-40, -80:20:-40, -80:20:-40, ...
                        -3:2:1, -0.2:0.2:0.2, 0:5:10, -3:2:1};
allTotalCurrentsYTickLocs = {-3:2:1, -3:2:1, -0.2:0.2:0.2, -10:5:5, ...
                            -0.2:0.2:0.2, -5:5:10, -0.2:0.2:0.2, -0.2:0.2:0.2};
allComponentCurrentsYTickLocs = {-10:5:5, -5:5:5, -5:5:5, -5:5:5, ...
                            -5:5:10, -5:5:5, -5:5:5, -5:5:5};
dend2ITpropertiesYTickLocs = {-5:5:5, 0:0.5:1, 0:0.5:1, 0:0.5:1, ...
                            0:0.5:1, [1e-6, 1e-1], [1e-6, 1e-1], ...
                            [1e-6, 1e-2], [1e-6, 1e-2], [1e0, 1e1], ...
                            [], []};
somaVoltageYTickLocs = {-90:20:-50, [1e-6, 1e-2]};
essentialYTickLocs = {-90:20:-50, 0:5:10, -0.4:0.2:0, ...
                            -15:5:0, [1e-6, 1e-2]};
m2hFig5Width = 10.5; %4.7; %5.7;
m2hTauhWidth = 8;
m2hGababWidth = 4;
m2hFigHeight = 3;
m2hXLimits = timeToStabilize + [800, 2800]; %[800, 2000];
m2hYLimits = [1e-7, 1e0];
m2hYTickLocs = [1e-5, 1e-3, 1e-1];

voltageVsOpdTimeLimits1 = timeToStabilize + [800, 2800];
voltageVsOpdTimeLimits2 = timeToStabilize + [1000, 2000];
voltageVsOpdSiMs = 1;
voltageVsOpdFig5FigWidth = 5.5 * 2;
voltageVsOpdFig5FigHeight = 5 * 3;
voltageVsOpdFig5XLimits = [1e-7, 1e0];
voltageVsOpdFig5YLimits = [-95, -45];
voltageVsOpdFig5YTickLocs = [];
voltageVsOpdFig5ToAnnotate = true;

voltageVsOpdTauhFigWidth = 8;
voltageVsOpdTauhFigHeight = 11.5;
voltageVsOpdTauhXLimits = [1e-7, 1e0];
voltageVsOpdTauhYLimits = [-95, -25];
voltageVsOpdTauhYTickLocs = -95:5:-25;
voltageVsOpdTauhToAnnotate = true;

voltageVsOpdGababFigWidth = 3.5 * 2; %4.7 * 2;
voltageVsOpdGababFigHeight = 2.2 * 3;
voltageVsOpdGababXLimits = [1e-7, 1e0];
voltageVsOpdGababYLimits = [-95, -45];
voltageVsOpdGababYTickLocs = [];
voltageVsOpdGababToAnnotate = false;

pharmLabelsLong = {'{\it sim}Control'; '{\it sim}GAT1-Block'; ...
                    '{\it sim}GAT3-Block'; '{\it sim}Dual-Block'};
traceLabelsFig5 = pharmLabelsLong;
traceLabelsTauh = pharmLabelsLong;
traceLabelsGabab = {};
legendLocationFig5 = 'best';
legendLocationTauh = 'best';
legendLocationGabab = 'suppress';

figTypes = {'png', 'epsc'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Make sure NEURON scripts are up to date in figure05Dir
if updateScripts
    update_neuron_scripts(fitDirectory, figure05Dir);
end

%% Load sweep info
% Read from datalogPath
swpInfo = m3ha_load_sweep_info('Directory', figure02Dir);

%% Find NEURON parameter tables
if simulateIpscr || simulateTauhModes || simulateIpscVariation || ...
        plotEssential || plotAllVoltages || plotAllTotalCurrents || ...
        plotAllComponentCurrents || plotDend2ITproperties || ...
        plotSomaVoltage || plotM2hFig5 || plotM2hTauh || plotM2hGabab || ...
        plotVoltageVsOpdFig5 || plotVoltageVsOpdTauh || ...
        plotVoltageVsOpdGabab || ...
        createPlotMovieFig5 || createPlotMovieTauh || ...
        createPlotMovieGabab || createPhasePlotOnlyMovieFig5 || ...
        createPhasePlotOnlyMovieTauh || createPhasePlotOnlyMovieGabab || ...
        simulateNoITSoma

    % Find NEURON parameter tables
    [~, exampleParamPaths] = ...
        find_matching_files(exampleCellNames, 'Directory', figure05Dir, ...
                            'Suffix', paramFileSuffix, 'Extension', 'csv', ...
                            'Recursive', false, 'ForceCellOutput', true);

    % Extract file bases
    exampleParamFileBases = extract_fileparts(exampleParamPaths, 'base');

    % Extract example labels
    exampleLabels = extractBefore(exampleParamFileBases, ...
                                    ['_', paramFileSuffix]);


    % Create tauhMode suffixes
    tauhModeSuffixes = create_labels_from_numbers(tauhModesAll, ...
                                                    'Prefix', 'tauhmode');

    % Update labels for each type of simulation
    exampleLabelsIpscr = strcat(exampleLabels, '_ipscr');
    exampleLabelsModeAll = cellfun(@(x) strcat(exampleLabels, '_', x), ...
                                    tauhModeSuffixes, 'UniformOutput', false);
    exampleLabelsVaryAll = cellfun(@(x) strcat(exampleLabels, '_', x), ...
                                gababIpscSheetBases, 'UniformOutput', false);
    exampleLabelsNoITSoma = strcat(exampleLabels, '_ipscr_no_ITsoma');

    % Create output folder names
    outFoldersIpscr = fullfile(figure05Dir, exampleLabelsIpscr);
    outFoldersModeAll = cellfun(@(x) fullfile(figure06Dir, x), ...
                                exampleLabelsModeAll, 'UniformOutput', false);
    outFoldersVaryAll = cellfun(@(x) fullfile(figure07Dir, x), ...
                                exampleLabelsVaryAll, 'UniformOutput', false);
    outFoldersNoITSoma = fullfile(figure05Dir, exampleLabelsNoITSoma);
end

%% Simulate regular IPSC responses
if simulateIpscr
    cd(figure05Dir);
    check_dir(outFoldersIpscr);
    cellfun(@(x, y, z) simulate_ipscr(x, y, z, 0, dataModeIpscr, ...
                                    rowmodeIpscr, attemptNumberIpscr), ...
            exampleLabelsIpscr, exampleParamPaths, outFoldersIpscr);
end

%% Simulate tauhMode == 1, 2 and 3
if simulateTauhModes
    cd(figure05Dir);
    check_dir([outFoldersModeAll{:}]);
    for iMode = 1:numel(tauhModesAll)
        cellfun(@(x, y, z) simulate_ipscr(x, y, z, tauhModesAll(iMode), ...
                    dataModeIpscr, rowmodeIpscr, attemptNumberIpscr), ...
                exampleLabelsModeAll{iMode}, exampleParamPaths, ...
                outFoldersModeAll{iMode});
    end
end

%% Compute all GABAB IPSC parameters and plot them
if computeIpscVariation
    m3ha_compute_gabab_ipsc(gababIpscSetNumber, 'OutFolder', gababIpscDir);
end

%% Simulate IPSC variation
if simulateIpscVariation
    cd(figure05Dir);
    for iSheet = 1:numel(gababIpscSheetBases)
        % Construct full path to GABA-B IPSC parameters spreadsheet
        gababIpscSheetPath = fullfile(gababIpscDir, ...
                                [gababIpscSheetBases{iSheet}, '.csv']);

        % Read GABA-B IPSC parameters table
        gababTable = readtable(gababIpscSheetPath);

        % Convert to a scalar structure
        gababStruct = table2struct(gababTable, 'ToScalar', true);

        % Simulate for each cell
        cellfun(@(x, y, z) simulate_variation(x, y, z, gababStruct), ...
                exampleLabelsVaryAll{iSheet}, exampleParamPaths, ...
                outFoldersVaryAll{iSheet});
    end
end

%% Plot all voltages
if plotAllVoltages
    cellfun(@(x, y) plot_overlapped(x, y, 'allVoltages', ...
                    figure05Dir, figTypes, ...
                    overlappedFig5Width, ...
                    numel(allVoltagesYLimits) * overlappedFigHeightPerRow, ...
                    overlappedXLimits, allVoltagesYLimits, ...
                    colorMapPharm, allVoltagesYTickLocs), ...
            exampleLabelsIpscr, outFoldersIpscr);
end

%% Plot all currents
if plotAllTotalCurrents
    cellfun(@(x, y) plot_overlapped(x, y, 'allTotalCurrents', ...
                    figure05Dir, figTypes, ...
                    overlappedFig5Width, ...
                    numel(allTotalCurrentsYLimits) * overlappedFigHeightPerRow, ...
                    overlappedXLimits, allTotalCurrentsYLimits, ...
                    colorMapPharm, allTotalCurrentsYTickLocs), ...
            exampleLabelsIpscr, outFoldersIpscr);
end

%% Plot component currents
if plotAllComponentCurrents
    cellfun(@(x, y) plot_overlapped(x, y, 'allComponentCurrents', ...
                    figure05Dir, figTypes, ...
                    overlappedFig5Width, ...
                    numel(allComponentCurrentsYLimits) * overlappedFigHeightPerRow, ...
                    overlappedXLimits, allComponentCurrentsYLimits, ...
                    colorMapPharm, allComponentCurrentsYTickLocs), ...
            exampleLabelsIpscr, outFoldersIpscr);
end

%% Plot all T channel properties
if plotDend2ITproperties
    cellfun(@(x, y) plot_overlapped(x, y, 'dend2ITproperties', ...
                    figure05Dir, figTypes, ...
                    overlappedFig5Width, ...
                    numel(dend2ITpropertiesYLimits) * overlappedFigHeightPerRow, ...
                    overlappedXLimits, dend2ITpropertiesYLimits, ...
                    colorMapPharm, dend2ITpropertiesYTickLocs), ...
            exampleLabelsIpscr, outFoldersIpscr);
end

%% Plot only soma voltages
if plotSomaVoltage
    for iMode = 1:numel(tauhModesAll)
        cellfun(@(x, y) plot_overlapped(x, y, 'somaVoltage', ...
                        figure06Dir, figTypes, ...
                        overlappedTauhWidth, ...
                        numel(somaVoltageYLimits) * overlappedFigHeightPerRow, ...
                        overlappedXLimits, somaVoltageYLimits, ...
                        colorMapPharm, somaVoltageYTickLocs), ...
                exampleLabelsModeAll{iMode}, outFoldersModeAll{iMode});
    end
end

%% Plot essential plots
if plotEssential
    for iSheet = 1:numel(gababIpscSheetBases)
        essentialYLimitsTemp = essentialYLimits;
        essentialYTicksTemp = essentialYTickLocs;
        if contains(gababIpscSheetBases{iSheet}, ...
                    {'dual_to_gat3_to_gat1', 'vary_amp2'})
            essentialYLimitsTemp{2} = [0, 35];
            essentialYTicksTemp{2} = 0:15:30;
        end

        cellfun(@(x, y) plot_overlapped(x, y, 'essential', ...
                        figure07Dir, figTypes, ...
                        overlappedGababWidth, ...
                        numel(essentialYLimitsTemp) * overlappedFigHeightPerRow, ...
                        overlappedXLimits, essentialYLimitsTemp, ...
                        colorMapVary, essentialYTicksTemp), ...
                exampleLabelsVaryAll{iSheet}, outFoldersVaryAll{iSheet});
    end
end

%% Plot m2h in dendrite 2 against its steady state
if plotM2hFig5
    cellfun(@(x, y) plot_m2h(x, y, figure05Dir, figTypes, ...
                                m2hFig5Width, m2hFigHeight, ...
                                m2hXLimits, m2hYLimits, ...
                                m2hYTickLocs, colorMapPharm), ...
            exampleLabelsIpscr, outFoldersIpscr);
end
if plotM2hTauh
    for iMode = 1:numel(tauhModesAll)
        cellfun(@(x, y) plot_m2h(x, y, figure06Dir, figTypes, ...
                                m2hTauhWidth, m2hFigHeight, ...
                                m2hXLimits, m2hYLimits, ...
                                m2hYTickLocs, colorMapPharm), ...
                exampleLabelsModeAll{iMode}, outFoldersModeAll{iMode});
    end
end
if plotM2hGabab
    for iSheet = 1:numel(gababIpscSheetBases)
        cellfun(@(x, y) plot_m2h(x, y, figure07Dir, figTypes, ...
                                m2hGababWidth, m2hFigHeight, ...
                                m2hXLimits, m2hYLimits, ...
                                m2hYTickLocs, colorMapVary), ...
                exampleLabelsVaryAll{iSheet}, outFoldersVaryAll{iSheet});
    end
end

%% Plot voltage in soma against m2hDiff in dendrite 2
if plotVoltageVsOpdFig5
    cellfun(@(x, y) plot_voltage_vs_opd(x, y, figure05Dir, figTypes, ...
                    voltageVsOpdFig5FigWidth, voltageVsOpdFig5FigHeight, ...
                    voltageVsOpdTimeLimits1, voltageVsOpdSiMs, ...
                    voltageVsOpdFig5XLimits, voltageVsOpdFig5YLimits, ...
                    voltageVsOpdFig5YTickLocs, voltageVsOpdFig5ToAnnotate, ...
                    colorMapPharm, 'voltageVsOpd1'), ...
            exampleLabelsIpscr, outFoldersIpscr);
end
if plotVoltageVsOpdTauh
    for iMode = 1:numel(tauhModesAll)
        cellfun(@(x, y) plot_voltage_vs_opd(x, y, figure06Dir, figTypes, ...
                    voltageVsOpdTauhFigWidth, voltageVsOpdTauhFigHeight, ...
                    voltageVsOpdTimeLimits1, voltageVsOpdSiMs, ...
                    voltageVsOpdTauhXLimits, voltageVsOpdTauhYLimits, ...
                    voltageVsOpdTauhYTickLocs, voltageVsOpdTauhToAnnotate, ...
                    colorMapPharm, 'voltageVsOpd0'), ...
                exampleLabelsModeAll{iMode}, outFoldersModeAll{iMode});
    end
end
if plotVoltageVsOpdGabab
    for iSheet = 1:numel(gababIpscSheetBases)
        cellfun(@(x, y) plot_voltage_vs_opd(x, y, figure07Dir, figTypes, ...
                    voltageVsOpdGababFigWidth, voltageVsOpdGababFigHeight, ...
                    voltageVsOpdTimeLimits1, voltageVsOpdSiMs, ...
                    voltageVsOpdGababXLimits, voltageVsOpdGababYLimits, ...
                    voltageVsOpdGababYTickLocs, voltageVsOpdGababToAnnotate, ...
                    colorMapVary, 'voltageVsOpd1'), ...
                exampleLabelsVaryAll{iSheet}, outFoldersVaryAll{iSheet});
    end
end

%% Make movies of somatic voltage against m2hDiff in dendrite 2
if createPlotMovieFig5
    cellfun(@(x, y) plot_voltage_vs_opd(x, y, figure05Dir, figTypes, ...
                    voltageVsOpdFig5FigWidth, voltageVsOpdFig5FigHeight, ...
                    voltageVsOpdTimeLimits2, voltageVsOpdSiMs, ...
                    voltageVsOpdFig5XLimits, voltageVsOpdFig5YLimits, ...
                    voltageVsOpdFig5YTickLocs, voltageVsOpdFig5ToAnnotate, ...
                    colorMapPharm, 'voltageVsOpd2', ...
                    traceLabelsFig5, legendLocationFig5), ...
            exampleLabelsIpscr, outFoldersIpscr);
end
if createPlotMovieTauh
    for iMode = 1:numel(tauhModesAll)
        cellfun(@(x, y) plot_voltage_vs_opd(x, y, figure06Dir, figTypes, ...
                    voltageVsOpdTauhFigWidth, voltageVsOpdTauhFigHeight, ...
                    voltageVsOpdTimeLimits2, voltageVsOpdSiMs, ...
                    voltageVsOpdTauhXLimits, voltageVsOpdTauhYLimits, ...
                    voltageVsOpdTauhYTickLocs, voltageVsOpdTauhToAnnotate, ...
                    colorMapPharm, 'voltageVsOpd2', ...
                    traceLabelsTauh, legendLocationTauh), ...
                exampleLabelsModeAll{iMode}, outFoldersModeAll{iMode});
    end
end
if createPlotMovieGabab
    for iSheet = 1:numel(gababIpscSheetBases)
        cellfun(@(x, y) plot_voltage_vs_opd(x, y, figure07Dir, figTypes, ...
                    voltageVsOpdGababFigWidth, voltageVsOpdGababFigHeight, ...
                    voltageVsOpdTimeLimits2, voltageVsOpdSiMs, ...
                    voltageVsOpdGababXLimits, voltageVsOpdGababYLimits, ...
                    voltageVsOpdGababYTickLocs, voltageVsOpdGababToAnnotate, ...
                    colorMapVary, 'voltageVsOpd2', ...
                    traceLabelsGabab, legendLocationGabab), ...
                exampleLabelsVaryAll{iSheet}, outFoldersVaryAll{iSheet});
    end
end

%% Make movies comparing open discrepancy slope with voltage slope
if createPhasePlotOnlyMovieFig5
    cellfun(@(x, y) plot_voltage_vs_opd(x, y, figure05Dir, figTypes, ...
                    voltageVsOpdFig5FigWidth, voltageVsOpdFig5FigHeight, ...
                    voltageVsOpdTimeLimits2, voltageVsOpdSiMs, ...
                    voltageVsOpdFig5XLimits, voltageVsOpdFig5YLimits, ...
                    voltageVsOpdFig5YTickLocs, voltageVsOpdFig5ToAnnotate, ...
                    colorMapPharm, 'voltageVsOpd3', ...
                    traceLabelsFig5, legendLocationFig5), ...
            exampleLabelsIpscr, outFoldersIpscr);
end
if createPhasePlotOnlyMovieTauh
    for iMode = 1:numel(tauhModesAll)
        cellfun(@(x, y) plot_voltage_vs_opd(x, y, figure06Dir, figTypes, ...
                    voltageVsOpdTauhFigWidth, voltageVsOpdTauhFigHeight, ...
                    voltageVsOpdTimeLimits2, voltageVsOpdSiMs, ...
                    voltageVsOpdTauhXLimits, voltageVsOpdTauhYLimits, ...
                    voltageVsOpdTauhYTickLocs, voltageVsOpdTauhToAnnotate, ...
                    colorMapPharm, 'voltageVsOpd3', ...
                    traceLabelsTauh, legendLocationTauh), ...
                exampleLabelsModeAll{iMode}, outFoldersModeAll{iMode});
    end
end
if createPhasePlotOnlyMovieGabab
    for iSheet = 1:numel(gababIpscSheetBases)
        cellfun(@(x, y) plot_voltage_vs_opd(x, y, figure07Dir, figTypes, ...
                    voltageVsOpdGababFigWidth, voltageVsOpdGababFigHeight, ...
                    voltageVsOpdTimeLimits2, voltageVsOpdSiMs, ...
                    voltageVsOpdGababXLimits, voltageVsOpdGababYLimits, ...
                    voltageVsOpdGababYTickLocs, voltageVsOpdGababToAnnotate, ...
                    colorMapVary, 'voltageVsOpd3', ...
                    traceLabelsGabab, legendLocationGabab), ...
                exampleLabelsVaryAll{iSheet}, outFoldersVaryAll{iSheet});
    end
end

%% Simulate when there is no T current in the soma
if simulateNoITSoma
    cd(figure05Dir);
    check_dir(outFoldersNoITSoma);
    cellfun(@(x, y, z) simulate_ipscr(x, y, z, 0, dataModeIpscr, ...
                                    rowmodeIpscr, attemptNumberIpscr, ...
                                    newParamsNoITSoma), ...
            exampleLabelsNoITSoma, exampleParamPaths, outFoldersNoITSoma);
end

%% Archive all scripts for this run
if archiveScriptsFlag
    if simulateIpscr || plotAllVoltages || plotAllTotalCurrents || ...
            plotAllComponentCurrents || plotDend2ITproperties || ...
            plotM2hFig5 || plotVoltageVsOpdFig5 || createPlotMovieFig5 || ...
            createPhasePlotOnlyMovieFig5 || simulateNoITSoma
        archive_dependent_scripts(mfilename, 'OutFolder', figure05Dir);
    end
    if simulateTauhModes || plotSomaVoltage || plotM2hTauh || ...
            plotVoltageVsOpdTauh || createPlotMovieTauh || ...
            createPhasePlotOnlyMovieTauh
        archive_dependent_scripts(mfilename, 'OutFolder', figure06Dir);
    end
    if simulateIpscVariation || plotEssential || plotM2hGabab || ...
            plotVoltageVsOpdGabab || createPlotMovieGabab || ...
            createPhasePlotOnlyMovieGabab
        archive_dependent_scripts(mfilename, 'OutFolder', figure07Dir);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function simulate_ipscr(label, neuronParamsFile, outFolder, ...
                        tauhMode, dataMode, rowmode, attemptNumber, newParams)

if nargin < 8
    newParams = struct.empty;
end

% Simulate
m3ha_neuron_run_and_analyze(neuronParamsFile, ...
                        'OutFolder', outFolder, 'Prefix', label, ...
                        'BuildMode', 'active', 'SimMode', 'active', ...
                        'TauhMode', tauhMode, 'DataMode', dataMode, ...
                        'ColumnMode', 1, 'Rowmode', rowmode, ...
                        'AttemptNumber', attemptNumber, ...
                        'PlotAllFlag', false, 'PlotIndividualFlag', true, ...
                        'SaveSimOutFlag', true, 'NewParams', newParams);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function simulate_variation (label, neuronParamsFile, outFolder, gababParams)

% Simulate
m3ha_neuron_run_and_analyze(neuronParamsFile, ...
                        'OutFolder', outFolder, 'Prefix', label, ...
                        'BuildMode', 'active', 'SimMode', 'active', ...
                        'PlotAllFlag', false, 'PlotOverlappedFlag', true, ...
                        'SaveSimOutFlag', true, 'NoRealDataFlag', true, ...
                        gababParams);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_overlapped (expStr, directory, plotType, outFolder, figTypes, ...
                            figWidth, figHeight, xLimits, yLimits, ...
                            colorMap, yTickLocs)

% Create figure names
figPathBase = fullfile(outFolder, [expStr, '_', plotType]);
figPathBaseOrig = [figPathBase, '_orig'];

% Create the figure
fig = set_figure_properties('AlwaysNew', true);

% Plot traces
handles = ...
    m3ha_plot_simulated_traces('Directory', directory, 'ExpStr', expStr, ...
                'PlotType', plotType, 'FigHandle', fig, ...
                'FigTitle', 'suppress', 'XLabel', 'suppress', ...
                'XLimits', xLimits, 'YLimits', yLimits, ...
                'ColorMap', colorMap);

update_figure_for_corel(fig, 'YTickLocs', yTickLocs);

% Add a threshold line
switch plotType
    case {'essential', 'somaVoltage', 'dend2ITproperties'}
        % Must be consistent with m3ha_simulate_population.m
        opdThreshold = 1e-2;

        % Find the subplot of interest
        subPlots = handles.subPlots;
        switch plotType
            case 'essential'
                subplot(subPlots(5));
            case 'somaVoltage'
                subplot(subPlots(2));
            case 'dend2ITproperties'
                subplot(subPlots(9));
        end
        hold on;

        % Add a threshold line
        plot_horizontal_line(opdThreshold, 'ColorMap', 'DarkGreen', ...
                                'LineStyle', ':', 'LineWidth', 1);
    otherwise
        % Do nothing
end

% Save original figure
save_all_figtypes(fig, figPathBaseOrig, 'png');

% Plot a scale bar
hold on
plot_scale_bar('x', 'XBarUnits', 'ms', 'XBarLength', 200, ...
                'XPosNormalized', 0.1, 'YPosNormalized', 0.8);

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                'Width', figWidth, 'Height', figHeight, ...
                'AlignSubplots', true, ...
                'RemoveXTicks', true, 'RemoveXRulers', true);

% Save the figure
save_all_figtypes(fig, figPathBase, figTypes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_m2h (expStr, directory, outFolder, figTypes, ...
                    figWidth, figHeight, xLimits, yLimits, ...
                    yTickLocs, colorMap)

% Create a figure name
figPathBaseVoltVsOpd = fullfile(outFolder, [expStr, '_m2h']);
figPathBaseM2hOrig = [figPathBaseVoltVsOpd, '_orig'];

% Create the figure
figM2h = set_figure_properties('AlwaysNew', true);

% Plot traces
m3ha_plot_simulated_traces('Directory', directory, 'ExpStr', expStr, ...
                            'PlotType', 'm2h', 'FigHandle', figM2h, ...
                            'FigTitle', 'suppress', ...
                            'XLimits', xLimits, 'YLimits', yLimits, ...
                            'ColorMap', colorMap);

% Save the figure
save_all_figtypes(figM2h, figPathBaseM2hOrig, 'png');

% Plot a scale bar
hold on
plot_scale_bar('x', 'XBarUnits', 'ms', 'XBarLength', 200, ...
                'XPosNormalized', 0.1, 'YPosNormalized', 0.8);

% Update figure for CorelDraw
update_figure_for_corel(figM2h, 'Units', 'centimeters', ...
                'Width', figWidth, 'Height', figHeight, ...
                'AlignSubplots', true, 'YTickLocs', yTickLocs, ...
                'RemoveXTicks', true, 'RemoveXRulers', true);

% Save the figure
save_all_figtypes(figM2h, figPathBaseVoltVsOpd, figTypes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_voltage_vs_opd (expStr, directory, outFolder, figTypes, ...
                                figWidth, figHeight, timeLimits, siMs, ...
                                xLimits, yLimits, yTickLocs, ...
                                toAnnotate, colorMap, plotType, ...
                                traceLabels, legendLocation)

% Hard-coded constants
MS_PER_S = 1000;

% Hard-coded parameters
nSamples = floor(diff(timeLimits) / siMs);
tVecToMatch = create_time_vectors(nSamples, 'TimeStart', timeLimits(1), ...
                            'SamplingIntervalMs', siMs, 'TimeUnits', 'ms');
plotMarkerSize = 3;
% fiSeconds = (siMs / MS_PER_S) * 100;            % play at 1/100 x speed
fiSeconds = (siMs / MS_PER_S) * 50;            % play at 1/50 x speed

% Create figure path bases
figPathBaseVoltVsOpd = fullfile(outFolder, [expStr, '_voltageVsOpd']);
figPathBaseVoltVsOpdOrig = [figPathBaseVoltVsOpd, '_orig'];
figPathBaseVoltVsOpdCompressed = [figPathBaseVoltVsOpd, '_compressed'];
figPathBaseVoltVsOpdMovie = [figPathBaseVoltVsOpd, '_movie'];
figPathBaseVoltVsOpdAligned = [figPathBaseVoltVsOpd, '_movie_aligned'];
figPathBaseVoltVsOpdSimple = [figPathBaseVoltVsOpd, '_movie_simple'];
figPathBaseConcavityVsSlope = ...
    fullfile(outFolder, [expStr, '_concavityVsSlope']);
figPathBaseConcavityVsSlopeMovie = ...
    [figPathBaseConcavityVsSlope, '_movie'];
figPathBaseConcavityVsSlopeAligned = ...
    [figPathBaseConcavityVsSlope, '_movie_aligned'];
figPathBaseConcavityVsSlopeSimple = ...
    [figPathBaseConcavityVsSlope, '_movie_simple'];

% Create and save plot
switch plotType
case {'voltageVsOpd0', 'voltageVsOpd1'}
    % Create the figure
    figVoltVsOpdOrig = set_figure_properties('AlwaysNew', true);

    % Plot traces
    handles = ...
        m3ha_plot_simulated_traces('Directory', directory, 'ExpStr', expStr, ...
                            'PlotType', plotType, ...
                            'FigHandle', figVoltVsOpdOrig, ...
                            'TimeLimits', timeLimits, ...
                            'XLimits', xLimits, 'YLimits', yLimits, ...
                            'ColorMap', colorMap, 'LineStyle', 'none', ...
                            'tVecs', tVecToMatch, 'Marker', '.');

    % Save the figure
    save_all_figtypes(figVoltVsOpdOrig, figPathBaseVoltVsOpdOrig, 'png');

    % Update figure for CorelDraw
    update_figure_for_corel(figVoltVsOpdOrig, 'AlignSubplots', true);
    update_figure_for_corel(figVoltVsOpdOrig, 'Units', 'centimeters', ...
                            'Width', figWidth, 'Height', figHeight, ...
                            'RemoveLabels', true, 'RemoveTitles', true, ...
                            'YTickLocs', yTickLocs, ...
                            'PlotMarkerSize', plotMarkerSize);
    if ~toAnnotate
        update_figure_for_corel(figVoltVsOpdOrig, 'RemoveCircles', true);
    end

    % Plot a scale bar
    subplot(handles.ax(1)); hold on;
    plot_scale_bar('x', 'XBarUnits', 'ms', 'XBarLength', 200, ...
                    'XPosNormalized', 0.1, 'YPosNormalized', 0.8);

    % Save the figure
    save_all_figtypes(figVoltVsOpdOrig, figPathBaseVoltVsOpd, figTypes);

    % Update figure for CorelDraw
    update_figure_for_corel(figVoltVsOpdOrig, 'Units', 'centimeters', ...
                            'Width', figWidth, 'Height', figHeight/2);

    % Save the figure
    save_all_figtypes(figVoltVsOpdOrig, figPathBaseVoltVsOpdCompressed, ...
                        figTypes);
case {'voltageVsOpd2', 'voltageVsOpd3'}    
    % Create the figure
    figMovie1 = set_figure_properties('AlwaysNew', true);

    % Decide on figure title
    if contains(expStr, 'D101310_aft_ipscr')
        figTitle = ['Voltage and T channel gating trajectories are ', ...
                    'different following distinct GABA_B IPSCs'];
    elseif contains(expStr, 'dual_vary_tau')
        figTitle = ['Voltage and T channel gating trajectories ', ...
                    'change as {\it sim}Dual-Block waveform time ', ...
                    'constants are varied'];
    else
        figTitle = '';
    end

    % Plot traces
    handles = m3ha_plot_simulated_traces('Directory', directory, ...
                    'ExpStr', expStr, 'PlotType', plotType, ...
                    'FigHandle', figMovie1, 'TimeLimits', timeLimits, ...
                    'XLimits', xLimits, 'YLimits', yLimits, ...
                    'ColorMap', colorMap, 'LineStyle', 'none', ...
                    'tVecs', tVecToMatch, 'Marker', '.', ...
                    'PlotSelected', true, 'FigTitle', figTitle);

    % Plot legend on first subplot
    if ~strcmp(legendLocation, 'suppress')
        ax1 = handles.ax(1);
        plots = handles.ax1Stuff.traces.plotsData;
        legend(ax1, plots, traceLabels, 'Location', legendLocation);
    end

    % Decide on movie file base
    switch plotType
    case 'voltageVsOpd2'
        fileBase = figPathBaseVoltVsOpdMovie;
    case 'voltageVsOpd3'
        fileBase = figPathBaseConcavityVsSlopeMovie;
    end

    % Create movie
    create_plot_movie(figMovie1, fiSeconds, 'FileBase', fileBase);

    % Create the figure
    figMovie2 = set_figure_properties('AlwaysNew', true);

    % Plot traces
    handles = m3ha_plot_simulated_traces('Directory', directory, ...
                    'ExpStr', expStr, 'PlotType', plotType, ...
                    'FigHandle', figMovie2, 'TimeLimits', timeLimits, ...
                    'XLimits', xLimits, 'YLimits', yLimits, ...
                    'ColorMap', colorMap, 'LineStyle', 'none', ...
                    'tVecs', tVecToMatch, 'Marker', '.', ...
                    'PlotSelected', true, 'FigTitle', figTitle);

    % Plot legend on first subplot
    if ~strcmp(legendLocation, 'suppress')
        ax1 = handles.ax(1);
        plots = handles.ax1Stuff.traces.plotsData;
        legend(ax1, plots, traceLabels, 'Location', legendLocation);
    end

    % Decide on movie file base
    switch plotType
    case 'voltageVsOpd2'
        fileBase = figPathBaseVoltVsOpdAligned;
    case 'voltageVsOpd3'
        fileBase = figPathBaseConcavityVsSlopeAligned;
    end

    if contains(expStr, 'D101310_aft_ipscr') || contains(expStr, 'dual_vary_tau')
        textArrows = handles.textArrows;
        if contains(expStr, 'D101310_aft_ipscr')
            frameNumForTexts = 740 * ones(size(textArrows));
        else
            frameNumForTexts = 680 * ones(size(textArrows));
        end
        objectFrameStartPair = {textArrows, frameNumForTexts};
    else
        objectFrameStartPair = {};
    end
    
    % Create movie
    create_plot_movie(figMovie2, fiSeconds, 'FileBase', fileBase, ...
                        'AlignToSelected', true, ...
                        'ObjectFrameStartPair', objectFrameStartPair);

    % Create the figure
    figMovie3 = set_figure_properties('AlwaysNew', true);

    % Decide on figure title
    if contains(expStr, 'D101310_aft_ipscr')
        figTitle = ['Voltage and T channel gating trajectories are ', ...
                    'different following {\it sim}GAT3-Block ', ...
                    'vs {\it sim}Dual-Block'];
    elseif contains(expStr, 'dual_vary_tau')
        figTitle = ['Voltage and T channel gating trajectories ', ...
                'are different for two {\it sim}Dual-Block time constants'];
    else
        figTitle = '';
    end

    % Decide on colors to keep
    if contains(expStr, 'D101310_aft_ipscr')
        nColors = 4;
        iColorsToKeep = 3:4;
        darkPercentage = [];
    elseif contains(expStr, 'dual_vary_tau')
        % nColors = 17;
        % iColorsToKeep = 15:16;
        nColors = 9;
        iColorsToKeep = 7:8;
        darkPercentage = 25;
    end

    % Decide on the color map
    colorMapOrig = decide_on_colormap(colorMap, nColors, ...
                        'DarkPercentage', darkPercentage, ...
                        'ForceCellOutput', true);

    % Plot traces
    handles = m3ha_plot_simulated_traces('Directory', directory, ...
                    'ExpStr', expStr, 'PlotType', plotType, ...
                    'FigHandle', figMovie3, 'TimeLimits', timeLimits, ...
                    'XLimits', xLimits, 'YLimits', yLimits, ...
                    'ColorMap', colorMapOrig, 'LineStyle', 'none', ...
                    'tVecs', tVecToMatch, 'Marker', '.', ...
                    'PlotSelected', true, 'FigTitle', figTitle);

    % Plot legend on first subplot
    if ~strcmp(legendLocation, 'suppress')
        ax1 = handles.ax(1);
        plots = handles.ax1Stuff.traces.plotsData;
        legend(ax1, plots, traceLabels, 'Location', legendLocation);
    end

    % Decide on movie file base
    switch plotType
    case 'voltageVsOpd2'
        fileBase = figPathBaseVoltVsOpdSimple;
    case 'voltageVsOpd3'
        fileBase = figPathBaseConcavityVsSlopeSimple;
    end

    % Remove line objects
    lines = findobj(figMovie3, 'Type', 'Line');
    colors = extract_fields(lines, 'Color');
    if contains(expStr, 'D101310_aft_ipscr') || contains(expStr, 'dual_vary_tau')
        colorMapFaded = decide_on_colormap(colorMapOrig, nColors, ...
                            'FadePercentage', 70, 'ForceCellOutput', true);
        toRemove = ~ismember_custom(colors, [colorMapOrig(iColorsToKeep); ...
                                            colorMapFaded(iColorsToKeep)]);

        delete(lines(toRemove));

        if contains(expStr, 'dual_vary_tau')
            % Remove color bar
            colorBar = findobj(gcf, 'Type', 'ColorBar');
            delete(colorBar);

            % Create new legend
            % traces = handles.ax5Stuff.traces.plotsData;
            traces = handles.ax1Stuff.traces.plotsData;
            legend(traces(iColorsToKeep), ...
                    {'\tau = 2.0 sec', '\tau = 2.3 sec'}, ...
                    'Location', 'southeast');
        end

        textArrows = handles.textArrows;
        if contains(expStr, 'D101310_aft_ipscr')
            frameNumForTexts = 740 * ones(size(textArrows));
        else
            frameNumForTexts = 660 * ones(size(textArrows));
        end
        objectFrameStartPair = {textArrows, frameNumForTexts};

        % Draw now
        drawnow;

        % Create movie
        create_plot_movie(figMovie3, fiSeconds, 'FileBase', fileBase, ...
                            'AlignToSelected', true, ...
                            'ObjectFrameStartPair', objectFrameStartPair);
    end
otherwise
    error('plotType unrecognized!')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Set pharm legend to autoupdate
if contains(expStr, 'D101310_aft_ipscr')
    lgds = findobj(gcf, 'Type', 'Legend');
    set(lgds(1), 'AutoUpdate', 'on');
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
