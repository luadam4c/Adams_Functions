% m3ha_plot_figure05.m
%% Plots Figure 05 & Figure 06 for the GAT Blocker paper
%
% Requires:
%       cd/archive_dependent_scripts.m
%       cd/check_dir.m
%       cd/create_labels_from_numbers.m
%       cd/create_time_vectors.m
%       cd/extract_fileparts.m
%       cd/find_matching_files.m
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
% 2020-04-13 Added plotVoltageVsOpd
% 2020-05-04 Added createPlotMovie

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

plotM2h = false; %true;

plotVoltageVsOpd = true;
createPlotMovieFig5 = true;
createPlotMovieFig6 = false;

simulateNoITSoma = false; %true;

archiveScriptsFlag = true;

% Directories
parentDirectory = fullfile('/media', 'adamX', 'm3ha');
figure02Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure02');
figure05Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure05');
figure06Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure06');
gababIpscDir = figure06Dir;
fitDirectory = fullfile(parentDirectory, 'optimizer4gabab');

% Files
sweepInfoFile = 'dclampdatalog_take4.csv';
datalogPath = fullfile(figure02Dir, sweepInfoFile);
paramFileSuffix = 'params';

% Analysis settings
% exampleCellNames = {'D101310'; 'C101210'};
% exampleCellNames = {'D101310'; 'M101210'};
exampleCellNames = {'D101310'};
% exampleCellNames = {'D101310'; 'G101310'};

% Must be consistent with m3ha_compute_gabab_ipsc.m
gababIpscSheetBases = {'gababipsc_gat3_vary_amp2', ...
                        'gababipsc_gat3_vary_amp', ...
                        'gababipsc_dual_vary_amp', ...
                        'gababipsc_gat3_vary_tau', ...
                        'gababipsc_dual_vary_tau', ...
                        'gababipsc_vary_dual_to_gat3_to_gat1', ...
                        'gababipsc_original'};

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
% tauhModesAll = 4:5;
% tauhModesAll = 1:5;
% tauhModesAll = 6:7;
tauhModesAll = 1:7;

newParamsNoITSoma = {'pcabarITSoma', 0};

% The following must be consistent with singleneuron4compgabab.hoc
timeToStabilize = 2000;         % padded time (ms) to make sure initial value 
                                %   of simulations are stabilized

% Plot settings
colorMapPharm = [];                 % use m3ha default
colorMapVary = @jet;                % rainbow colors

overlappedFigWidth = 4.7; %5.7;
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
                            [1e-8, 1e0], [1e-8, 1e0], [1e-1, 1e2]};
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
                            [1e-6, 1e-2], [1e-6, 1e-2], [1e0, 1e1]};
somaVoltageYTickLocs = {-90:20:-50, [1e-6, 1e-2]};
essentialYTickLocs = {-90:20:-50, 0:5:10, -0.4:0.2:0, ...
                            -15:5:0, [1e-6, 1e-2]};
m2hFigWidth = 4.7; %5.7;
m2hFigHeight = 3;
m2hXLimits = timeToStabilize + [800, 2800]; %[800, 2000];
m2hYLimits = [1e-7, 1e0];
m2hYTickLocs = [1e-5, 1e-3, 1e-1];

voltageVsOpdTimeLimits1 = timeToStabilize + [800, 2800];
voltageVsOpdTimeLimits2 = timeToStabilize + [1000, 2000];
voltageVsOpdSiMs = 1;
voltageVsOpdFig5FigWidth = 5.5 * 2;
voltageVsOpdFig5FigHeight = 5 * 2;
voltageVsOpdFig5XLimits = [1e-7, 1e0];
voltageVsOpdFig5YLimits = [-95, -45];
voltageVsOpdFig5YTickLocs = [];
voltageVsOpdFig5ToAnnotate = true;

voltageVsOpdFig6FigWidth = 4.7 * 2;
voltageVsOpdFig6FigHeight = 2.2 * 2;
voltageVsOpdFig6XLimits = [1e-7, 1e0];
voltageVsOpdFig6YLimits = [-95, -45];
voltageVsOpdFig6YTickLocs = [];
voltageVsOpdFig6ToAnnotate = false;

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
        plotSomaVoltage || plotM2h || plotVoltageVsOpd || simulateNoITSoma
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
    outFoldersVaryAll = cellfun(@(x) fullfile(figure06Dir, x), ...
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
    m3ha_compute_gabab_ipsc(gababIpscDir);
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
                    overlappedFigWidth, ...
                    numel(allVoltagesYLimits) * overlappedFigHeightPerRow, ...
                    overlappedXLimits, allVoltagesYLimits, ...
                    colorMapPharm, allVoltagesYTickLocs), ...
            exampleLabelsIpscr, outFoldersIpscr);
end

%% Plot all currents
if plotAllTotalCurrents
    cellfun(@(x, y) plot_overlapped(x, y, 'allTotalCurrents', ...
                    figure05Dir, figTypes, ...
                    overlappedFigWidth, ...
                    numel(allTotalCurrentsYLimits) * overlappedFigHeightPerRow, ...
                    overlappedXLimits, allTotalCurrentsYLimits, ...
                    colorMapPharm, allTotalCurrentsYTickLocs), ...
            exampleLabelsIpscr, outFoldersIpscr);
end

%% Plot component currents
if plotAllComponentCurrents
    cellfun(@(x, y) plot_overlapped(x, y, 'allComponentCurrents', ...
                    figure05Dir, figTypes, ...
                    overlappedFigWidth, ...
                    numel(allComponentCurrentsYLimits) * overlappedFigHeightPerRow, ...
                    overlappedXLimits, allComponentCurrentsYLimits, ...
                    colorMapPharm, allComponentCurrentsYTickLocs), ...
            exampleLabelsIpscr, outFoldersIpscr);
end

%% Plot all T channel properties
if plotDend2ITproperties
    cellfun(@(x, y) plot_overlapped(x, y, 'dend2ITproperties', ...
                    figure05Dir, figTypes, ...
                    overlappedFigWidth, ...
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
                        overlappedFigWidth, ...
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
                        figure06Dir, figTypes, ...
                        overlappedFigWidth, ...
                        numel(essentialYLimitsTemp) * overlappedFigHeightPerRow, ...
                        overlappedXLimits, essentialYLimitsTemp, ...
                        colorMapVary, essentialYTicksTemp), ...
                exampleLabelsVaryAll{iSheet}, outFoldersVaryAll{iSheet});
    end
end

%% Plot m2h in dendrite 2 against its steady state
if plotM2h
    cellfun(@(x, y) plot_m2h(x, y, figure05Dir, figTypes, ...
                                m2hFigWidth, m2hFigHeight, ...
                                m2hXLimits, m2hYLimits, ...
                                m2hYTickLocs, colorMapPharm), ...
            exampleLabelsIpscr, outFoldersIpscr);
    for iMode = 1:numel(tauhModesAll)
        cellfun(@(x, y) plot_m2h(x, y, figure06Dir, figTypes, ...
                                m2hFigWidth, m2hFigHeight, ...
                                m2hXLimits, m2hYLimits, ...
                                m2hYTickLocs, colorMapPharm), ...
                exampleLabelsModeAll{iMode}, outFoldersModeAll{iMode});
    end

    for iSheet = 1:numel(gababIpscSheetBases)
        cellfun(@(x, y) plot_m2h(x, y, figure06Dir, figTypes, ...
                                m2hFigWidth, m2hFigHeight, ...
                                m2hXLimits, m2hYLimits, ...
                                m2hYTickLocs, colorMapVary), ...
                exampleLabelsVaryAll{iSheet}, outFoldersVaryAll{iSheet});
    end
end

%% Plot voltage in soma against m2hDiff in dendrite 2
if plotVoltageVsOpd
    cellfun(@(x, y) plot_voltage_vs_opd(x, y, figure05Dir, figTypes, ...
                    voltageVsOpdFig5FigWidth, voltageVsOpdFig5FigHeight, ...
                    voltageVsOpdTimeLimits1, voltageVsOpdSiMs, ...
                    voltageVsOpdFig5XLimits, voltageVsOpdFig5YLimits, ...
                    voltageVsOpdFig5YTickLocs, voltageVsOpdFig5ToAnnotate, ...
                    colorMapPharm, createPlotMovieFig5, voltageVsOpdTimeLimits2), ...
            exampleLabelsIpscr, outFoldersIpscr);
%{
    for iMode = 1:numel(tauhModesAll)
        cellfun(@(x, y) plot_voltage_vs_opd(x, y, figure06Dir, figTypes, ...
                    voltageVsOpdFig6FigWidth, voltageVsOpdFig6FigHeight, ...
                    voltageVsOpdTimeLimits, voltageVsOpdSiMs, ...
                    voltageVsOpdFig6XLimits, voltageVsOpdFig6YLimits, ...
                    voltageVsOpdFig6YTickLocs, voltageVsOpdFig6ToAnnotate, ...
                    colorMapPharm, createPlotMovieFig6), ...
                exampleLabelsModeAll{iMode}, outFoldersModeAll{iMode});
    end
    for iSheet = 1:numel(gababIpscSheetBases)
        cellfun(@(x, y) plot_voltage_vs_opd(x, y, figure06Dir, figTypes, ...
                    voltageVsOpdFig6FigWidth, voltageVsOpdFig6FigHeight, ...
                    voltageVsOpdTimeLimits, voltageVsOpdSiMs, ...
                    voltageVsOpdFig6XLimits, voltageVsOpdFig6YLimits, ...
                    voltageVsOpdFig6YTickLocs, voltageVsOpdFig6ToAnnotate, ...
                    colorMapVary, createPlotMovieFig6), ...
                exampleLabelsVaryAll{iSheet}, outFoldersVaryAll{iSheet});
    end
%}
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
    archive_dependent_scripts(mfilename, 'OutFolder', figure05Dir);
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
                                toAnnotate, colorMap, ...
                                createPlotMovie, timeLimitsMovie)

% Hard-coded constants
MS_PER_S = 1000;

% Hard-coded parameters
nSamples = floor(diff(timeLimits) / siMs);
tVecToMatch = create_time_vectors(nSamples, 'TimeStart', timeLimits(1), ...
                            'SamplingIntervalMs', siMs, 'TimeUnits', 'ms');
plotMarkerSize = 3;
fiSeconds = (siMs / MS_PER_S) * 100;            % play at 1/100 x speed

% Create a figure name
figPathBaseVoltVsOpd = fullfile(outFolder, [expStr, '_voltageVsOpd']);
figPathBaseVoltVsOpdOrig = [figPathBaseVoltVsOpd, '_orig'];
figPathBaseVoltVsOpdCompressed = [figPathBaseVoltVsOpd, '_compressed'];
figPathBaseVoltVsOpdMovie = [figPathBaseVoltVsOpd, '_movie'];

% Create the figure
figVoltVsOpdOrig = set_figure_properties('AlwaysNew', true);

% Plot traces
handles = ...
    m3ha_plot_simulated_traces('Directory', directory, 'ExpStr', expStr, ...
                        'PlotType', 'voltageVsOpd1', ...
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
save_all_figtypes(figVoltVsOpdOrig, figPathBaseVoltVsOpdCompressed, figTypes);

% Create a plot movie if requested
if createPlotMovie
    % Create the figure
    figMovie = set_figure_properties('AlwaysNew', true);

    % Plot traces
    handles = m3ha_plot_simulated_traces('Directory', directory, ...
                    'ExpStr', expStr, 'PlotType', 'voltageVsOpd2', ...
                    'FigHandle', figMovie, 'TimeLimits', timeLimitsMovie, ...
                    'XLimits', xLimits, 'YLimits', yLimits, ...
                    'ColorMap', colorMap, 'LineStyle', 'none', ...
                    'tVecs', tVecToMatch, 'Marker', '.');

    % Create movie
    create_plot_movie(figMovie, fiSeconds, ...
                        'FileBase', figPathBaseVoltVsOpdMovie);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
