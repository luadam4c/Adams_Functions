% m3ha_plot_figure04.m
%% Plots Figure 04 for the GAT Blocker paper
%
% Requires:
%       cd/update_neuron_scripts.m
% TODO
%       cd/check_dir.m
%       cd/extract_fileparts.m
%       cd/find_matching_files.m
%       cd/m3ha_load_sweep_info.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/plot_scale_bar.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/update_figure_for_corel.m

% File History:
% 2019-12-29 Modified from m3ha_plot_figure03.m

%% Hard-coded parameters
% Flags
updateScripts = false; %true;
simulateIpscr = false; %true;
plotOverlapped = true;
plotM2h = true;

% Directories
parentDirectory = fullfile('/media', 'adamX', 'm3ha');
figure02Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure02');
figure04Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure04');
fitDirectory = fullfile(parentDirectory, 'optimizer4gabab');

% Files
sweepInfoFile = 'dclampdatalog_take4.csv';
datalogPath = fullfile(figure02Dir, sweepInfoFile);
paramFileSuffix = 'params';

% Analysis settings
exampleCellNames = {'D101310'; 'C101210'};
% exampleCellNames = {'D101310'};
% exampleCellNames = {'C101210'};

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
attemptNumberIpscr = 1;             % attempt number for IPSC response
                                    %   1 - Use 4 traces @ 200% gIncr 
                                    %           for this data mode
                                    %   2 - Use all traces @ 200% gIncr 
                                    %           for this data mode
                                    %   3 - Use all traces for this data mode
                                    %   4 - Use 1 trace for each pharm x gIncr 
                                    %           pair for this data mode
                                    %   5 - Use 4 traces @ 400% gIncr 
                                    %       for this data mode

% Plot settings
overlappedFigWidth = 6;
overlappedFigHeightPerSubplot = 1.5;
overlappedXLimits = [2800, 4000];
overlappedYLimits = [];
m2hFigWidth = 6;
m2hFigHeight = 3;
m2hXLimits = [2800, 4500];
m2hYLimits = [1e-5, 0.1];

figTypes = {'png', 'epsc2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Make sure NEURON scripts are up to date in figure04Dir
if updateScripts
    update_neuron_scripts(fitDirectory, figure04Dir);
end

%% Load sweep info
% Read from datalogPath
swpInfo = m3ha_load_sweep_info('Directory', figure02Dir);

%% Find NEURON parameter tables
if simulateIpscr || plotOverlapped || plotM2h
    % Find NEURON parameter tables
    [~, exampleParamPaths] = ...
        find_matching_files(exampleCellNames, 'Directory', figure04Dir, ...
                            'Suffix', paramFileSuffix, 'Extension', 'csv', ...
                            'Recursive', false);

    % Extract file bases
    exampleParamFileBases = extract_fileparts(exampleParamPaths, 'base');

    % Extract example labels
    exampleLabels = extractBefore(exampleParamFileBases, ...
                                    ['_', paramFileSuffix]);

    % Update labels for each type of simulation
    exampleLabelsIpscr = strcat(exampleLabels, '_ipscr');

    % Create and check output folders
    outFoldersIpscr = fullfile(figure04Dir, exampleLabelsIpscr);
    check_dir(outFoldersIpscr);
end

%% Simulate IPSC responses
if simulateIpscr
    cellfun(@(x, y, z) simulate_ipscr(x, y, z, dataModeIpscr, rowmodeIpscr, ...
                                attemptNumberIpscr), ...
            exampleLabelsIpscr, exampleParamPaths, outFoldersIpscr);
end

%% Plot overlapped traces
if plotOverlapped
    % Plot all voltages
    cellfun(@(x, y) plot_overlapped(x, y, 'allvoltages', ...
                    figure04Dir, figTypes, ...
                    overlappedFigWidth, 7 * overlappedFigHeightPerSubplot, ...
                    overlappedXLimits, overlappedYLimits), ...
            exampleLabelsIpscr, outFoldersIpscr);

    % Plot all currents
    cellfun(@(x, y) plot_overlapped(x, y, 'allcurrents', ...
                    figure04Dir, figTypes, ...
                    overlappedFigWidth, 10 * overlappedFigHeightPerSubplot, ...
                    overlappedXLimits, overlappedYLimits), ...
            exampleLabelsIpscr, outFoldersIpscr);

    % Plot all T channel properties
    cellfun(@(x, y) plot_overlapped(x, y, 'dend2ITproperties', ...
                    figure04Dir, figTypes, ...
                    overlappedFigWidth, 5 * overlappedFigHeightPerSubplot, ...
                    overlappedXLimits, overlappedYLimits), ...
            exampleLabelsIpscr, outFoldersIpscr);
end

%% Plot m2h in dendrite 2 against its steady state
if plotM2h
    cellfun(@(x, y) plot_m2h(x, y, figure04Dir, figTypes, ...
                                m2hFigWidth, m2hFigHeight, ...
                                m2hXLimits, m2hYLimits), ...
            exampleLabelsIpscr, outFoldersIpscr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function simulate_ipscr(label, neuronParamsFile, outFolder, ...
                        dataMode, rowmode, attemptNumber)

% Simulate
m3ha_neuron_run_and_analyze(neuronParamsFile, ...
                        'OutFolder', outFolder, 'Prefix', label, ...
                        'BuildMode', 'active', 'SimMode', 'active', ...
                        'DataMode', dataMode, 'ColumnMode', 1, ...
                        'Rowmode', rowmode, 'AttemptNumber', attemptNumber, ...
                        'PlotAllFlag', false, 'PlotIndividualFlag', true, ...
                        'SaveSimOutFlag', true);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_overlapped(expStr, directory, plotType, outFolder, figTypes, ...
                            figWidth, figHeight, xLimits, yLimits)

% Create figure names
figPathBaseOrig = fullfile(outFolder, [expStr, '_', plotType, '_orig']);
figPathBase = fullfile(outFolder, [expStr, '_', plotType]);

% Create the figure
fig = set_figure_properties('AlwaysNew', true);

% Plot traces
m3ha_plot_simulated_traces('Directory', directory, 'ExpStr', expStr, ...
                'PlotType', plotType, 'FigHandle', fig, ...
                'FigTitle', 'suppress', 'XLabel', 'suppress', ...
                'XLimits', xLimits, 'YLimits', yLimits);

% Save the figure
save_all_figtypes(fig, figPathBaseOrig, 'png');

% Plot a scale bar
hold on
plot_scale_bar('x', 'XBarUnits', 'ms', 'XBarLength', 400, ...
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

function plot_m2h(expStr, directory, outFolder, figTypes, ...
                    figWidth, figHeight, xLimits, yLimits)

% Create a figure name
figPathBaseM2hOrig = fullfile(outFolder, [expStr, '_m2h_orig']);
figPathBaseM2h = fullfile(outFolder, [expStr, '_m2h']);

% Create the figure
figM2h = set_figure_properties('AlwaysNew', true);

% Plot traces
m3ha_plot_simulated_traces('Directory', directory, 'ExpStr', expStr, ...
                'PlotType', 'm2h', 'FigHandle', figM2h, ...
                'FigTitle', 'suppress', ...
                'XLimits', xLimits, 'YLimits', yLimits);

% Save the figure
save_all_figtypes(figM2h, figPathBaseM2hOrig, 'png');

% Plot a scale bar
hold on
plot_scale_bar('x', 'XBarUnits', 'ms', 'XBarLength', 400, ...
                'XPosNormalized', 0.1, 'YPosNormalized', 0.8);

% Update figure for CorelDraw
update_figure_for_corel(figM2h, 'Units', 'centimeters', ...
                'Width', figWidth, 'Height', figHeight, ...
                'AlignSubplots', true, ...
                'RemoveXTicks', true, 'RemoveXRulers', true);

% Save the figure
save_all_figtypes(figM2h, figPathBaseM2h, figTypes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
