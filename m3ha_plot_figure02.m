% m3ha_plot_figure02.m
%% Plots Figure 02 for the GAT Blocker paper
%
% Requires:
%       cd/all_files.m
%       cd/copy_into.m
%       cd/extract_columns.m
%       cd/extract_fileparts.m
%       cd/extract_substrings.m
%       cd/m3ha_create_cell_info_table.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_load_sweep_info.m
%       cd/plot_scale_bar.m
%       cd/plot_traces.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/update_figure_for_corel.m

% File History:
% 2019-11-25 Created by Adam Lu

%% Hard-coded parameters
% Directories
parentDirectory = fullfile('/media', 'adamX', 'm3ha');
figure02Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure02');
matFilesDir = fullfile(parentDirectory, 'data_dclamp', 'take4', 'matfiles');

% Files
datalogPath = fullfile(figure02Dir, 'dclampdatalog_take4.csv');
cellNamePattern = '[A-Z][0-9]{6}';

% Flags
saveCellInfo = false;
copyExampleFiles = false; %true;
plotExamplesFlag = false; %true;
plotBoxPlotsFlag = false; %true;
plotBarPlotsFlag = true;

% Analysis settings
% Note: must be consistent with m3ha_compute_ltsburst_statistics.m
measuresOfInterest = {'ltsOnsetTime'; 'ltsProbability'; 'spikesPerLts'; ...
                    'burstOnsetTime'; 'burstProbability'; 'spikesPerBurst'};

% Plot settings
exampleGincr = 200;
exampleHeight = 11;
exampleWidth = 8.5;             % in centimeters
exampleLineWidth = 0.5;         % in centimeters
exampleXlimits = [800, 2200];
exampleYlimits = [-100, 20];
pharmLabels = {'Control', 'GAT1 Block', 'GAT3 Block', 'Dual Block'};
figTypes = {'png', 'epsc2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load sweep info
% Read from datalogPath
swpInfo = m3ha_load_sweep_info('Directory', figure02Dir);

% Look for all import logs
[~, exampleLogPaths] = all_files('Directory', figure02Dir, ...
                        'Suffix', 'imported_files', 'Extension', 'txt');

% Extract example cell names
exampleCellNames = ...
    extract_substrings(exampleLogPaths, 'Regexp', cellNamePattern);

%% Create cell info table
if saveCellInfo
    cellInfo = m3ha_create_cell_info_table('SwpInfo', swpInfo);
end

%% Copy example traces
if copyExampleFiles
    % Copy files from the log files
    cellfun(@(x) copy_files_from_log(x, matFilesDir, figure02Dir), ...
            exampleLogPaths);
end

%% Plot example traces
if plotExamplesFlag
    % Plot example traces for each cell
    handles = cellfun(@(x) m3ha_plot_example_traces(x, figure02Dir, ...
                                            swpInfo, exampleGincr, ...
                                            exampleXlimits, exampleYlimits, ...
                                            pharmLabels, exampleLineWidth, ...
                                            exampleHeight, exampleWidth, ...
                                            figTypes), ...
                        exampleCellNames);

end

%% Compute statistics
if plotBoxPlotsFlag || plotBarPlotsFlag
    % Compute statistics for all features
    statsTable = ...
        m3ha_compute_ltsburst_statistics('PharmConditions', num2cell(1:4), ...
                                'GIncrCondition', num2cell([100; 200; 400]));

    % Restrict to measures of interest
    statsTable = statsTable(measuresOfInterest, :);
end

%% Plot 2D box plots
if plotBoxPlotsFlag
    
end

%% Plot 3D bar plots
if plotBarPlotsFlag
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function copy_files_from_log(logFile, dirFrom, dirTo)

% Read the log file
logTable = readtable(logFile, 'ReadVariableNames', false, 'Delimiter', ' ');

% Extract .mat file names
matFiles = logTable.Var1;

% Construct full paths to .mat files
matPaths = fullfile(dirFrom, matFiles);

% Copy files
copy_into(matPaths, dirTo);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = m3ha_plot_example_traces (cellName, figure02Dir, ...
                                            swpInfo, gIncrOfInterest, ...
                                            xLimits, yLimits, ...
                                            pharmLabels, lineWidth, ...
                                            figHeight, figWidth, figTypes)
%% Plots example traces
% TODO: Pull out as its own function with only cellName as the required argument

NS_PER_US = 1000;

%% Hard-coded parameters
conductanceLabel = 'Conductance (nS)';
conductanceYLimits = [0, 10];
timeBarLength = 100;
timeBarUnits = 'ms';

% Look for all .mat files
[~, matPathsAll] = all_files('Directory', figure02Dir, ...
                        'KeyWord', cellName, 'Extension', 'mat');

% Extract file bases
matBasesAll = extract_fileparts(matPathsAll, 'base');

% Get all the conductance amplitudes
gIncrAll = swpInfo{matBasesAll, 'grow'};

% Restrict to the conductance amplitude of interest
exampleBases = matBasesAll(gIncrAll == gIncrOfInterest);

% Get all the pharm conditions
pharmAll = swpInfo{exampleBases, 'prow'};

% Sort according to pharm condition
[~, origInd] = sort(pharmAll);
exampleBases = exampleBases(origInd);

% Import traces
data = m3ha_import_raw_traces(exampleBases, 'SweepInfoAll', swpInfo, ...
                                'OutFolder', figure02Dir);

% Extract time, voltage and conductance vectors
[tVecs, vVecs, gVecs] = extract_columns(data, [1:2, 4]);

% Create figure for voltage traces
figV = set_figure_properties('AlwaysNew', true);

% Plot the voltage traces
tracesV = plot_traces(tVecs, vVecs, ...
            'PlotMode', 'parallel', 'LineWidth', lineWidth, ...
            'LinkAxesOption', 'xy', 'FigTitle', 'suppress', ...
            'XLimits', xLimits, 'Ylimits', yLimits, ...
            'XLabel', 'suppress', 'YLabel', pharmLabels, ...
            'LegendLocation', 'suppress', ...
            'FigHandle', figV);

% Get voltage subplots
axV = tracesV.subPlots;

% Remove x axis
set(axV, 'XTick', []);

% Plot a time bar
plot_scale_bar('x', 'BarLength', timeBarLength, 'BarUnits', timeBarUnits);

% Update figure for CorelDraw
update_figure_for_corel(figV, 'Units', 'centimeters', ...
                        'Height', figHeight, 'Width', figWidth);

% Create figure base
figBase = sprintf('%s_gincr_%g_examples', cellName, gIncrOfInterest);
figPathBase = fullfile(figure02Dir, figBase);

% Save the figure
save_all_figtypes(figV, figPathBase, figTypes);

% Create figure for conductance traces
%   Note: prevent the legend from updating
figG = set_figure_properties('AlwaysNew', true, 'defaultLegendAutoUpdate', 'off');

% Convert conductance to nS
% TODO: use convert_units.m
gVecs = cellfun(@(x) x .* NS_PER_US, gVecs, 'UniformOutput', false);

% Plot the conductance traces
tracesG = plot_traces(tVecs, gVecs, ...
            'PlotMode', 'overlapped', 'LineWidth', lineWidth, ...
            'FigTitle', 'suppress', ...
            'XLimits', xLimits, 'YLimits', conductanceYLimits, ...
            'XLabel', 'suppress', 'YLabel', conductanceLabel, ...
            'TraceLabels', pharmLabels, ...
            'LegendLocation', 'northeast', ...
            'FigHandle', figG);

% Get conductance trace axes
axG = tracesG.subPlots;

% Remove x axis
set(axG, 'XTick', []);

% Plot a time bar
plot_scale_bar('x', 'BarLength', timeBarLength, 'BarUnits', timeBarUnits);

% TODO: match_axes_size.m 
%   match_axes_size(axG, axV(1), 'width')

% Get the axes width for the voltage plot
voltageSubPlotPosition = get(axV(1), 'Position');

% Update axes width to be consistent with the voltage plot
set_axes_properties('AxesHandle', axG, 'Width', voltageSubPlotPosition(3));

% Update figure for CorelDraw
update_figure_for_corel(figG, 'Units', 'centimeters', ...
                        'Height', figHeight / 4, 'Width', figWidth);

% Create figure base
figBase = sprintf('%s_gincr_%g_conductance', cellName, gIncrOfInterest);
figPathBase = fullfile(figure02Dir, figBase);

% Save the figure
save_all_figtypes(figG, figPathBase, figTypes);

% Output figure handles
handles.figG = figG;
handles.figV = figV;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

barColor = 'k';
barLineWidth = 1;
xPosNormalized = 0.8;
yPosBarNormalized = 0.1;
yPosTextNormalized = 0.05;
xPosBar = xLimits(1) + xPosNormalized * diff(xLimits);
yPosBar = yLimits(1) + yPosBarNormalized * diff(yLimits);
xPosText = xLimits(1) + xPosNormalized * diff(xLimits);
yPosText = yLimits(1) + yPosTextNormalized * diff(yLimits);
xLimitsBar = xPosBar + [0, barLengthMs];
plot_horizontal_line(yPosBar, 'XLimits', xLimitsBar, ...
                    'Color', barColor, 'LineWidth', barLineWidth, ...
                    'AxesHandle', subPlots(end));
text(xPosText, yPosText, barText);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%