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
plotExamplesFlag = true;
plotBoxPlotsFlag = true;
plotBarPlotsFlag = true;

% Plot settings
exampleGincr = 200;
exampleHeight = 4;
exampleWidth = 3;
exampleLineWidth = 0.5;
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

%% Compute statistics
% TODO

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
    figs = cellfun(@(x) m3ha_plot_example_traces(x, figure02Dir, ...
                                            swpInfo, exampleGincr, ...
                                            exampleXlimits, exampleYlimits, ...
                                            pharmLabels, exampleLineWidth, ...
                                            exampleHeight, exampleWidth), ...
                    exampleCellNames);

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

function fig = m3ha_plot_example_traces (cellName, figure02Dir, ...
                                            swpInfo, gIncrOfInterest, ...
                                            xLimits, ylimits, ...
                                            yLabels, lineWidth, ...
                                            figHeight, figWidth)
%% Plots example traces
% TODO: Pull out as its own function with only cellName as the required argument

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
data = m3ha_import_raw_traces(exampleBases, 'Directory', figure02Dir);

% Extract time and voltage vectors
[tVecs, vVecs] = extract_columns(data, 1:2);

% Create figure
fig = set_figure_properties('AlwaysNew', true);

% Plot the traces
plot_traces(tVecs, vVecs, ...
            'PlotMode', 'parallel', 'LineWidth', lineWidth, ...
            'LinkAxesOption', 'xy', 'FigTitle', 'suppress', ...
            'XLimits', xLimits, 'Ylimits', ylimits, ...
            'XLabel', 'suppress', 'YLabel', yLabels, ...
            'LegendLocation', 'suppress', ...
            'FigHandle', fig);

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'inches', ...
                        'Height', figHeight, 'Width', figWidth);

% Create figure base
figBase = sprintf('%s_gincr_%g_examples', cellName, gIncrOfInterest);
figPathBase = fullfile(figure02Dir, figBase);

% Save the figure
save_all_figtypes(fig, figPathBase, {'png', 'epsc2'});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%