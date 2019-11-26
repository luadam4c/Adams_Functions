% m3ha_plot_figure02.m
%% Plots Figure 02 for the GAT Blocker paper
%
% Requires:
%       cd/all_files.m
%       cd/copy_into.m

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

% Sizes
exampleHeight = 1;
exampleWidth = 3;

% Flags
saveCellInfo = false;
copyExampleFiles = false; %true;
plotExamplesFlag = true;
plotBoxPlotsFlag = true;
plotBarPlotsFlag = true;

% Plot settings
stimStartWindow = 1000;
responseWindow = [800, 2000];
gincrExamples = 200;
pharmLabels = {'Control', 'GAT1', 'GAT3', 'Dual'};
figTypes = {'png', 'epsc2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load sweep info
% Read from datalogPath
swpInfo = readtable(datalogPath, 'HeaderLines', 1);

% Extract file bases as row names
fileNames = swpInfo.fnrow;
fileBases = extractBefore(fileNames, '.');
swpInfo.Properties.RowNames = fileBases;

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
    % Look for all .mat files
    [~, matPathsAll] = all_files('Directory', figure02Dir, ...
                            'KeyWord', exampleCellNames{1}, 'Extension', 'mat');

    % Extract file bases
    matBasesAll = extract_fileparts(matPathsAll, 'base');

    % Get all the conductance amplitudes
    gIncrAll = swpInfo(matBasesAll, 'grow');

    % Restrict to the conductance amplitude of interest
    exampleBases = matBasesAll(gIncrAll == gincrExamples);

    % Import traces
    data = m3ha_import_raw_traces(exampleBases, 'Directory', matFilesDir, ...
                                'StimStartWindow', stimStartWindow, ...
                                'ResponseWindow', responseWindow);

    % Extract time and voltage vectors
    [tVecs, vVecs] = extract_columns(data, 1:2);

    % Create figure
    fig = set_figure_properties;

    % Plot the traces
    plot_traces(tVecs, vVecs, ...
                'PlotMode', 'parallel', 'LineWidth', 2, ...
                'LinkAxesOption', 'xy', 'FigTitle', 'suppress', ...
                'XLabel', 'suppress', 'YLabel', pharmLabels, ...
                'LegendLocation', 'suppress', ...
                'FigHandle', fig);

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

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%