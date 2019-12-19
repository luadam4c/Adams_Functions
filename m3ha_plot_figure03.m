% m3ha_plot_figure03.m
%% Plots Figure 03 for the GAT Blocker paper
%
% Requires:
% TODO
%       cd/m3ha_load_sweep_info.m

% File History:
% 2019-12-18 Created by Adam Lu

%% Hard-coded parameters
% Directories
parentDirectory = fullfile('/media', 'adamX', 'm3ha');
figure02Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure02');
figure03Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure03');
matFilesDir = fullfile(parentDirectory, 'data_dclamp', 'take4', 'matfiles');

% Files
datalogPath = fullfile(figure02Dir, 'dclampdatalog_take4.csv');

% Flags
estimatePassiveParams = true;

% Analysis settings
exampleCellNames = {'D101310', 'C101210'};

% Data mode
dataMode = 2;

% Plot settings
figTypes = {'png', 'epsc2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load sweep info
% Read from datalogPath
swpInfo = m3ha_load_sweep_info('Directory', figure02Dir);

%% Perform curve fitting to estimate passive parameters
if estimatePassiveParams
    cellfun(@(x) estimate_passive_params_for_one_cell(x, ...
                    swpInfo, dataMode, matFilesDir), ...
            exampleCellNames);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function estimate_passive_params_for_one_cell(cellName, swpInfo, ...
                                        dataMode, matFilesDir)

% Select sweeps for this cell and this data mode
[swpInfo, fileBasesToUse] = ...
    m3ha_select_sweeps('SwpInfo', swpInfo, 'DataMode', dataMode, ...
                    'CellNames', cellName);

% Construct full paths to data files
dataFilePaths = fullfile(matFilesDir, strcat(fileBasesToUse, '.mat'));

% Load all sweeps for this cell and this data mode


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
