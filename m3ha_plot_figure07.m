% m3ha_plot_figure07.m
%% Plots Figure 07 for the GAT Blocker paper
%
% Requires:
%       cd/archive_dependent_scripts.m
%       cd/m3ha_network_compare_ipsc.m
%       cd/m3ha_network_plot_essential.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/update_figure_for_corel.m

% File History:
% 2020-01-30 Modified from m3ha_plot_figure05.m

%% Hard-coded parameters
% Flags
plotIpscComparison = true;
plot2CellExamples = true;
plot2CellPopulation = false; %true;

archiveScriptsFlag = false; %true;

% Directories
parentDirectory = fullfile('/media', 'adamX', 'm3ha');
figure07Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure07');
networkDirectory = fullfile(parentDirectory, 'network_model');
iterName = '20200129T1001_using_bestparams_20200126_singleneuronfitting101';

% Files

% Analysis settings
exampleCellNames = {'D101310'; 'M101210'};
gIncr = 200;            % Original dynamic clamp gIncr value

% Plot settings
ipscFigWidth = 8.5;
ipscFigHeight = 4;
exampleFigWidth = 8.5;
exampleFigHeight = 8.5;

figTypes = {'png', 'epsc2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Find the directory for this iteration
iterDir = fullfile(networkDirectory, iterName);

%% Find example files and directories
if plotIpscComparison || plot2CellExamples
    % Find example network directories
    [~, exampleDirs] = ...
        cellfun(@(x) all_subdirs('Directory', iterDir, 'Keyword', x), ...
                exampleCellNames, 'UniformOutput', false);
end

%% Plots figures for comparing dynamic clamp ipsc
if plotIpscComparison
    cellfun(@(x, y) plot_ipsc_comparison(x, iterName, gIncr, y, ...
                                        figure07Dir, figTypes, ...
                                        ipscFigWidth, ipscFigHeight), ...
            exampleCellNames, exampleDirs);
end

%% Plots example 2-cell networks
if plot2CellExamples
    for pharm = 1:4
        cellfun(@(x, y) plot_2cell_examples(x, iterName, gIncr, pharm, y, ...
                                        figure07Dir, figTypes, ...
                                        exampleFigWidth, exampleFigHeight), ...
                exampleCellNames, exampleDirs);
    end
end

%% Plots quantification over all 2-cell networks
if plot2CellPopulation
    % TODO
end

%% Archive all scripts for this run
if archiveScriptsFlag
    archive_dependent_scripts(mfilename, 'OutFolder', figure07Dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_ipsc_comparison(cellName, iterName, gIncr, ...
                                inFolder, outFolder, ...
                                figTypes, figWidth, figHeight)
% Plot an IPSC comparison plot

% Create a gIncr string
gIncrStr = ['gIncr', num2str(gIncr)];

% Create figure names
figPathBase = fullfile(outFolder, [cellName, '_', iterName, '_', gIncrStr, ...
                                    '_gabab_ipsc_comparison']);
figPathBaseOrig = [figPathBase, '_orig'];

% Create the figure
fig = set_figure_properties('AlwaysNew', true);

% Plot comparison
m3ha_network_compare_ipsc('SaveNewFlag', false, 'InFolder', inFolder, ...
                            'AmpScaleFactor', gIncr);

% Save original figure
save_all_figtypes(fig, figPathBaseOrig, 'png');

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Width', figWidth, 'Height', figHeight, ...
                        'AlignSubplots', true);

% Save the figure
save_all_figtypes(fig, figPathBase, figTypes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_2cell_examples(cellName, iterName, gIncr, pharm, ...
                                inFolder, outFolder, ...
                                figTypes, figWidth, figHeight)
% Plot 2-cell network examples

% Create a gIncr string
gIncrStr = ['gIncr', num2str(gIncr)];
pharmStr = ['pharm', num2str(pharm)];

% Create figure names
figPathBase = fullfile(outFolder, [cellName, '_', iterName, '_', gIncrStr, ...
                                    '_', pharmStr, '_example']);
figPathBaseOrig = [figPathBase, '_orig'];

% Create the figure
fig = set_figure_properties('AlwaysNew', true);

% Plot example
m3ha_network_plot_essential('SaveNewFlag', false, 'InFolder', inFolder, ...
                            'AmpScaleFactor', gIncr, 'PharmCondition', pharm);

% Save original figure
save_all_figtypes(fig, figPathBaseOrig, 'png');

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Width', figWidth, 'Height', figHeight, ...
                        'AlignSubplots', true);

% Save the figure
save_all_figtypes(fig, figPathBase, figTypes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
