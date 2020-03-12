% m3ha_plot_figure04.m
%% Plots Figure 04 for the GAT Blocker paper
%
% Requires:
%       cd/archive_dependent_scripts.m
%       cd/argfun.m
%       cd/create_label_from_sequence.m
%       cd/create_labels_from_numbers.m
%       cd/find_matching_files.m
%       cd/force_string_end.m
%       cd/m3ha_compute_statistics.m
%       cd/m3ha_extract_cell_name.m
%       cd/m3ha_load_sweep_info.m
%       cd/m3ha_plot_bar3.m
%       cd/m3ha_plot_violin.m
%       cd/m3ha_select_sweeps.m

% File History:
% 2020-01-02 Modified from m3ha_plot_figure02.m
% 2020-03-10 Reordered measuresOfInterest
% 2020-03-10 Updated pharm labels

%% Hard-coded parameters
% Flags
plotViolinPlotsFlag = true;
plotBarPlotsFlag = false; %true;
archiveScriptsFlag = true;

% Directories
parentDirectory = fullfile('/media', 'adamX', 'm3ha');
figure02Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure02');
figure04Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure04');
fitDirName = 'optimizer4gabab';

% rankDirName = '20191229_ranked_singleneuronfitting0-91';
% rankNumsToUse = [1, 2, 5, 7, 8, 9, 10, 13, 17, 34];

% rankDirName = '20200103_ranked_singleneuronfitting0-94';
% rankNumsToUse = 1:11;

% rankDirName = '20200131_ranked_singleneuronfitting0-102';
% rankNumsToUse = 1:11;

% rankDirName = '20200203_ranked_manual_singleneuronfitting0-102';
% rankNumsToUse = [1, 2, 5:10, 12:25, 29, 33];

% rankDirName = '20200203_ranked_manual_singleneuronfitting0-102';
% rankNumsToUse = [1, 2, 4:10, 12:25, 29, 33];

rankDirName = '20200207_ranked_manual_singleneuronfitting0-102';
rankNumsToUse = 1:23;

% Files
datalogPath = fullfile(figure02Dir, 'dclampdatalog_take4.csv');

% Analysis settings
% Note: must be consistent with m3ha_compute_statistics.m
measuresOfInterest = {'ltsProbability'; 'ltsOnsetTime'; ...
                        'spikesPerLts'; 'ltsAmplitude'; 'ltsMaxSlope'; ...
                        'burstProbability'; 'burstOnsetTime'; 'spikesPerBurst'; ...
                        'ltsConcavity'; 'ltsProminence'; 'ltsWidth'; ...
                        'spikeMaxAmp'; 'spikeMinAmp'; ...
                        'spikeFrequency'; 'spikeAdaptation'; ...
                        'ltsTimeJitter'; 'burstTimeJitter'};
dataMode = 1;           % data mode:
                        %   0 - all data
                        %   1 - all of g incr = 100%, 200%, 400% 
                        %   2 - same g incr but exclude 
                        %       cell-pharm-g_incr sets 
                        %       containing problematic sweeps

% Plot settings
pharmAll = [1; 2; 3; 4];
pharmLabelsLong = {'{\it d}Control', '{\it d}GAT1-Block', ...
                    '{\it d}GAT3-Block', '{\it d}Dual-Block'};
pharmLabelsShort = {'{\it d}Con', '{\it d}GAT1', ...
                    '{\it d}GAT3', '{\it d}Dual'};
if dataMode == 0
    gIncrAll = [25; 50; 100; 200; 400; 800];
    gIncrLabels = {'25%', '50%', '100%', '200%', '400%', '800%'};
elseif dataMode == 1 || dataMode == 2
    gIncrAll = [100; 200; 400];
    gIncrLabels = {'100%', '200%', '400%'};
end
conditionLabel2D = 'pharm_1-4_gincr_200_rec';
pCond2D = num2cell(pharmAll);
gCond2D = 200;
conditionLabel3D = 'pharm_1-4_gincr_all_rec';
pCond3D = num2cell(pharmAll);
gCond3D = num2cell(gIncrAll);

% Note: Use the default in m3ha_plot_violin.m and m3ha_plot_bar3.m
% violinFigHeight = 3;            % in centimeters
% violinFigWidth = 3.4;             % in centimeters
% violinRelativeBandWidth = 0.1;    % bandwidth relative to data range
% medianColor = rgb('GreenYellow'); % color of median circle
% medianSize = 6;                   % size of median circle in points
% bar3FigHeight = 4.3;              % in centimeters
% bar3FigWidth = 4.3;               % in centimeters

figTypes = {'png', 'epsc'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Locate the fit directory
fitDirectory = fullfile(parentDirectory, fitDirName);

% Locate the ranked directory
rankDirectory = fullfile(fitDirectory, rankDirName);

% Read from datalogPath
swpInfo = m3ha_load_sweep_info('Directory', figure02Dir);

% Create rank number prefixes
rankPrefixes = create_labels_from_numbers(rankNumsToUse, ...
                                    'Prefix', 'rank_', 'Suffix', '_');

% Create a cell choice string
cellsStr = [rankDirName, '_rank', create_label_from_sequence(rankNumsToUse)];

% Update condition labels
[conditionLabel2D, conditionLabel3D] = ...
    argfun(@(x) force_string_end(x, ['_', cellsStr]), ...
            conditionLabel2D, conditionLabel3D);

% Find png files matching the rank prefixes
[~, pngPaths] = find_matching_files(rankPrefixes, 'PartType', 'prefix', ...
                        'Directory', rankDirectory, 'Extension', 'png', ...
                        'ExtractDistinct', false);

% Extract the cell names
cellNames = m3ha_extract_cell_name(pngPaths, 'FromBaseName', true);

% Select sweeps based on data mode
swpInfo = m3ha_select_sweeps('SwpInfo', swpInfo, 'Verbose', false, ...
                                'DataMode', dataMode, 'CellNames', cellNames);

%% Plot 2D violin plots
if plotViolinPlotsFlag
    % Construct stats table path
    stats2dPath = fullfile(figure04Dir, strcat(conditionLabel2D, '_stats.mat'));

    % Compute statistics if not done already
    if ~isfile(stats2dPath)
        % Compute statistics for all features
        disp('Computing statistics for violin plots ...');
        statsTable = m3ha_compute_statistics('SwpInfo', swpInfo, ...
                                                'PharmConditions', pCond2D, ...
                                                'GIncrConditions', gCond2D, ...
                                                'DataMode', dataMode);

        % Generate labels
        conditionLabel = conditionLabel2D;
        pharmLabels = pharmLabelsShort;

        % Save stats table
        save(stats2dPath, 'statsTable', 'pharmLabels', ...
                            'conditionLabel', '-v7.3');
    end

    % Plot violin plots
    m3ha_plot_violin(stats2dPath, 'RowsToPlot', measuresOfInterest, ...
                    'OutFolder', figure04Dir);
end

%% Plot 3D bar plots
if plotBarPlotsFlag
    % Construct stats table path
    stats3dPath = fullfile(figure04Dir, strcat(conditionLabel3D, '_stats.mat'));

    % Compute statistics if not done already
    if ~isfile(stats3dPath)
        % Compute statistics for all features
        disp('Computing statistics for 3D bar plots ...');
        statsTable = m3ha_compute_statistics('SwpInfo', swpInfo, ...
                                                'PharmConditions', pCond3D, ...
                                                'GIncrConditions', gCond3D, ...
                                                'DataMode', dataMode);

        % Generate labels
        conditionLabel = conditionLabel3D;
        pharmLabels = pharmLabelsLong;

        % Save stats table
        save(stats3dPath, 'statsTable', 'pharmLabels', ...
                        'gIncrLabels', 'conditionLabel', '-v7.3');
    end

    % Plot bar plots
    m3ha_plot_bar3(stats3dPath, 'RowsToPlot', measuresOfInterest, ...
                    'OutFolder', figure04Dir);
end

% Archive all scripts for this run
if archiveScriptsFlag
    archive_dependent_scripts(mfilename, 'OutFolder', figure04Dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
