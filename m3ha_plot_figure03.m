% m3ha_plot_figure03.m
%% Plots Figure 03 for the GAT Blocker paper
%
% Requires:
%       cd/find_matching_files.m
%       cd/find_passive_params.m
%       cd/m3ha_correct_unbalanced_bridge.m
%       cd/m3ha_load_sweep_info.m
%       cd/m3ha_parse_mat.m
%       cd/m3ha_select_sweeps.m
%       cd/m3ha_specs_for_datamode.m

% File History:
% 2019-12-18 Created by Adam Lu

%% Hard-coded parameters
% Directories
parentDirectory = fullfile('/media', 'adamX', 'm3ha');
figure02Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure02');
figure03Dir = fullfile(parentDirectory, 'manuscript', 'figures', 'Figure03');
matFilesDir = fullfile(parentDirectory, 'data_dclamp', 'take4', 'matfiles');

% Files
sweepInfoFile = 'dclampdatalog_take4.csv';
initialSlopesFile = 'initial_slopes_nSamplesForPlot_2_threeStdMainComponent.mat';
datalogPath = fullfile(figure02Dir, sweepInfoFile);
initialSlopesPath = fullfile(figure03Dir, initialSlopesFile);

% Flags
estimatePassiveParams = true;
plotCurveFit = true;

% Analysis settings
exampleCellNames = {'D101310', 'C101210'};

% Data mode
dataMode = 2;

% Output files
passiveLogSuffix = 'dclampPassiveLog';

% Plot settings
figTypes = {'png', 'epsc2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load sweep info
% Read from datalogPath
swpInfo = m3ha_load_sweep_info('Directory', figure02Dir);

%% Perform curve fitting to estimate passive parameters
if estimatePassiveParams
    cellfun(@(x) estimate_passive_params_for_one_cell(x, ...
                    swpInfo, dataMode, matFilesDir, ...
                    initialSlopesPath, figure03Dir, passiveLogSuffix), ...
            exampleCellNames);
end

%% Plot curve fits
if plotCurveFit
    % Find passive log paths
    [~, passiveLogPaths] = ...
        find_matching_files(exampleCellNames, 'Directory', figure03Dir, ...
                            'Suffix', passiveLogSuffix, 'Extension', 'mat');


    cellfun(@(x, y) plot_curve_fit(x, y, figure03Dir, figTypes), ...
            exampleCellNames, passiveLogPaths);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function estimate_passive_params_for_one_cell(cellName, swpInfo, ...
                                            dataMode, matFilesDir, ...
                                            initialSlopesPath, outFolder, ...
                                            outMatSuffix)

% Hard-coded parameters
%   Note: must be consistent with m3ha_estimate_passive_params.m
cpWin = [95, 115];          % window in which the current pulse would lie (ms) 
                            %       (Supposed to be 100-110 ms but 
                            %           there will be offset)
cprWin = [95, 260];         % window in which the current pulse response 
                            %   would lie (ms)
plotFlag = true;

% Decide on an output file base
outFileBase = cellName;
outMatPath = fullfile(outFolder, [outFileBase, '_', outMatSuffix, '.mat']);

% Set suffix and title modification according to dataMode
[suffix, titleMod] = m3ha_specs_for_datamode(dataMode);
fprintf('Using fit mode == %d ... \n', dataMode);

% Select sweeps for this cell and this data mode
[swpInfo, dataFileBases] = ...
    m3ha_select_sweeps('SwpInfo', swpInfo, 'DataMode', dataMode, ...
                    'CellNames', cellName);

% Extract parameters
actIhold = swpInfo{dataFileBases, 'actIhold'};

% Construct full paths to data files
dataFilePaths = fullfile(matFilesDir, strcat(dataFileBases, '.mat'));

% Load vectors from data matfiles
%   restricted to the current pulse response window
[~, parsedData] = m3ha_parse_mat(dataFilePaths, 'LoadWindow', cprWin);
tvec0 = parsedData.tvec0;
ivec0s = parsedData.ivec0s;
vvec0s = parsedData.vvec0s;
ivec1s = parsedData.ivec1s;

% Fix current pulse response traces that may have 
%   out-of-balance bridges
vvec0s = m3ha_correct_unbalanced_bridge(dataFileBases, vvec0s, ...
                                        ivec0s, initialSlopesPath);

% Analyze passive parameters such as input resistance (MOhm)
fprintf('ANALYZING passive parameters for %s ...\n', outFileBase);
[passiveParams, fitResults, fitObject, ...
    goodnessOfFit, algorithmInfo, decision, tVecFitted, vVecFitted] = ...
    find_passive_params(tvec0, ivec0s, vvec0s, ...
                         'HoldCurrent', actIhold, ...
                         'PulseWindow', cpWin, ...
                         'PulseResponseWindow', cprWin, ...
                         'PlotFlag', plotFlag, ...
                         'OutFolder', outFolder, ...
                         'FileBase', outFileBase, ...
                         'Ivec1s', ivec1s, ...
                         'Suffix', suffix, 'TitleMod', titleMod);

save(outMatPath, 'cellName', 'dataMode', 'matFilesDir', 'initialSlopesPath', ...
        'outFolder', 'cpWin', 'cprWin', 'outMatSuffix', 'plotFlag', ...
        'outFileBase', 'outMatPath', 'suffix', 'titleMod', ...
        'swpInfo', 'dataFileBases', 'actIhold', 'dataFilePaths', ...
        'tvec0', 'ivec0s', 'vvec0s', 'ivec1s', ...
        'passiveParams', 'fitResults', 'fitObject', ...
        'goodnessOfFit', 'algorithmInfo', 'decision', ...
        'tVecFitted', 'vVecFitted', ...
        '-v7.3');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_curve_fit(cellName, passiveLogPath, outFolder, figTypes)

% Decide on the figure name
figPathBase = fullfile(outFolder, [cellName, '_curve_fit']);

% Read the passive log file
m = matfile(passiveLogPath);

% Load variables needed for plotting
tVecFitted = m.tVecFitted;
vVecFitted = m.vVecFitted;
fitObject = m.fitObject;
fitResults = m.fitResults;
goodnessOfFit = m.goodnessOfFit;
passiveParams = m.passiveParams;

% Create the figure
fig = set_figure_properties('AlwaysNew', true);

% Plot curve fit
plot_cfit_pulse_response(tVecFitted, vVecFitted, ...
                        'FitObject', fitObject, ...
                        'FitResults', fitResults, ...
                        'GoodnessOfFit', goodnessOfFit, ...
                        'PassiveParams', passiveParams, ...
                        'LegendLocation', 'suppress');

% Save the figure
save_all_figtypes(fig, figPathBase, figTypes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
