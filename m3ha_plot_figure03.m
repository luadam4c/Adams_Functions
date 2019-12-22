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
%       cd/plot_ball_stick.m
%       cd/plot_scale_bar.m

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
estimatePassiveParams = false; %true;
plotCurveFit = true;

% Analysis settings
exampleCellNames = {'D101310', 'C101210'};

% Data mode
dataMode = 2;

% Output files
passiveLogSuffix = 'dclampPassiveLog';

% Plot settings
somaColor = rgb('DarkGreen');
dendriteColor = rgb('DarkOrange');
curveFigWidth = 8.5;
curveFigHeight = 4;
curveXLimits = [0, 60];
curveYLimits = [-1.2, 0];
geomFigWidth = 4;
geomFigHeight = 7;
geomXLimits = [-100, 100];
geomYLimits = [-100, 100];

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

%% Plot curve fitting results
if plotCurveFit
    % Find passive log paths
    [~, passiveLogPaths] = ...
        find_matching_files(exampleCellNames, 'Directory', figure03Dir, ...
                            'Suffix', passiveLogSuffix, 'Extension', 'mat');


    cellfun(@(x, y) plot_curve_fit_results(x, y, figure03Dir, figTypes, ...
                            somaColor, dendriteColor, ...
                            curveFigWidth, curveFigHeight, ...
                            curveXLimits, curveYLimits, ...
                            geomFigWidth, geomFigHeight, ...
                            geomXLimits, geomYLimits), ...
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

function plot_curve_fit_results(cellName, passiveLogPath, ...
                outFolder, figTypes, somaColor, dendriteColor, ...
                curveFigWidth, curveFigHeight, curveXLimits, curveYLimits, ...
                geomFigWidth, geomFigHeight, geomXLimits, geomYLimits)

% Decide on the figure name
figPathBaseCurveFit = fullfile(outFolder, [cellName, '_curve_fit']);
figPathBaseGeom = fullfile(outFolder, [cellName, '_geometry']);

% Read the passive log file
m = matfile(passiveLogPath);

% Load variables needed for plotting
tVecFitted = m.tVecFitted;
vVecFitted = m.vVecFitted;
fitObject = m.fitObject;
fitResults = m.fitResults;
goodnessOfFit = m.goodnessOfFit;
passiveParams = m.passiveParams;

%% Plot curve fit
% Create the figure
figCurveFit = set_figure_properties('AlwaysNew', true);

% Plot curve fit
plot_cfit_pulse_response(tVecFitted, vVecFitted, ...
                        'FitObject', fitObject, ...
                        'FitResults', fitResults, ...
                        'GoodnessOfFit', goodnessOfFit, ...
                        'PassiveParams', passiveParams, ...
                        'LegendLocation', 'suppress', ...
                        'Component1Color', somaColor, ...
                        'Component2Color', dendriteColor);

% Set axis limits
xlim(curveXLimits)
ylim(curveYLimits)

% Set title
title(['Curve Fit for ', cellName]);

% Remove texts
update_figure_for_corel(figCurveFit, 'RemoveTexts', true);

% Plot a scale bar
plot_scale_bar('xy', 'XBarUnits', 'ms', 'YBarUnits', 'mV', ...
                'XBarLength', 5, 'YBarLength', 0.1, ...
                'XPosNormalized', 0.6, 'YPosNormalized', 0.2);

% Update figure for CorelDraw
update_figure_for_corel(figCurveFit, 'Units', 'centimeters', ...
                        'Width', curveFigWidth, 'Height', curveFigHeight, ...
                        'RemoveRulers', true, 'RemoveLabels', true);

% Save the figure
save_all_figtypes(figCurveFit, figPathBaseCurveFit, figTypes);

%% Plot ball-and-stick model
% Create the figure
[figGeom, axGeom] = create_subplots(2, 1, 'AlwaysNew', true);

% Top subplot
subplot(axGeom(1));

% Plot ball-and-stick model
plot_ball_stick('GeomParams', passiveParams, ...
                'BallEdgeColor', somaColor, 'StickEdgeColor', dendriteColor, ...
                'FaceColor', 'none', 'LineWidth', 1);

% Set axis limits
xlim(geomXLimits)
ylim(geomYLimits)

% Set title
title(['Ball-stick model for ', cellName]);

% Plot a scale bar
plot_scale_bar('xy', 'XBarUnits', 'um', 'YBarUnits', 'um', ...
                'XBarLength', 30, 'YBarLength', 30, ...
                'XPosNormalized', 0.6, 'YPosNormalized', 0.2);

% Bottom subplot
subplot(axGeom(2));

% Plot final geometry
plot_ball_stick('GeomParams', passiveParams, ...
                'BallCurvature', [0, 0], 'NStickSegments', 2, ...
                'BallEdgeColor', somaColor, 'StickEdgeColor', dendriteColor, ...
                'FaceColor', 'none', 'LineWidth', 1);

% Set title
title(['Final geometry for ', cellName]);

% Link axes
linkaxes(axGeom, 'xy');

% Update figure for CorelDraw
update_figure_for_corel(figGeom, 'Units', 'centimeters', ...
                'Width', geomFigWidth, 'Height', geomFigHeight, ...
                'RemoveRulers', true);

% Save the figure
save_all_figtypes(figGeom, figPathBaseGeom, figTypes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
