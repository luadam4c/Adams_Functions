% function [output1] = m3ha_simulate_population (reqarg1, varargin)
%% Generates simulated IPSC responses that can be compared with recorded data
% Usage: [output1] = m3ha_simulate_population (reqarg1, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
%
% Arguments:
%       reqarg1     - TODO: Description of reqarg1
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/all_files.m
%       cd/all_subdirs.m
%       cd/argfun.m
%       cd/combine_param_tables.m
%       cd/compute_rms_error.m
%       cd/compute_weighted_average.m
%       cd/copy_into.m
%       cd/create_time_stamp.m
%       cd/create_labels_from_numbers.m
%       cd/create_label_from_sequence.m
%       cd/extract_columns.m
%       cd/extract_elements.m
%       cd/extract_vars.m
%       cd/find_matching_files.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
%       cd/is_field.m
%       cd/load_neuron_outputs.m
%       cd/m3ha_compute_and_plot_violin.m
%       cd/m3ha_compute_statistics.m
%       cd/m3ha_extract_cell_name.m
%       cd/m3ha_extract_iteration_string.m
%       cd/m3ha_extract_sweep_name.m
%       cd/m3ha_find_decision_point.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_load_sweep_info.m
%       cd/m3ha_plot_bar3.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/m3ha_plot_violin.m
%       cd/plot_horizontal_line.m
%       cd/plot_violin.m
%       cd/print_cellstr.m
%       cd/regroup_cell_of_cells.m
%       cd/renamevars_custom.m
%       cd/restrict_values.m
%       cd/vertcat_spreadsheets.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/test_difference.m
%       cd/update_neuron_scripts.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-12-11 Created by Adam Lu
% 2019-12-26 Completed
% 2019-12-27 Added HH channels
% 2020-02-18 Now computes open probability discrepancies over fit window only
% 2020-02-23 Now saves open probability stats
% 2020-03-10 Reordered measuresOfInterest
% 2020-03-10 Updated pharm labels
% 2020-03-11 Now plots violin plots for all gIncr
% 2020-04-09 Now plots essential plots
% 2020-04-09 Now finds special cases
% 2020-04-10 Now computes summary cell info table
% 2020-04-11 Now plots correlations
% 2020-04-12 Added m2hMaxError and m2hLogMaxError for open probability
% 2020-04-13 Fixed bugs in plotEssential and findSpecialCases
% 2020-04-27 Added plotIndividualFlag
% 2019-04-28 Changed timeToStabilize from 2000 to 3000

%% Hard-coded parameters
% Flags
updateScriptsFlag = false; %true;
chooseBestNeuronsFlag = false; %true;
simulateFlag = false; %true;
combineFeatureTablesFlag = false; %true;
computeOpenProbabilityFlag = false; %true;
plotIndividualFlag = false; %true;
plotEssentialFlag = false; %true;
plotOpenProbabilityFlag = true;
findSpecialCasesFlag = false; %true;
computeCellInfoTableFlag = false; %true;
plotCorrelationsFlag = false; %true;
plotViolinPlotsFlag = false; %true;
plotBarPlotsFlag = false; %true;
archiveScriptsFlag = true;

% yLabelLogOpd = 'rms(m^2h - m_{inf}^2h_{inf})';
% yLabelLogOpd = 'max(abs(m^2h - m_{inf}^2h_{inf}))';
% yLabelLogOpd = 'max(m^2h / m_{inf}^2h_{inf})';
% yLabelLogOpd = 'max(log(m^2h / m_{inf}^2h_{inf}))';
% logOpdStr = 'openProbabilityDiscrepancy';
% logOpdStr = 'm2hMaxLogRatio';
opdThreshold = 1e-2;
logOpdThreshold = log10(opdThreshold);

opdMeasureStrs = {'m2hLogMaxError'; 'vSlopeAtDecision'; 'm2hDiffSlopeAtDecision'};
yLabelsOpd = {'max(log(m^2h - m_{inf}^2h_{inf}))'; ...
                'Voltage Slope at Decision Point'; ...
                'Discrepancy Slope at Decision Point'};
yLimitsOpd = {[-4.5, 0.5]; [0, 0.35]; [0, 0.08]};
opdThresholds = {logOpdThreshold; []; []};

% logOpdStr = 'm2hDiffSlopeAtDecision';
% yLabelLogOpd = 'Discrepancy Slope at Decision Point';
% yLimitsOpd = [];
% logOpdThreshold = [];

% Simulation parameters
useHH = true;           % whether to use Hudgin-Huxley Na+ and K+ channels
useCvode = false;       % whether to use variable integration time steps
buildMode = 'active';
simMode = 'active';
dataMode = 1; %0;           % data mode:
                        %   0 - all data
                        %   1 - all of g incr = 100%, 200%, 400% 
                        %   2 - same g incr but exclude 
                        %       cell-pharm-g_incr sets 
                        %       containing problematic sweeps
attemptNumber = 3;      %   1 - Use 4 traces @ 200% gIncr 
                        %           for this data mode
                        %   2 - Use all traces @ 200% gIncr 
                        %           for this data mode
                        %   3 - Use all traces for this data mode
                        %   4 - Use 1 trace for each pharm x gIncr 
                        %           for this data mode
                        %   5 - Use 4 traces @ 400% gIncr 
                        %       for this data mode
                        %   6 - Same as 4 but prioritize least vHold
                        %   7 - Same as 1 but prioritize least vHold
                        %   8 - Same as 5 but prioritize least vHold

% Directory names
parentDirectoryTemp = '/media/adamX/m3ha';
fitDirName = 'optimizer4gabab';
defaultOutFolderStr = 'population';

% File names
simStr = 'sim';
paramsSuffix = 'params';
ltsParamsSuffix = 'ltsParams';
simLtsParamsSuffix = strcat(simStr, '_', ltsParamsSuffix);
simSwpInfoSuffix = strcat(simStr, '_swpInfo');
simCellInfoSuffix = strcat(simStr, '_cellInfo');
simOutExtension = 'out';

TIME_COL_REC = 1;

% Column numbers for simulated data
%   Note: Must be consistent with singleneuron4compgabab.hoc
TIME_COL_SIM = 1;
VOLT_COL_SIM = 2;
IT_M_DEND2 = 47;
IT_MINF_DEND2 = 48;
IT_H_DEND2 = 49;
IT_HINF_DEND2 = 50;

% Note: The following must be consistent with m3ha_parse_dclamp_data.m
condVarStrs = {'cellidrow', 'prow', 'vrow', 'grow', 'swpnrow', ...
                'gabab_amp', 'gabab_Trise', 'gabab_TfallFast', ...
                'gabab_TfallSlow', 'gabab_w'};
pharmAll = [1; 2; 3; 4];          
pharmLabelsLong = {'{\it s}Control', '{\it s}GAT1-Block', ...
                    '{\it s}GAT3-Block', '{\it s}Dual-Block'};
pharmLabelsShort = {'{\it s}Con', '{\it s}GAT1', ...
                    '{\it s}GAT3', '{\it s}Dual'};
if dataMode == 0 || dataMode == 3
    gIncrAll = [25; 50; 100; 200; 400; 800];
    gIncrLabels = {'25%', '50%', '100%', '200%', '400%', '800%'};
elseif dataMode == 1 || dataMode == 2
    gIncrAll = [100; 200; 400];
    gIncrLabels = {'100%', '200%', '400%'};
end
conditionLabels2D = [create_labels_from_numbers(gIncrAll, ...
                    'Prefix', 'pharm_1-4_gincr_', 'Suffix', ['_', simStr]); ...
                    'pharm_1-4_gincr_pooled'];
pConds2D = repmat({num2cell(pharmAll)}, numel(gIncrAll) + 1, 1);
gConds2D = [num2cell(gIncrAll); {gIncrAll}];
conditionLabel3D = 'pharm_1-4_gincr_all_sim';
pCond3D = num2cell(pharmAll);
gCond3D = num2cell(gIncrAll);
stats3dSuffix = strcat(simStr, '_', conditionLabel3D, '_stats.mat');

% The following must be consistent with singleneuron4compgabab.hoc
timeToStabilize = 3000;         % padded time (ms) to make sure initial value 
                                %   of simulations are stabilized
ipscrWindow = timeToStabilize + [0, 2800];     % only simulate up to that time
fitWindowIpscr = timeToStabilize + [1000, 2800];  
                                % the time window (ms) where all 
                                %   recorded LTS would lie

% Plot settings
% Note: must be consistent with m3ha_compute_statistics.m
measuresOfInterest = {'ltsProbability'; 'ltsOnsetTime'; ...
                    'spikesPerLts'; 'ltsAmplitude'; 'ltsMaxSlope'; ...
                    'burstProbability'; 'burstOnsetTime'; 'spikesPerBurst'; ...
                    'ltsConcavity'; 'ltsProminence'; 'ltsWidth'; ...
                    'spikeMaxAmp'; 'spikeMinAmp'; ...
                    'spikeFrequency'; 'spikeAdaptation'; ...
                    'ltsTimeJitter'; 'burstTimeJitter'};
openProbFigWidth = 5;       % (cm)
openProbFigHeight = 3;      % (cm)

% For summary cell info table
measuresOfInterestOrig = {'ltsProbability'; 'ltsPeakTime'; ...
                    'spikesPerPeak'; 'ltsPeakValue'; 'maxSlopeValue'; ...
                    'burstProbability'; 'burstTime'; 'spikesPerBurst'; ...
                    'maxSpikeAmp'; 'minSpikeAmp'; ...
                    'spikeFrequency'; 'spikeAdaptation'};
measuresOfInterestNew = {'ltsProbability'; 'ltsLatency'; ...
                    'spikesPerLts'; 'ltsPeakValue'; 'ltsMaxSlope'; ...
                    'burstProbability'; 'burstTime'; 'spikesPerBurst'; ...
                    'spikeMaxAmp'; 'spikeMinAmp'; ...
                    'spikeFrequency'; 'spikeAdaptation'};
measuresIsLogNew = false(12, 1);
paramsOfInterest = {'diamSoma'; 'LDend'; 'diamDend'; 'gpas'; ...
                    'pcabarITSoma'; 'pcabarITDend1'; 'pcabarITDend2'; ...
                    'ghbarIhSoma'; 'ghbarIhDend1'; 'ghbarIhDend2'; ...
                    'gkbarIKirSoma'; 'gkbarIKirDend1'; 'gkbarIKirDend2'; ...
                    'gkbarIASoma'; 'gkbarIADend1'; 'gkbarIADend2'; ...
                    'gnabarINaPSoma'; 'gnabarINaPDend1'; 'gnabarINaPDend2'};
paramIslog = [false; false; false; true; ...
                true; true; true; ...
                true; true; true; ...
                true; true; true; ...
                true; true; true; ...
                true; true; true];

% Compare with m3ha_plot_figure05.m
overlappedXLimits = timeToStabilize + [800, 2800]; %[800, 2000];
essentialYLimits = {[-110, -40]; [0, 20]; []; ...
                            [-20, 5]; [1e-8, 1e0]};
essentialYTickLocs = {-90:20:-50; 0:5:20; []; ...
                            -15:5:0; [1e-7, 1e-1]};

individualXLimits = timeToStabilize + [800, 2800];
individualYLimits = [];
individualYTickLocs = [];

% TODO: Make optional argument
% outFolder = '20191227_population_rank1-10_useHH_true';
% outFolder = fullfile(parentDirectoryTemp, fitDirName, ...
%         '20191230_population_singleneuronfitting0-91_rank1-2,5,7-10,13,17,34');
% prefix = '20191227_population';
% rankNumsToUse = [];
% maxRankToSim = 11;
% rankNumsToUse = [1, 2, 5, 6, 8, 9, 10, 11, 23, 34];
% rankDirName = '20191227_ranked_singleneuronfitting0-90';
% rankDirName = '20191229_ranked_singleneuronfitting0-91';
% rankNumsToUse = [1, 2, 5, 7, 8, 9, 10, 13, 17, 34];
% rankDirName = '20200103_ranked_singleneuronfitting0-94';
% outFolder = fullfile(parentDirectoryTemp, fitDirName, ...
%                     '20200106_population_rank1-11_dataMode1_attemptNumber3');
% rankNumsToUse = 1:11;
% rankDirName = '20200103_ranked_singleneuronfitting0-94';
% outFolder = '20200203_population_rank1-2,5-10,12-25,29,33_dataMode1_attemptNumber3';
% rankDirName = '20200203_ranked_manual_singleneuronfitting0-102';
% rankNumsToUse = [1, 2, 5:10, 12:25, 29, 33];
% outFolder = '20200204_population_rank1-2,5-10,12-25,29,33_dataMode1_attemptNumber3_vtraub-65';
% rankDirName = '20200203_ranked_manual_singleneuronfitting0-102';
% rankNumsToUse = [1, 2, 5:10, 12:25, 29, 33];
% outFolder = '20200204_population_rank1-2,4-10,12-25,29,33_dataMode1_attemptNumber3_vtraub-65';
% rankDirName = '20200203_ranked_manual_singleneuronfitting0-102';
% rankNumsToUse = [1, 2, 4:10, 12:25, 29, 33];
% rankNumsOpenProbability = [];   % same as rankNumsToUse
% rankNumsOpenProbability = [6, 9];

% outFolder = fullfile(parentDirectoryTemp, fitDirName, ......
%                     '20200208_population_rank1-23_dataMode1_attemptNumber3');
% rankDirName = '20200207_ranked_manual_singleneuronfitting0-102';
% rankNumsToUse = 1:23;
% rankNumsOpenProbability = 1:23;

% outFolder = fullfile(parentDirectoryTemp, fitDirName, ......
%       '20200423_population_rank1-31_dataMode1_attemptNumber3_cvode_off');
% rankDirName = '20200207_ranked_manual_singleneuronfitting0-102';
% rankNumsToUse = 1:31;
% rankNumsOpenProbability = 1:31;

% outFolder = fullfile(parentDirectoryTemp, fitDirName, ......
%       '20200425_population_rank1-23_dataMode1_attemptNumber3_cvode_off');
% rankDirName = '20200207_ranked_manual_singleneuronfitting0-102';
% rankNumsToUse = 1:23;
% rankNumsOpenProbability = 1:23;

% outFolder = fullfile(parentDirectoryTemp, fitDirName, ......
%         '20200426_population_rank1-23_dataMode1_attemptNumber3_cvode_on_same_as_20200208');
% rankDirName = '20200207_ranked_manual_singleneuronfitting0-102';
% rankNumsToUse = 1:23;
% rankNumsOpenProbability = 1:23;

% outFolder = fullfile(parentDirectoryTemp, fitDirName, ......
%         '20200426_population_rank1-23_dataMode1_attemptNumber3_vtraub-52');
% rankDirName = '20200207_ranked_manual_singleneuronfitting0-102';
% rankNumsToUse = 1:23;
% rankNumsOpenProbability = 1:23;

% outFolder = fullfile(parentDirectoryTemp, fitDirName, ......
%         '20200427_population_rank1-23_dataMode1_attemptNumber3_noHH_cvode_on');
% rankDirName = '20200207_ranked_manual_singleneuronfitting0-102';
% rankNumsToUse = 1:23;
% rankNumsOpenProbability = 1:23;

% outFolder = fullfile(parentDirectoryTemp, fitDirName, ......
%         '20200427_population_rank1-23_dataMode1_attemptNumber3_noHH_cvode_off');
% rankDirName = '20200207_ranked_manual_singleneuronfitting0-102';
% rankNumsToUse = 1:23;
% rankNumsOpenProbability = 1:23;

% outFolder = fullfile(parentDirectoryTemp, fitDirName, ......
%         '20200428_population_rank1-23_dataMode1_attemptNumber3_noHH_cvode_off');
% rankDirName = '20200207_ranked_manual_singleneuronfitting0-102';
% rankNumsToUse = 1:23;
% rankNumsOpenProbability = 1:23;

outFolder = fullfile(parentDirectoryTemp, fitDirName, ......
        '20200429_population_rank1-31_dataMode1_attemptNumber3_noHH_cvode_on');
rankDirName = '20200207_ranked_manual_singleneuronfitting0-102';
rankNumsToUse = 1:31;
rankNumsOpenProbability = 1:31;

% outFolder = fullfile(parentDirectoryTemp, fitDirName, ......
%         '20200429_population_rank1-31_dataMode1_attemptNumber3_noHH_cvode_off');
% rankDirName = '20200207_ranked_manual_singleneuronfitting0-102';
% rankNumsToUse = 1:31;
% rankNumsOpenProbability = 1:31;

% outFolder = fullfile(parentDirectoryTemp, fitDirName, ......
%         '20200430_population_rank1-31_dataMode1_attemptNumber3_useHH_cvode_on');
% rankDirName = '20200207_ranked_manual_singleneuronfitting0-102';
% rankNumsToUse = 1:31;
% rankNumsOpenProbability = 1:31;

% outFolder = fullfile(parentDirectoryTemp, fitDirName, ......
%         '20200430_population_rank1-31_dataMode1_attemptNumber3_useHH_cvode_off');
% rankDirName = '20200207_ranked_manual_singleneuronfitting0-102';
% rankNumsToUse = 1:31;
% rankNumsOpenProbability = 1:31;

prefix = '';
figTypes = {'png', 'epsc'};

filtWidthMs = 30;
itm2hDiffLowerLimit = 1e-9;
itm2hDiffLeftBound = 1e-7;

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% Deal with arguments
% Check number of required arguments
if nargin < 1    % TODO: 1 might need to be changed
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'reqarg1');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, reqarg1, varargin{:});
param1 = iP.Results.param1;

% Keep unmatched arguments for the TODO() function
otherArguments = iP.Unmatched;
%}

%% Preparation
% Locate the home directory
% parentDirectory = m3ha_locate_homedir;
parentDirectory = parentDirectoryTemp;

% Locate the fit directory
fitDirectory = fullfile(parentDirectory, fitDirName);

% Locate the ranked directory
rankDirectory = fullfile(fitDirectory, rankDirName);

% Decide on the ranking numbers of cells to simulate
if isempty(rankNumsToUse)
    rankNumsToUse = 1:maxRankToSim;
end

% Decide on the ranking numbers of cells to compute open probability
if isempty(rankNumsOpenProbability)
    rankNumsOpenProbability = rankNumsToUse;
end

% Decide on output folder
if isempty(outFolder)
    % Create a rank string
    rankStr = ['rank', create_label_from_sequence(rankNumsToUse)];

    % Create a data mode string
    dataModeStr = ['dataMode', num2str(dataMode)];
    attemptNumberStr = ['attemptNumber', num2str(attemptNumber)];

    % Create output folder name
    outFolderName = strcat(create_time_stamp('FormatOut', 'yyyymmdd'), ...
                            '_', defaultOutFolderStr, '_', rankStr, ...
                            '_', dataModeStr, '_', attemptNumberStr);

    % Create full path to output folder
    outFolder = fullfile(fitDirectory, outFolderName);
end

% Create output folder for correlations
outFolderCorrMeas2Param = fullfile(outFolder, 'correlations-measures-to-params');
outFolderCorrParam2Param = fullfile(outFolder, 'correlations-params-to-params');

% Check if output folders exists
check_dir({outFolder, outFolderCorrMeas2Param, outFolderCorrParam2Param});

% Decide on output prefix
if isempty(prefix)
    % Extract output folder base name
    prefix = extract_fileparts(outFolder, 'dirbase');
end

% Construct full paths
simSwpInfoPath = fullfile(outFolder, [prefix, '_', simSwpInfoSuffix, '.csv']);
conditionLabels2D = strcat(prefix, '_', conditionLabels2D);
stats3dPath = fullfile(outFolder, [prefix, '_', stats3dSuffix, '.mat']);

% Create path to cell info table
simCellInfoPath = replace(simSwpInfoPath, simSwpInfoSuffix, ...
                            simCellInfoSuffix);

%% Make sure NEURON scripts are up to date in outFolder
if updateScriptsFlag
    update_neuron_scripts(fitDirectory, outFolder);
end

%% Choose the best cells and the best parameters for each cell
if chooseBestNeuronsFlag
    % Extract the corresponding cell names and iteration strings from 
    %   given rank numbers
    [cellNames, iterStrs] = ...
        extract_from_rank_numbers(rankNumsToUse, rankDirectory);

    % Find the parameter file directories
    [~, paramDirs] = cellfun(@(x) all_subdirs('Directory', rankDirectory, ...
                                        'RegExp', x, 'MaxNum', 1), ...
                            iterStrs, 'UniformOutput', false);

    % Find the parameter files for each cell
    [~, paramPaths] = cellfun(@(x, y) all_files('Directory', x, ...
                            'Keyword', y, 'Suffix', 'params', 'MaxNum', 1), ...
                            paramDirs, cellNames, 'UniformOutput', false);

    % Copy the parameter files into 
    copy_into(paramPaths, outFolder);
end

%% Find candidate NEURON parameter files
if simulateFlag || computeCellInfoTableFlag
    % Decide on candidate parameters files
    [~, paramPaths] = all_files('Directory', outFolder, 'Recursive', false, ...
                                'Suffix', paramsSuffix, 'Extension', 'csv');

    % Remove the ltsParams files
    paramPaths = paramPaths(~contains(paramPaths, 'ltsParams'));

    % Extract the cell names
    cellNames = m3ha_extract_cell_name(paramPaths);
end

%% Simulate
if simulateFlag
    % Change to the output directory
    cd(outFolder);

    % Find all simulated LTS stats spreadsheets
    [~, simLtsParamPaths] = ...
        cellfun(@(x) find_matching_files(x, 'Directory', outFolder, ...
                        'Suffix', simLtsParamsSuffix, 'Extension', 'csv', ...
                        'Recursive', false, 'ReturnEmpty', true), ...
                cellNames, 'UniformOutput', false);

    % Don't simulate again if LTS stats spreadsheets already exist
    alreadyDone = cellfun(@isfile, simLtsParamPaths);
    paramPaths = paramPaths(~alreadyDone);
    cellNames = cellNames(~alreadyDone);

    % Display message
    fprintf('All sweeps from the following cells will be simulated: \n');
    print_cellstr(cellNames, 'OmitBraces', true, 'Delimiter', '\n');

    % Run simulations for each parameter file
    cellfun(@(x) m3ha_neuron_run_and_analyze (x, 'DataMode', dataMode, ...
                    'IpscrWindow', ipscrWindow, ...
                    'FitWindowIpscr', fitWindowIpscr, ...
                    'BuildMode', buildMode, 'SimMode', simMode, ...
                    'UseHH', useHH, 'UseCvode', useCvode, ...
                    'AttemptNumber', attemptNumber, ...
                    'SaveSimOutFlag', true, 'SaveLtsInfoFlag', true), ...
            paramPaths);

    % Find all simulated LTS stats spreadsheets
    [~, simLtsParamPaths] = ...
        all_files('Directory', outFolder, 'Recursive', true, ...
                    'Suffix', simLtsParamsSuffix, 'Extension', 'csv');

    % Copy over the spreadsheets
    copy_into(simLtsParamPaths, outFolder);
end

%% Combine LTS & burst feature tables
if combineFeatureTablesFlag
    % Display message
    fprintf('Combining LTS & burst statistics ... \n');

    % Find all simulated LTS stats spreadsheets
    %   Note: Ignore stuff in backup folders
    [~, simLtsParamPaths] = ...
        all_files('Directory', outFolder, 'Suffix', simLtsParamsSuffix, ...
                    'Extension', 'csv', 'Recursive', false);

    % Combine the spreadsheets
    simSwpInfo = vertcat_spreadsheets(simLtsParamPaths);

    % Rename variables
    simSwpInfo = renamevars_custom(simSwpInfo, 'fileBase', 'simFileBase');

    % Extract the simulation file bases
    simFileBase = simSwpInfo.simFileBase;

    % Extract the original file bases
    fileBase = extractBefore(simFileBase, '_sim');

    % Make the original sweep names the row names
    simSwpInfo = addvars(simSwpInfo, fileBase, 'Before', 1);
    simSwpInfo.Properties.RowNames = fileBase;

    % Load original sweep info
    origSwpInfo = m3ha_load_sweep_info;

    % Extract the condition info
    condInfo = origSwpInfo(:, condVarStrs);

    % Join the condition info
    simSwpInfo = join(simSwpInfo, condInfo, 'Keys', 'Row');

    % Save the simulated sweep info table
    writetable(simSwpInfo, simSwpInfoPath, 'WriteRowNames', true);
end

%% Compute the rms error between open probability and its steady state
if computeOpenProbabilityFlag
    % Make sure the simulated sweep info table exists
    if ~isfile(simSwpInfoPath)
        error('Save a simulated sweep info table first!');
    end

    % Display message
    fprintf('Computing open probability discrepancies ... \n');

    % Read the simulated sweep info table
    simSwpInfo = readtable(simSwpInfoPath, 'ReadRowNames', true);

    % Extract the file bases
    fileBases = simSwpInfo.Properties.RowNames;

    % Locate the matching simulated output files
    [~, simOutPaths] = find_matching_files(fileBases, 'PartType', 'suffix', ...
                                        'Directory', outFolder, ...
                                        'Recursive', true, ...
                                        'Extension', simOutExtension);

    % Import one recorded trace for the time vector
    realData = m3ha_import_raw_traces(fileBases{1}, 'ImportMode', simMode, ...
                                'Verbose', false, 'OutFolder', outFolder);

    % Extract time vectors from recorded data
    tVecsRec = extract_columns(realData, TIME_COL_REC);

    % Load corresponding simulated data
    % If recorded data provided (tVecs not empty at this point),
    %   interpolate simulated data to match the time points of recorded data
    % Note: This is necessary because CVODE (variable time step method) 
    %       is applied in NEURON
    simData = load_neuron_outputs('FileNames', simOutPaths, ...
                                    'tVecs', tVecsRec, ...
                                    'TimeWindow', fitWindowIpscr);

    % Extract vectors from simulated data
    [tVecs, vVecs, m, minf, h, hinf] = ...
        extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                                    IT_M_DEND2, IT_MINF_DEND2, ...
                                    IT_H_DEND2, IT_HINF_DEND2]);

    % Compute dv/dt
    [dvdt, t1Vecs] = compute_derivative_trace(vVecs, tVecs);

    % Compute sampling interval
    siMs = compute_sampling_interval(tVecs);

    % Smooth dv/dt over filtWidthMs
    dvdtSmoothed = movingaveragefilter(dvdt, filtWidthMs, siMs);

    % Force as matrices
    [m, minf, h, hinf] = argfun(@force_matrix, m, minf, h, hinf);

    % Compute m2h and minf2hinf
    m2h = (m .^ 2) .* h;
    minf2hinf = (minf .^ 2) .* hinf;

    % Compute the rms error between m2h and minf2hinf
    m2hRmsError = compute_rms_error(m2h, minf2hinf, ...
                                    'ForceColumnOutput', true);

    % Compute the difference vector between m2h and minf2hinf
    m2hDiff = m2h - minf2hinf;

    % Compute the maximum absolute error between m2h and minf2hinf
    m2hMaxAbsError = force_column_vector(max(abs(m2hDiff), [], 1));

    % Compute the maximum error between m2h and minf2hinf
    m2hMaxError = force_column_vector(max(m2hDiff, [], 1));

    % Compute the maximum ratio between m2h and minf2hinf
    m2hRatioRaw = m2h ./ minf2hinf;
    m2hRatioRaw(isinf(m2hRatioRaw)) = NaN;
    m2hMaxRatio = force_column_vector(max(m2hRatioRaw, [], 1));

    % Compute the maximum log10 ratio between m2h and minf2hinf
    m2hLogRatioRaw = log10(m2h ./ minf2hinf);
    m2hLogRatioRaw(isinf(m2hLogRatioRaw)) = NaN;
    m2hMaxLogRatio = force_column_vector(max(m2hLogRatioRaw, [], 1));

    % Compute the log10 of everything
    m2hLogRmsError = log10(m2hRmsError);
    m2hLogMaxError = log10(m2hMaxError);
    m2hLogMaxAbsError = log10(m2hMaxAbsError);

    % Compute the slope of discrepancy at the decision point
    [indDecision, m2hDiffSlopeAtDecision] = ...
        m3ha_find_decision_point(m2hDiff, 'tVecs', tVecs, ...
                        'FiltWidthMs', filtWidthMs, ...
                        'Itm2hDiffLowerLimit', itm2hDiffLowerLimit, ...
                        'Itm2hDiffLeftBound', itm2hDiffLeftBound, ...
                        'OnlyIfReached', true);

    % Extract the slope value at the decision point
    vSlopeAtDecision = extract_elements(dvdtSmoothed, 'specific', ...
                                        'Index', indDecision - 1);

    % Add or update variable
    simSwpInfo = updatevars(simSwpInfo, m2hRmsError);
    simSwpInfo = updatevars(simSwpInfo, m2hMaxAbsError);
    simSwpInfo = updatevars(simSwpInfo, m2hMaxError);
    simSwpInfo = updatevars(simSwpInfo, m2hMaxRatio);
    simSwpInfo = updatevars(simSwpInfo, m2hLogRmsError);
    simSwpInfo = updatevars(simSwpInfo, m2hLogMaxError);
    simSwpInfo = updatevars(simSwpInfo, m2hLogMaxAbsError);
    simSwpInfo = updatevars(simSwpInfo, m2hMaxLogRatio);
    simSwpInfo = updatevars(simSwpInfo, m2hDiffSlopeAtDecision);
    simSwpInfo = updatevars(simSwpInfo, vSlopeAtDecision);

    % Resave the simulated sweep info table
    writetable(simSwpInfo, simSwpInfoPath, 'WriteRowNames', true);
end

%% Read LTS and open probability discrepancy info
if plotOpenProbabilityFlag || findSpecialCasesFlag
    % Display message
    fprintf('Reading LTS and open probability discrepancy info ... \n');

    % Make sure the simulated sweep info table exists
    if ~isfile(simSwpInfoPath)
        error('Save a simulated sweep info table first!');
    end

    % Create a rank string
    rankStrOP = ['rank', create_label_from_sequence(rankNumsOpenProbability)];

    % Decide on the suffix
    openProbSuffixes = strcat(simStr, '_', opdMeasureStrs, '_vs_hasLTS');

    % Construct figure paths
    openProbPathBases = fullfile(outFolder, strcat(prefix, '_', rankStrOP, ...
                                                    '_', openProbSuffixes));

    % Extract the corresponding cell names
    cellNamesOP = ...
        extract_from_rank_numbers(rankNumsOpenProbability, rankDirectory);

    % Read the simulated sweep info table
    simSwpInfo = readtable(simSwpInfoPath, 'ReadRowNames', true);

    % Restrict to the cells of interest
    simSwpInfoOP = simSwpInfo(contains(simSwpInfo.fileBase, cellNamesOP), :);

    % Read file bases 
    fileBasesOP = simSwpInfoOP.fileBase;

    % Extract cell names
    cellNamesEachFile = m3ha_extract_cell_name(fileBasesOP);

    % Read the LTS peak times
    ltsPeakTime = extract_vars(simSwpInfoOP, 'ltsPeakTime');

    % Read the open probability discrepancy measures
    opdMeasureValues = extract_vars(simSwpInfoOP, opdMeasureStrs);

    % Set infinite values as NaN
    opdMeasureValues = restrict_values(opdMeasureValues, 'Inf2NaN', true);

    % Determine whether each sweep has an LTS
    noLts = isnan(ltsPeakTime);
end

%% Plot all individual plots
if plotIndividualFlag
    % Display message
    fprintf('Plotting all individual fits ... \n');

    % Locate all simulated .out files
    [~, allSimOutPaths] = all_files('Directory', outFolder, ...
                                    'Recursive', true, 'Extension', 'out');

    % Median filter recorded trace only if not using HH2.mod
    toMedianFilter = ~useHH;

    % Plot simulated traces
    array_fun(@(a) m3ha_plot_individual(a, individualXLimits, ...
                            individualYLimits, individualYTickLocs, ...
                            toMedianFilter), ...
            allSimOutPaths);
end

%% Plot all essential plots
if plotEssentialFlag
    % Display message
    fprintf('Plotting all essential plots ... \n');

    % Locate all simulated .out files
    [~, allSimOutPaths] = all_files('Directory', outFolder, ...
                                    'Recursive', true, 'Extension', 'out');

    % Plot simulated traces
    array_fun(@(a) m3ha_plot_essential(a, overlappedXLimits, ...
                            essentialYLimits, essentialYTickLocs, ...
                            opdThreshold), ...
            allSimOutPaths);
end

%% Plot open probability discrepancy against LTS presence
if plotOpenProbabilityFlag
    % Display message
    fprintf('Plotting open probability discrepancies ... \n');

    % Compute and plot all open probability discrepancy measures
    cellfun(@(values, pathBase, threshold, yLimits, yLabel) ...
                compute_and_plot_opd(values, noLts, ...
                        cellNamesEachFile, cellNamesOP, pathBase, prefix, ...
                        openProbFigWidth, openProbFigHeight, figTypes, ...
                        threshold, yLimits, yLabel), ...
            opdMeasureValues, openProbPathBases, ...
            opdThresholds, yLimitsOpd, yLabelsOpd);
end 

%% Find and copy special cases
if findSpecialCasesFlag
    % Display message
    fprintf('Finding and copying special cases ... \n');

    % Create special cases directories
    specialDir = fullfile(outFolder, 'special_cases');
    [falseNegDir, falsePosDir] = ...
        argfun(@(s) fullfile(specialDir, s), ...
                'highDiscrepancyNoLts', 'lowDiscrepancyHasLts');
    check_dir({specialDir, falseNegDir, falsePosDir});

    % Determine whether it is a special case
    isFalseNeg = noLts & (logOpenProbabilityDiscrepancy >= logOpdThreshold);
    isFalsePos = ~noLts & (logOpenProbabilityDiscrepancy < logOpdThreshold);

    % Find corresponding file base
    [fileBasesFalseNeg, fileBasesFalsePos] = ...
        argfun(@(i) fileBasesOP(i), isFalseNeg, isFalsePos);

    % Add '_' to file base
    [fileBasesModFalseNeg, fileBasesModFalsePos] = ...
        argfun(@(fb) strcat(fb, '_'), fileBasesFalseNeg, fileBasesFalsePos);

    % Find corresponding essential plots
    [~, essentialPathsFalseNeg] = ...
        find_matching_files(fileBasesModFalseNeg, 'Directory', outFolder, ...
                                'Recursive', true, 'Suffix', 'essential', ...
                                'Extension', 'png');
    [~, essentialPathsFalsePos] = ...
        find_matching_files(fileBasesModFalsePos, 'Directory', outFolder, ...
                                'Recursive', true, 'Suffix', 'essential', ...
                                'Extension', 'png');

    % Copy into special cases directories
    copy_into(essentialPathsFalseNeg, falseNegDir);
    copy_into(essentialPathsFalsePos, falsePosDir);
end

%% Compute cell info table
if computeCellInfoTableFlag
    % Display message
    fprintf('Computing summary cell info table ... \n');

    % Make sure the simulated sweep info table exists
    if ~isfile(simSwpInfoPath)
        error('Save a simulated sweep info table first!');
    end

    % Create the cell info table by combining the parameter tables
    simCellInfoTable = ...
        combine_param_tables(paramPaths, 'NewRowNames', cellNames);

    % Read the simulated sweep info table
    simSwpInfo = readtable(simSwpInfoPath, 'ReadRowNames', true);

    % Extract the sweep names
    swpNames = simSwpInfo.Properties.RowNames;
    cellName = m3ha_extract_cell_name(swpNames);

    % Determine whether each sweep has an LTS or burst
    ltsProbability = ~isnan(simSwpInfo.ltsPeakTime);
    burstProbability = ~isnan(simSwpInfo.burstTime);

    % Add ltsProbability & burstProbability
    simSwpInfoOfInterest = ...
        addvars(simSwpInfo, ltsProbability, burstProbability, 'Before', 1);

    % Restrict to measures of interest
    simSwpInfoOfInterest = simSwpInfoOfInterest(:, measuresOfInterestOrig);

    % Add cell names
    simSwpInfoOfInterest = ...
        addvars(simSwpInfoOfInterest, cellName, 'Before', 1);

    % Average all columns, grouped by cellName
    %   Note: groupsummary() is available since R2018a
    cellMeasureTable = ...
        groupsummary(simSwpInfoOfInterest, 'cellName', @nanmean);

    % Create the variable names returned by groupsummary()
    summaryVarNames = strcat('fun1_', measuresOfInterestOrig);

    % Rename variables
    cellMeasureTable = renamevars_custom(cellMeasureTable, ...
                            summaryVarNames, measuresOfInterestNew);

    % Make cell names row names
    cellMeasureTable.Properties.RowNames = cellMeasureTable.cellName;

    % Join tables
    simCellInfoTable = join(simCellInfoTable, cellMeasureTable, 'Keys', 'Row');

    % Write to the table
    writetable(simCellInfoTable, simCellInfoPath, 'WriteRowNames', true);
end

%% Plot correlations
if plotCorrelationsFlag
    % Display message
    fprintf('Plotting correlation plots ... \n');

    % Make sure cell info table exists
    if ~isfile(simCellInfoPath)
        error('Save a simulated cell info table first!');
    end

    % Read the cell info table
    simCellInfoTable = readtable(simCellInfoPath, 'ReadRowNames', true);

    % Collect all x & y variables
    paramVecs = extract_vars(simCellInfoTable, paramsOfInterest);
    measureVecs = extract_vars(simCellInfoTable, measuresOfInterestNew);

    % Plot all correlations between measures and params
    plot_all_correlations(paramVecs, measureVecs, ...
                        paramsOfInterest, measuresOfInterestNew, ...
                        paramIslog, measuresIsLogNew, ...
                        outFolderCorrMeas2Param);

    % Plot all correlations between params and params
    plot_all_correlations(paramVecs, paramVecs, ...
                        paramsOfInterest, paramsOfInterest, ...
                        paramIslog, paramIslog, ...
                        outFolderCorrParam2Param);
end

%% Plot violin plots
if plotViolinPlotsFlag
    % Make sure the simulated sweep info table exists
    if ~isfile(simSwpInfoPath)
        error('Save a simulated sweep info table first!');
    end

    % Read the simulated sweep info table
    simSwpInfo = readtable(simSwpInfoPath, 'ReadRowNames', true);

    % Compute and plot violin plots
    cellfun(@(conditionLabel2D, pCond2D, gCond2D) ...
            m3ha_compute_and_plot_violin(outFolder, ...
                    pharmLabelsShort, conditionLabel2D, ...
                    'SwpInfo', simSwpInfo, 'DataMode', dataMode, ...
                    'PharmConditions', pCond2D, 'GIncrConditions', gCond2D, ...
                    'RowsToPlot', measuresOfInterest, ...
                    'OutFolder', outFolder), ...
            conditionLabels2D, pConds2D, gConds2D);
end

%% Plot bar plots
if plotBarPlotsFlag
    % Make sure the simulated sweep info table exists
    if ~isfile(simSwpInfoPath)
        error('Save a simulated sweep info table first!');
    end

    % Read the simulated sweep info table
    simSwpInfo = readtable(simSwpInfoPath, 'ReadRowNames', true);
    % Compute statistics if not done already
    if ~isfile(stats3dPath)
        % Load sweep info
        simSwpInfo = readtable(simSwpInfoPath, 'ReadRowNames', true);

        % Compute statistics for all features
        disp('Computing statistics for 3D bar plots ...');
        statsTable = m3ha_compute_statistics('SwpInfo', simSwpInfo, ...
                                                'PharmConditions', pCond3D, ...
                                                'GIncrConditions', gCond3D, ...
                                                'DataMode', dataMode);
        % Generate a condition label
        pharmLabels = pharmLabelsLong;
        conditionLabel = conditionLabel3D;

        % Save stats table
        save(stats3dPath, 'statsTable', 'pharmLabels', ...
                        'gIncrLabels', 'conditionLabel', '-v7.3');
    end

    % Plot all 3D bar plots
    m3ha_plot_bar3(stats3dPath, 'RowsToPlot', measuresOfInterest);
end

% Archive all scripts for this run
if archiveScriptsFlag
    archive_dependent_scripts(mfilename, 'OutFolder', outFolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function myTable = updatevars (myTable, varValue)
%% Replace a variable in a table or add it if it doesn't exist
% TODO: Pull out as its own function

varName = inputname(2);
if is_field(myTable, varName)
    myTable.(varName) = varValue;
else
    myTable = addvars(myTable, varValue, 'NewVariableNames', varName);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cellNames, iterStrs] = ...
                extract_from_rank_numbers (rankNumsToUse, rankDirectory)

% Create rank number prefixes
rankPrefixes = create_labels_from_numbers(rankNumsToUse, ...
                                    'Prefix', 'rank_', 'Suffix', '_');

% Find png files matching the rank prefixes
[~, pngPaths] = find_matching_files(rankPrefixes, 'PartType', 'prefix', ...
                        'Directory', rankDirectory, 'Extension', 'png', ...
                        'ExtractDistinct', false);

% Extract the cell names
cellNames = m3ha_extract_cell_name(pngPaths, 'FromBaseName', true);

% Extract the iteration numbers
iterStrs = m3ha_extract_iteration_string(pngPaths, 'FromBaseName', true);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m3ha_plot_individual (simOutPath, xLimits, yLimits, yTickLocs, ...
                                toMedianFilter)
%% Plot an essential plot for a simulation

% Extract the sweep name
sweepName = m3ha_extract_sweep_name(simOutPath, 'FromBaseName', true);

% Make the sweep the the figure title
figTitle = replace(['Fits for ', sweepName], '_', '\_');

% Create figure names
figPathBase = replace(simOutPath, '.out', '_individual.png');

% Create the figure
fig = set_figure_properties('AlwaysNew', true);

% Plot essential traces
handles = ...
    m3ha_plot_simulated_traces('PlotType', 'individual', ...
                            'ToMedianFilter', toMedianFilter, ...
                            'FileNames', simOutPath, ...
                            'XLimits', xLimits, 'YLimits', yLimits, ...
                            'YTickLocs', yTickLocs, ...
                            'FigHandle', fig, 'FigTitle', figTitle);

% Save original figure
save_all_figtypes(fig, figPathBase, 'png');

% Close figure
close(fig);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m3ha_plot_essential (simOutPath, xLimits, ...
                                yLimits, yTickLocs, opdThreshold)
%% Plot an essential plot for a simulation

% Extract the sweep name
sweepName = m3ha_extract_sweep_name(simOutPath, 'FromBaseName', true);

% Make the sweep the the figure title
figTitle = replace(['Essential plots for ', sweepName], '_', '\_');

% Create figure names
figPathBase = replace(simOutPath, '.out', '_essential.png');

% Create the figure
fig = set_figure_properties('AlwaysNew', true);

% Plot essential traces
handles = ...
    m3ha_plot_simulated_traces('PlotType', 'essential', ...
                            'FileNames', simOutPath, ...
                            'XLimits', xLimits, 'YLimits', yLimits, ...
                            'YTickLocs', yTickLocs, ...
                            'FigHandle', fig, 'FigTitle', figTitle);

% Get all subplots
subPlots = handles.subPlots;

% Add a threshold line
subplot(subPlots(5)); hold on;
plot_horizontal_line(opdThreshold, 'ColorMap', 'DarkGreen', ...
                        'LineStyle', ':', 'LineWidth', 1);

% Save original figure
save_all_figtypes(fig, figPathBase, 'png');

% Close figure
close(fig);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_all_correlations (xVarVecs, yVarVecs, ...
                                xVarNames, yVarNames, ...
                                xVarIsLog, yVarIsLog, outFolder)
%% Plots all correlations between two sets of vectors
% TODO: Pull out as its own function
% Requires:
% TODO
% force_column_cell
% all_ordered_pairs
% plot_correlation

% Force as column cell arrays
[xVarVecs, yVarVecs] = argfun(@force_column_cell, xVarVecs, yVarVecs);

% Create all ordered pairs
allVecPairs = all_ordered_pairs({xVarVecs, yVarVecs});
allVarPairs = all_ordered_pairs({xVarNames, yVarNames});
allIsLogPairs = all_ordered_pairs({xVarIsLog, yVarIsLog});
allFigPaths = cellfun(@(b) fullfile(outFolder, ...
                                    [b{2}, '_vs_', b{1}, '.png']), ...
                    allVarPairs, 'UniformOutput', false);

% Plot correlations
% cellfun(@(a, b, c, d) plot_correlation(a{1}, a{2}, b{1}, b{2}, ...
%                                           c(1), c(2), d), ...
%         allVecPairs, allVarPairs, allIsLogPairs, allFigPaths);
% array_fun(@(a, b, c, d) plot_correlation(a{1}, a{2}, b{1}, b{2}, ...
%                                           c(1), c(2), d), ...
%         allVecPairs, allVarPairs, allIsLogPairs, allFigPaths);
parfor iPair = 1:numel(allVecPairs)
    vecPair = allVecPairs{iPair};
    varPair = allVarPairs{iPair};
    figPath = allFigPaths{iPair};
    isLogPair = allIsLogPairs{iPair};
    plot_correlation(vecPair{1}, vecPair{2}, varPair{1}, varPair{2}, ...
                    isLogPair(1), isLogPair(2), figPath);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = plot_correlation (xVec, yVec, xVar, yVar, ...
                                        xIsLog, yIsLog, figPath)
%% Plots correlation between two vectors

% Add correlation value to figure name
addCorrValueToFigName = true;
addCovValueToFigName = true;

% Find directory of figure path
figDir = extract_fileparts(figPath, 'directory');
if isempty(figDir)
    figDir = pwd;
end

% Remove the pairs of values with NaN
neitherNaN = ~isnan(xVec) & ~isnan(yVec);
xVecNoNaN = xVec(neitherNaN);
yVecNoNaN = yVec(neitherNaN);

% Convert to log scale if requested 
if xIsLog
    xVecToCorr = log10(xVecNoNaN);
else
    xVecToCorr = xVecNoNaN;
end
if yIsLog
    yVecToCorr = log10(yVecNoNaN);
else
    yVecToCorr = yVecNoNaN;
end

% Compute the correlation coefficient
corrValue = corr2(xVecToCorr, yVecToCorr);

% Compute the covariance value
covMatrix = cov(xVecToCorr, yVecToCorr);
covValue = covMatrix(1, 2);

% Add covariance value to figure name if requested
if addCovValueToFigName
    figNameOld = extract_fileparts(figPath, 'name');
    figNameNew = sprintf(['cov_%.2f_%s'], covValue, figNameOld);
    figPath = replace(figPath, figNameOld, figNameNew);
end

% Add covariance value to figure name if requested
if addCorrValueToFigName
    figNameOld = extract_fileparts(figPath, 'name');
    figNameNew = sprintf(['corr_%.2f_%s'], corrValue, figNameOld);
    figPath = replace(figPath, figNameOld, figNameNew);
end

% Decide on the text color
if abs(corrValue) > 0.6 && abs(corrValue) ~= 1
    isSignificant = true;
    textColor = 'r';
else
    isSignificant = false;
    textColor = 'k';
end

% Create the figure
fig = set_figure_properties('AlwaysNew', true);

% Plot scatter plot
scatter(xVecNoNaN, yVecNoNaN, 'o', 'LineWidth', 2);

% Change axis scale if requested
if xIsLog || yIsLog
    [xScale, yScale] = argfun(@islog2scale, xIsLog, yIsLog);
    set(gca, 'XScale', xScale, 'YScale', yScale);
end

% Show correlation coefficient
text(0.1, 0.95, ...
    ['Correlation coefficient: ', num2str(corrValue, 3)], ...
    'Units', 'normalized', 'Color', textColor); 

% Titles and axis labels
title(['Correlation of ', yVar, ' vs. ', xVar]);
xlabel(xVar);
ylabel(yVar);

% Save original figure
save_all_figtypes(fig, figPath, 'png');

if isSignificant
    % Create directory for significant
    sigDir = fullfile(figDir, 'significant');
    check_dir(sigDir);

    % Create path to significant figure
    figName = extract_fileparts(figPath, 'name');
    figPathSig = fullfile(sigDir, figName);

    % Save original figure
    save_all_figtypes(fig, figPathSig, 'png');
end

% Close figure
close(fig);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compute_and_plot_opd (opdMeasure, noLts, ...
                    cellNamesEachFile, cellNames, pathBase, ...
                    prefix, figWidth, figHeight, figTypes, ...
                    threshold, yLimits, yLabel)

% Find the indices for each cell with or without LTS
[indEachCellWithLTS, indEachCellWithNoLTS] = ...
    argfun(@(condition) ...
                cellfun(@(c) condition & strcmp(cellNamesEachFile, c), ...
                        cellNames, 'UniformOutput', false), ...
            ~noLts, noLts);

% Group by cell
[logOpdByCellWithLTS, logOpdByCellWithNoLTS] = ...
    argfun(@(indEachCell) cellfun(@(ind) opdMeasure(ind), indEachCell, ...
                                'UniformOutput', false), ...
            indEachCellWithLTS, indEachCellWithNoLTS);

% Extract the log open probability discrepancy for all traces
[logOpdWithLTSAllTraces, logOpdWithNoLTSAllTraces] = ...
    argfun(@(ind) opdMeasure(ind), ...
            ~noLts, noLts);

% averageFun = @(x) compute_weighted_average(x, 'AverageMethod', 'geometric', ...
%                                             'IgnoreNan', true);
averageFun = @(x) compute_weighted_average(x, 'AverageMethod', 'arithmetic', ...
                                            'IgnoreNan', true);

% Compute the average (arithmetic mean) log open probability discrepancy
%   for each cell with or without LTS
[logOpdWithLTS, logOpdWithNoLTS] = ...
    argfun(@(x) cellfun(averageFun, x), ...
            logOpdByCellWithLTS, logOpdByCellWithNoLTS);

% Separate into two groups
twoGroups = {logOpdWithNoLTS; logOpdWithLTS};
twoGroupsByCell = {logOpdByCellWithNoLTS; logOpdByCellWithLTS};
twoGroupsAllTraces = {logOpdWithNoLTSAllTraces; logOpdWithLTSAllTraces};

% Test for differences
opdMeasureDifferences = test_difference(twoGroups, 'IsPaired', true, ...
                                        'SaveFlag', true, 'FileBase', pathBase);

% Contruct path base for all traces plot
pathBaseAllTraces = [pathBase, '_all_traces'];

%{
% Plot violin plots
plot_open_probability_discrepancy(twoGroups, 'violin', ...
                pathBase, prefix, figWidth, figHeight, figTypes, ...
                threshold, yLimits, yLabel, {});
plot_open_probability_discrepancy(twoGroupsAllTraces, 'violin', ...
                pathBaseAllTraces, prefix, figWidth, figHeight, figTypes, ...
                threshold, yLimits, yLabel, {});

% Plot jitter plots
plot_open_probability_discrepancy(twoGroupsByCell, 'byCell', ...
                pathBaseAllTraces, prefix, figWidth, figHeight, figTypes, ...
                threshold, yLimits, yLabel, cellNames);
%}

%% Do this for each cell
% Contruct path base for traces from each cell
pathBaseEachCell = strcat(pathBase, '_', cellNames);

% Reorganize cell arrays
twoGroupsEachCell = regroup_cell_of_cells(twoGroupsByCell);

% Test for differences
cellfun(@(twoGroups, pathBase) ...
            test_difference(twoGroups, 'IsPaired', false, ...
                        'SaveFlag', true, 'FileBase', pathBase), ...
        twoGroupsEachCell, pathBaseEachCell);

% Plot violin plots
cellfun(@(twoGroups, pathBase, figTitle) ...
            plot_open_probability_discrepancy(twoGroups, 'violin', ...
                    pathBase, figTitle, figWidth, figHeight, 'png', ...
                    threshold, yLimits, yLabel, {}), ...
        twoGroupsEachCell, pathBaseEachCell, cellNames);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_open_probability_discrepancy (twoGroups, plotType, ...
                        pathBase, figTitle, figWidth, figHeight, figTypes, ...
                        logOpdThreshold, yLimits, yLabel, cellLabels)

% Hard-coded parameters
xTickLabels = {'No LTS', 'With LTS'};
colorMapViolin = {'Black', 'DarkGreen'};

% Create figure
fig = set_figure_properties('AlwaysNew', true);

% Plot violin plot
switch plotType
    case 'violin'
        violins = plot_violin(twoGroups, 'ColorMap', colorMapViolin, ...
                            'XTickLabels', xTickLabels, 'YLabel', yLabel);
    case 'byCell'
        % Add to path base
        pathBase = [pathBase, '_byCell'];

        % Plot groups as a grouped jitter plot
        jitters = plot_grouped_jitter(twoGroups, 'XTickLabels', xTickLabels, ...
                            'YLabel', yLabel, 'GroupingLabels', cellLabels);
    otherwise
        error('plotType unrecognized!');
end

% Set y axis limits
if ~isempty(yLimits)
    ylim(yLimits);
end

% Plot threshold
if ~isempty(logOpdThreshold)
    hold on
    plot_horizontal_line(logOpdThreshold, 'ColorMap', 'r', 'LineWidth', 2);
end

% Create a title
title(replace(figTitle, '_', '\_'));

% Save the figure
pathBaseOrig = [pathBase, '_orig'];
save_all_figtypes(fig, pathBaseOrig, 'png');

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                'Width', figWidth, 'Height', figHeight, 'RemoveLegends', true);

% Save figure
save_all_figtypes(fig, pathBase, figTypes);

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%% Compile
% Display message
fprintf('Compiling all .mod files ... \n');

% Compile or re-compile .mod files in the fitting directory
compile_mod_files(fitDirectory);

%% Select recorded data
% Display message
fprintf('Selecting recorded sweeps for all cells ... \n');

% Locate the data directory
dataDir = fullfile(parentDirectory, dataDirName);

% Construct full paths to other directories used 
%   and previously analyzed results under dataDir
[matFilesDir, specialCasesDir] = ...
    argfun(@(x) fullfile(dataDir, x), matFilesDirName, specialCasesDirName);

% Select the sweep indices that will be simulated
swpInfo = m3ha_select_sweeps('SwpInfo', swpInfo, 'DataMode', dataMode, ...
                                'CasesDir', specialCasesDir);

% Select the raw traces to import for each cell to fit
[fileNamesToFit, rowConditionsToFit] = ...
    arrayfun(@(x) m3ha_select_raw_traces(rowmodeAcrossTrials, ...
                    columnMode, attemptNumberAcrossTrials, ...
                    x, swpInfo, cellInfo), ...
            cellNamesToFit, 'UniformOutput', false);

% Separate into two groups
twoGroups = {logOpenProbabilityDiscrepancy(noLts); ...
                logOpenProbabilityDiscrepancy(~noLts)};

logOpenProbabilityDiscrepancyDifferences = test_difference(twoGroups);

% Compute the average (geometric mean) open probability discrepancy
%   for each cell with or without LTS
[opdWithLTS, opdWithNoLTS] = ...
    argfun(@(indEachCell) ...
            cellfun(@(ind) compute_weighted_average(...
                    openProbabilityDiscrepancy(ind), ...
                    'AverageMethod', 'geometric'), indEachCell), ...
            indEachCellWithLTS, indEachCellWithNoLTS);

% Compute the average (arithmetic mean) log open probability discrepancy
%   for each cell with or without LTS
[logOpdWithLTS, logOpdWithNoLTS] = ...
    argfun(@(indEachCell) ...
            cellfun(@(ind) compute_weighted_average(...
                    opdMeasure(ind), ...
                    'AverageMethod', 'arithmetic', ...
                    'IgnoreNan', true), indEachCell), ...
            indEachCellWithLTS, indEachCellWithNoLTS);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
