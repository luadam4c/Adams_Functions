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
%       cd/compute_rms_error.m
%       cd/copy_into.m
%       cd/create_time_stamp.m
%       cd/create_labels_from_numbers.m
%       cd/create_label_from_sequence.m
%       cd/extract_columns.m
%       cd/find_matching_files.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
%       cd/is_field.m
%       cd/load_neuron_outputs.m
%       cd/m3ha_compute_statistics.m
%       cd/m3ha_extract_cell_name.m
%       cd/m3ha_extract_iteration_string.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_load_sweep_info.m
%       cd/m3ha_plot_bar3.m
%       cd/m3ha_plot_violin.m
%       cd/plot_violin.m
%       cd/print_cellstr.m
%       cd/renamevars.m
%       cd/vertcat_spreadsheets.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-12-11 Created by Adam Lu
% 2019-12-26 Completed
% 2019-12-27 Added HH channels
% 

%% Hard-coded parameters
% Flags
chooseBestNeuronsFlag = false; %true;
simulateFlag = false; %true;
combineFeatureTablesFlag = false; %true;
computeOpenProbabilityFlag = false; %true;
plotOpenProbabilityFlag = false; %true;
plotViolinPlotsFlag = true;
plotBarPlotsFlag = false; %true;
archiveScriptsFlag = true;

% Simulation parameters
useHH = true;           % whether to use Hudgin-Huxley Na+ and K+ channels
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
ltsParamsSuffix = '_ltsParams';
simLtsParamsSuffix = strcat(simStr, '_ltsParams');
simSwpInfoSuffix = strcat(simStr, '_swpInfo');
openProbSuffix = strcat(simStr, '_openProbabilityDiscrepancy_vs_hasLTS');
simOutExtension = 'out';

TIME_COL_REC = 1;

% Column numbers for simulated data
%   Note: Must be consistent with singleneuron4compgabab.hoc
IT_M_DEND2 = 47;
IT_MINF_DEND2 = 48;
IT_H_DEND2 = 49;
IT_HINF_DEND2 = 50;

% Note: The following must be consistent with m3ha_parse_dclamp_data.m
condVarStrs = {'cellidrow', 'prow', 'vrow', 'grow', 'swpnrow', ...
                'gabab_amp', 'gabab_Trise', 'gabab_TfallFast', ...
                'gabab_TfallSlow', 'gabab_w'};
pharmAll = [1; 2; 3; 4];          
pharmLabelsLong = {'{\it s}-Control', '{\it s}-GAT1 Block', ...
                    '{\it s}-GAT3 Block', '{\it s}-Dual Block'};
pharmLabelsShort = {'{\it s}-Con', '{\it s}-GAT1', ...
                    '{\it s}-GAT3', '{\it s}-Dual'};
if dataMode == 0 || dataMode == 3
    gIncrAll = [25; 50; 100; 200; 400; 800];
    gIncrLabels = {'25%', '50%', '100%', '200%', '400%', '800%'};
elseif dataMode == 1 || dataMode == 2
    gIncrAll = [100; 200; 400];
    gIncrLabels = {'100%', '200%', '400%'};
end
conditionLabel2D = 'pharm_1-4_gincr_200_sim';
pCond2D = num2cell(pharmAll);
gCond2D = 200;
stats2dSuffix = strcat(simStr, '_', conditionLabel2D, '_stats.mat');
conditionLabel3D = 'pharm_1-4_gincr_all_sim';
pCond3D = num2cell(pharmAll);
gCond3D = num2cell(gIncrAll);
stats3dSuffix = strcat(simStr, '_', conditionLabel3D, '_stats.mat');

% Plot settings
% Note: must be consistent with m3ha_compute_statistics.m
measuresOfInterest = {'ltsAmplitude'; 'ltsMaxSlope'; ...
                    'ltsConcavity'; 'ltsProminence'; ...
                    'ltsWidth'; 'ltsOnsetTime'; 'ltsTimeJitter'; ...
                    'ltsProbability'; 'spikesPerLts'; ...
                    'spikeMaxAmp'; 'spikeMinAmp'; ...
                    'spikeFrequency'; 'spikeAdaptation'
                    'burstOnsetTime'; 'burstTimeJitter'; ...
                    'burstProbability'; 'spikesPerBurst'};
openProbFigWidth = 5;       % (cm)
openProbFigHeight = 3;      % (cm)

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

outFolder = '20200208_population_rank1-23_dataMode1_attemptNumber3';
figTypes = {'png', 'epsc'};
rankDirName = '20200207_ranked_manual_singleneuronfitting0-102';
rankNumsToUse = 1:23;
ipscrWindow = [2000, 4800];     % only simulate up to that time
fitWindowIpscr = [3000, 4800];  % the time window (ms) where all 
                                %   recorded LTS would lie
prefix = '';

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

% Check if output folder exists
check_dir(outFolder);

% Decide on output prefix
if isempty(prefix)
    % Extract output folder base name
    prefix = extract_fileparts(outFolder, 'dirbase');
end

% Construct full paths
simSwpInfoPath = fullfile(outFolder, [prefix, '_', simSwpInfoSuffix, '.csv']);
stats2dPath = fullfile(outFolder, [prefix, '_', stats2dSuffix, '.mat']);
stats3dPath = fullfile(outFolder, [prefix, '_', stats3dSuffix, '.mat']);
openProbPathBase = fullfile(outFolder, [prefix, '_', openProbSuffix]);
openProbPathBaseOrig = [openProbPathBase, '_orig'];

%% Choose the best cells and the best parameters for each cell
if chooseBestNeuronsFlag
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

%% Simulate
if simulateFlag
    % Decide on candidate parameters files
    [~, paramPaths] = all_files('Directory', outFolder, 'Suffix', 'params', ...
                                'Recursive', false);

    % Remove the ltsParams files
    paramPaths = paramPaths(~contains(paramPaths, 'ltsParams'));

    % Extract the cell names
    cellNames = m3ha_extract_cell_name(paramPaths);

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
                    'UseHH', useHH, 'AttemptNumber', attemptNumber, ...
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
    simSwpInfo = renamevars(simSwpInfo, 'fileBase', 'simFileBase');

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
    simData = load_neuron_outputs('FileNames', simOutPaths, 'tVecs', tVecsRec);

    % Extract vectors from simulated data
    [m, minf, h, hinf] = ...
        extract_columns(simData, [IT_M_DEND2, IT_MINF_DEND2, ...
                                    IT_H_DEND2, IT_HINF_DEND2]);

    % Force as matrices
    [m, minf, h, hinf] = argfun(@force_matrix, m, minf, h, hinf);

    % Compute m2h and minf2hinf
    m2h = (m .^ 2) .* h;
    minf2hinf = (minf .^ 2) .* hinf;

    % Compute the rms error between m2h and minf2hinf
    m2hRmsError = compute_rms_error(m2h, minf2hinf, ...
                                            'ForceColumnOutput', true);

    % Compute the maximum error between m2h and minf2hinf
    m2hMaxError = ...
        force_column_vector(max(abs(m2h - minf2hinf), [], 1));

    % Compute the maximum ratio between m2h and minf2hinf
    m2hRatioRaw = m2h ./ minf2hinf;
    m2hRatioRaw(isinf(m2hRatioRaw)) = NaN;
    m2hMaxRatio = force_column_vector(max(m2hRatioRaw, [], 1));

    % Add or replace variable
    simSwpInfo = replacevars(simSwpInfo, m2hRmsError);
    simSwpInfo = replacevars(simSwpInfo, m2hMaxError);
    simSwpInfo = replacevars(simSwpInfo, m2hMaxRatio);

    % Resave the simulated sweep info table
    writetable(simSwpInfo, simSwpInfoPath, 'WriteRowNames', true);
end

%% Plot open probability discrepancy against LTS presence
if plotOpenProbabilityFlag
    % Make sure the simulated sweep info table exists
    if ~isfile(simSwpInfoPath)
        error('Save a simulated sweep info table first!');
    end

    % Display message
    fprintf('Plotting open probability discrepancies ... \n');

    % Read the simulated sweep info table
    simSwpInfo = readtable(simSwpInfoPath, 'ReadRowNames', true);

    % Read the LTS peak times
    if ~is_field(simSwpInfo, 'ltsPeakTime')
        error('ltsPeakTime does not exist!');
    else
        ltsPeakTime = simSwpInfo.ltsPeakTime;
    end

    % Read the LTS peak times
    if ~is_field(simSwpInfo, 'm2hMaxRatio')
        error('m2hMaxRatio does not exist yet!');
    else
        % openProbabilityDiscrepancy = simSwpInfo.m2hRmsError;
        % openProbabilityDiscrepancy = simSwpInfo.m2hMaxError;
        openProbabilityDiscrepancy = simSwpInfo.m2hMaxRatio;
    end

    % Determine whether each sweep has an LTS
    noLts = isnan(ltsPeakTime);

    % Separate into two groups
    twoGroups = {openProbabilityDiscrepancy(noLts); ...
                    openProbabilityDiscrepancy(~noLts)};

    % Create figure
    fig = set_figure_properties('AlwaysNew', true);

    % Plot violin plot
    violins = plot_violin(twoGroups, 'XTickLabels', {'No LTS', 'With LTS'}, ...
                            'YLabel', 'max(m^2h / m_{inf}^2h_{inf})');
                            % 'YLabel', 'max(abs(m^2h - m_{inf}^2h_{inf}))');
                            % 'YLabel', 'rms(m^2h - m_{inf}^2h_{inf})');

    % Create a title
    title(sprintf('Open Probability Discrepancy for %s', ...
            replace(prefix, '_', '\_')));

    % Set y axis to be on a log scale
    set(gca, 'YScale', 'log');

    % Save the figure
    save_all_figtypes(fig, openProbPathBaseOrig, 'png');

    % Update figure for CorelDraw
    update_figure_for_corel(fig, 'Units', 'centimeters', ...
                    'Width', openProbFigWidth, 'Height', openProbFigHeight);

    % Save figure
    save_all_figtypes(fig, openProbPathBase, figTypes);

end 

%% Plot violin plots
if plotViolinPlotsFlag
    % Make sure the simulated sweep info table exists
    if ~isfile(simSwpInfoPath)
        error('Save a simulated sweep info table first!');
    end

    % Read the simulated sweep info table
    simSwpInfo = readtable(simSwpInfoPath, 'ReadRowNames', true);
    % Compute statistics if not done already
    if ~isfile(stats2dPath)
        % Load sweep info
        simSwpInfo = readtable(simSwpInfoPath, 'ReadRowNames', true);

        % Compute statistics for all features
        disp('Computing statistics for violin plots ...');
        statsTable = m3ha_compute_statistics('SwpInfo', simSwpInfo, ...
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

    % Plot all 2D violin plots
    m3ha_plot_violin(stats2dPath, 'RowsToPlot', measuresOfInterest);
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

function myTable = replacevars (myTable, varValue)
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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
