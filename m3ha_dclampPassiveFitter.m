function [passiveParams, fitResults, fitObject, goodnessOfFit, ...
            algorithmInfo, decision, allResults] = ...
                    m3ha_dclampPassiveFitter (fitMode, ...
                        infolder, outFolder, plotpassiveflag, groupmode)
%% Estimates passive parameters for each cell from dclamp data recorded by Mark & Christine
% Usage: [passiveParams, fitResults, fitObject, goodnessOfFit, ...
%           algorithmInfo, decision, allResults] = ...
%                   m3ha_dclampPassiveFitter (fitMode, ...
%                       infolder, outFolder, plotpassiveflag, groupmode)
% Side effects:
%       Saves results in .mat files
%
% Arguments: 
%       fitMode     - 0 - all data
%                   - 1 - all of g incr = 100%, 200%, 400%
%                   - 2 - all of g incr = 100%, 200%, 400% 
%                   but exclude cell-pharm-g_incr sets containing problematic sweeps
%       infolder    - (opt) the directory that contains the matfile to read
%                   must be a directory
%                   default == //media/adamX/m3ha/data_dclamp/take4/
%       outFolder   - (opt) the directory to output graphs 
%                       (different subdirectories will be created for each fitMode)
%                   must be a directory
%                   default == //media/adamX/m3ha/data_dclamp/take4/
%       plotpassiveflag    - (opt) whether to plot graphs
%                   default == true
%       groupmode   must be one of the following:
%                       'cell_actVhold'
%                       'cell'
%                   default == 'cell'
% 
% Requires:
%       "infolder"/dclampdatalog_take4.mat
%       cd/check_dir.m
%       cd/find_passive_params.m
%       cd/find_in_strings.m
%       cd/m3ha_correct_unbalanced_bridge.m
%       cd/m3ha_find_ind_to_fit.m
%       cd/m3ha_locate_homedir.m
%       cd/m3ha_parse_mat.m
%       cd/m3ha_specs_for_fitmode.m
%
% Used by:
%       cd/m3ha_dclampDataExtractor.m

% File History:
% 2016-11-01 Adapted from PlotHistogramsRefineThreshold.m
% 2016-11-08 Changed from using Vhold to using actVhold for binning
% 2016-11-10 Nows saves set info linearly in _set and 2-dimensionally in _cv
% 2016-11-10 Added fn_set & fn_cv
% 2016-12-01 Added groupmode, made default to be grouping by cell only
% 2016-12-04 Changed current pulse response to last just 150 ms 
%               (cprWin is changed from [95, 500] to [95, 260])
% 2016-12-04 Added logheaderSwpInfo && logvariablesSwpInfo
% 2016-12-04 Added Rmemb, Gmemb, Gsoma, Gdend to be saved in params
% 2016-12-04 Added rmse_R && rmse_F, the root-mean-squared errors of the 
%               rising and falling phase, respectively, for each sweep
% 2017-12-21 SpecsForFitmode() -> specs_for_fitmode()
% 2018-10-03 Changed tabs to spaces
% 2018-10-15 Multiple changes
% TODO: Input parser

%% Flags
debugflag = 0;
correctDcStepsFlag = 1;

%% The .mat file for sweep info; assumed to be in infolder
sweepInfoFile = 'dclampdatalog_take4.mat';
initialSlopesFile = 'initial_slopes_nSamplesForPlot_2_threeStdMainComponent.mat';
outSheetName = 'dclampPassiveParams_byCells.xlsx';
outMatNameByCells = 'dclampPassiveLog_byCells.mat';
outMatNameBySets = 'dclampPassiveLog_bySets.mat';

%% The data directory containing .mat data files; assumed to be in infolder
dataDirectoryName = 'matfiles';

%% Fixed parameters used in the experiments
VholdBC = [-62.5; -67.5; -72.5];    
                            % Possible holding potential bin centers 
                            %   (LJP-corrected)
cpWin = [95, 115];          % window in which the current pulse would lie (ms) 
                            %       (Supposed to be 100-110 ms but 
                            %           there will be offset)
cprWin = [95, 260];         % window in which the current pulse response 
                            %   would lie (ms)

%% Headers and variables for plotting
% For parameters
logheaderParams = {'Coefficient of first exponential, C0 (mV)', ...
            'Time constant of first exponential, tau0 (ms)', ...
            'Coefficient of second exponential, C1 (mV)', ...
            'Time constant of second exponential, tau1 (ms)', ...
            'Average current pulse amplitude (pA)', ...
            'Voltage difference at steady state (mV)', ...
            'Input resistance of setup (MOhm)', ...
            'alpha1 = sqrt(tau0/tau1 - 1)', ...
            'Initial guess for electrotonic length, L_init', ...
            'Electrotonic length, L', ...
            'Dendritic to somatic conductance ratio, rho', ...
            'Input resistance of cell, R_N (Mohm)', ...
            'Somatic resistance, R_soma (Mohm)', ...
            'Dendritic resistance, Rdend (Mohm)', ...
            'Input conductance of cell, G_N (uS)', ...
            'Somatic conductance, G_soma (uS)', ...
            'Dendritic conductance, G_dend (uS)', ...
            'Membrane time constant, taum (ms)', ...
            'Specific membrane resistivity (Ohm-cm^2)', ...
            'Radius of soma (um)', ...
            'Diameter of dendrite (um)', ...
            'Radius of dendrite (um)', ...
            'Length of dendrite (um)'};

logvariablesParams = {'C0', ...
            'tau0', ...
            'C1', ...
            'tau1', ...
            'cpa_mean', ...
            'dvss', ...
            'Rinput', ...
            'alpha1', ...
            'L_init', ...
            'L', ...
            'rho', ...
            'Rmemb', ...
            'Rsoma', ...
            'Rdend', ...
            'Gmemb', ...
            'Gsoma', ...
            'Gdend', ...
            'taum', ...
            'Rm', ...
            'rad_soma', ...
            'diam_dend', ...
            'rad_dend', ...
            'length_dend'};

% logheaderSwpInfo && logvariablesSwpInfo depend on groupmode
logheaderSwpInfoCellActVHold = {'Data filename', 'Cell #', 'Set #', ...
            'Rough holding potential (mV)', ...
            'Current pulse amplitude (pA)', 'Pulse width (ms)', ...
            'Overall change in membrane potential recorded (mV)'}; %, ...
%            'RMSE (mV) in the rising phase', 'RMSE (mV) in the falling phase'};

logvariablesSwpInfoCellActVHold = {'dataFileName', 'cellId', 'setNumber', ...
            'holdPotential', ...
            'pulseAmplitude', 'pulseWidth', 'voltageChange'}; %, 
%            'rmse_R_row', 'rmse_F_row'};

logheaderSwpInfoCell = {'Data filename', 'Cell #', 'Set #', ...
            'Current pulse amplitude (pA)', 'Pulse width (ms)', ...
            'Overall change in membrane potential recorded (mV)'}; %, ...
%            'RMSE (mV) in the rising phase', 'RMSE (mV) in the falling phase'};

logvariablesSwpInfoCell = {'dataFileName', 'cellId', 'setNumber', ...
            'pulseAmplitude', 'pulseWidth', 'voltageChange'}; %, 
%            'rmse_R_row', 'rmse_F_row'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% If not compiled, add directories to search path for required functions
if ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for find_passive_params.m, etc.
    addpath(fullfile(functionsDirectory, 'Adams_Functions')); 
end

% Locate the home directory
homeDirectory = m3ha_locate_homedir;

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Check arguments
% TODO
if nargin < 1
    error('A fitMode is required, type ''help m3ha_dclampPassiveFitter'' for usage');
elseif isempty(fitMode) || ~isnumeric(fitMode) || ~(fitMode == 0 || fitMode == 1 || fitMode == 2)
    error('fitMode out of range!, type ''help m3ha_dclampPassiveFitter'' for usage');
elseif nargin >= 2 && ~isdir(infolder)
    error('infolder must be a directory!');
elseif nargin >= 3 && ~isdir(outFolder)
    error('outFolder must be a directory!');
elseif nargin >= 4 && ~islogical(plotpassiveflag)
    error('plotpassiveflag must be either true (1) or false (0)!')
elseif nargin >= 5 && ~ischar(groupmode)
    error('groupmode must be either ''cell_actVhold'' or ''cell''!')
elseif nargin >= 5 && ischar(groupmode) && (~strcmp(groupmode, 'cell_actVhold') || ~strcmp(groupmode, 'cell'))
    error('groupmode must be either ''cell_actVhold'' or ''cell''!')
end

% Set defaults for optional arguments
if nargin < 2
    infolder = fullfile(homeDirectory, '/data_dclamp/take4/');
end
if nargin < 3
    if debugflag
        outFolder = fullfile(homeDirectory, '/data_dclamp/take4/debug/');
    else
        outFolder = fullfile(homeDirectory, '/data_dclamp/take4/');
    end
end
if nargin < 4
    plotpassiveflag = 1;
end
if nargin < 5
    groupmode = 'cell';
end

% Check if output directories exist
check_dir(outFolder);

%% Preparation
% Construct the full path to the data directory and other files
dataDirectory = fullfile(infolder, dataDirectoryName);
initialSlopesPath = fullfile(infolder, initialSlopesFile);
sweepInfoPath = fullfile(infolder, sweepInfoFile);

% Set suffix and title modification according to fitMode
[suffix, titleMod] = m3ha_specs_for_fitmode(fitMode);
fprintf('Using fit mode == %d ... \n', fitMode);

% Set logheaderSwpInfo && logvariablesSwpInfo depending on groupmode
if strcmp(groupmode, 'cell_actVhold')
    logheaderSwpInfo = logheaderSwpInfoCellActVHold;
    logvariablesSwpInfo = logvariablesSwpInfoCellActVHold;
elseif strcmp(groupmode, 'cell')
    logheaderSwpInfo = logheaderSwpInfoCell;
    logvariablesSwpInfo = logvariablesSwpInfoCell;
end

% Read in sweep info
sweepInfo = matfile(sweepInfoPath);
fprintf('Using sweep info matfile == %s ... \n', sweepInfoPath);

% Extract from sweep info
fnrow = sweepInfo.fnrow;
cellidrow = sweepInfo.cellidrow;
actVhold = sweepInfo.actVhold;
actIhold = sweepInfo.actIhold;

% Restrict data to fit if fitMode > 0
if fitMode > 0
    % Find indices of fnrow in dclampdatalog_take4.mat that will be used for fitting
    fnrow_old = sweepInfo.fnrow;            % file names for each sweep, used for indtofit
    cellidrow_old = sweepInfo.cellidrow;    % Cell ID # for each sweep, used for indtofit
    prow_old = sweepInfo.prow;              % Pharm condition for each sweep, used for indtofit
    grow_old = sweepInfo.grow;              % G incr for each sweep, used for indtofit
    indtofit = m3ha_find_ind_to_fit(fnrow_old, cellidrow_old, prow_old, ...
                                    grow_old, fitMode, infolder);

    % Restrict vectors
    fnrow = fnrow(indtofit);
    cellidrow = cellidrow(indtofit);
    actVhold = actVhold(indtofit);
    actIhold = actIhold(indtofit);
end

% Count the total number of sweeps used
nTotalSweeps = numel(fnrow);

% Construct full paths to the data files
dataFilesAll = construct_fullpath(fnrow, 'Directory', dataDirectory);

%% Do the passive fitting
if strcmp(groupmode, 'cell')        % Group the sweeps by cellidrow
    % Get all the unique cell IDs recorded 
    cellId = transpose(unique(cellidrow));

    % Count the number of cells to analyze 
    nCells = length(cellId);

    % Preparation for each cell
    cellName = cell(nCells, 1);
    fileBase = cell(nCells, 1);
    dataFileNames = cell(nCells, 1);
    dataFilePaths = cell(nCells, 1);
    parfor iCell = 1:nCells
        % Get the current cell ID
        cellIdThis = cellId(iCell);

        % Find all the sweep indices with this cell ID
        indThis = find(cellidrow == cellIdThis);

        % Get the data file names for this cell
        dataNamesThis = fnrow(indThis);

        % Get the cell name from the first 7 characters of the first file name
        cellNameThis = dataNamesThis{1}(1:7);

        % Store the cell name for this cell
        cellName{iCell} = cellNameThis;

        % Construct output file base for this cell
        fileBase{iCell} = ['Cell', num2str(cellIdThis), '_', cellNameThis];

        % Store the data file names for this cell
        dataFileNames{iCell} = dataNamesThis;

        % Store the data file paths for this cell
        dataFilePaths{iCell} = dataFilesAll(indThis);
    end

    % Extract passive parameters
    passiveParams = cell(nCells, 1);
    fitResults = cell(nCells, 1);
    fitObject = cell(nCells, 1);
    goodnessOfFit = cell(nCells, 1);
    algorithmInfo = cell(nCells, 1);
    decision = cell(nCells, 1);
    allResults = cell(nCells, 1);
    parfor iCell = 1:nCells
%    for iCell = 1:nCells
        % Extract from cell arrays
        cellIdThis = cellId(iCell);
        fileBaseThis = fileBase{iCell};
        dataNamesThis = dataFileNames{iCell};
        dataPathsThis = dataFilePaths{iCell};

        % Print message
        fprintf(['The current cellID # is ', num2str(cellIdThis), '  ...\n']);

        % Start timer
        tic;

        % Count the total number of sweeps for this cell
        nSwps = numel(dataPathsThis);
        fprintf(['The total number of sweeps for this cell is ', ...
                    num2str(nSwps), ' ...\n']);

        % Find passive parameters
        if nSwps > 0
            % Load vectors from data matfiles
            %   restricted to the current pulse response window
            fprintf('LOADING .mat files for the cell %s ...\n', fileBaseThis);
            [~, parsedData] = ...
                m3ha_parse_mat(dataPathsThis, 'LoadWindow', cprWin);
            tvec0 = parsedData.tvec0;
            ivec0s = parsedData.ivec0s;
            vvec0s = parsedData.vvec0s;
            ivec1s = parsedData.ivec1s;

            % Fix current pulse response traces that may have 
            %   out-of-balance bridges if requested
            if correctDcStepsFlag
                vvec0s = ...
                    m3ha_correct_unbalanced_bridge(dataNamesThis, vvec0s, ...
                                                ivec0s, initialSlopesPath);
            end

            % Analyze passive parameters such as input resistance (MOhm)
            fprintf('ANALYZING passive parameters for %s ...\n', fileBaseThis);
            [passiveParams{iCell}, fitResults{iCell}, fitObject{iCell}, ...
                goodnessOfFit{iCell}, algorithmInfo{iCell}, ...
                decision{iCell}, allResults{iCell}] = ...
                find_passive_params (tvec0, ivec0s, vvec0s, ...
                                     'HoldCurrent', actIhold, ...
                                     'PulseWindow', cpWin, ...
                                     'PulseResponseWindow', cprWin, ...
                                     'PlotFlag', plotpassiveflag, ...
                                     'OutFolder', outFolder, ...
                                     'FileBase', fileBase{iCell}, ...
                                     'Ivec1s', ivec1s, ...
                                     'Suffix', suffix, 'TitleMod', titleMod);
        end

        % End timer and print extra line
        toc;
        fprintf('\n');
    end
elseif strcmp(groupmode, 'cell_actVhold')           % TODO: Fix everything 
% TODO TODO TODO
    % Group the sweeps by cellidrow and vrow
    maxCellNumber = max(cellidrow);                % total number of cells recorded 
                                % (might not have sweeps of interest)
    nvholds = length(VholdBC);                % total number of Vhold conditions 
                                % (might not have sweeps of interest)
    files_cv = cell(maxCellNumber, nvholds);            % stores files for each set of conditions
    ct = 0;                            % counts the set number in use
    for v = 1:nvholds
        if v == 1
            v_ind = find(actVhold > VholdBC(v) - 2.5);
        elseif v == 2
            v_ind = find(actVhold > VholdBC(v) - 2.5 & actVhold <= VholdBC(v) + 2.5);
        elseif v == 3
            v_ind = find(actVhold <= VholdBC(v) + 2.5);
        end
        for c = 1:maxCellNumber
            indThisCellNumber = find(cellidrow == c);
            cv_ind = intersect(indThisCellNumber, v_ind);
            if ~isempty(cv_ind)
                ct = ct + 1;
                dataFileNames{ct} = dataFilesAll(cv_ind)';
                celln_set(ct) = c;
                holdPotential(ct) = VholdBC(v);
            end
        end
    end
    numsets = ct;            % total number of sets used

    % Extract passive parameters
    fn_set = cell(1, numsets);    % stores file base for the set
    passiveParams = cell(1, numsets);    % stores passive parameters for the set
    params_L_F2_set = cell(1, numsets);    % stores passive parameters estimated 
                        %    from the long pulse coefficients of the pooled falling phase fits
    params_S_R2_set = cell(1, numsets);    % stores passive parameters estimated 
                        %    from the short pulse coefficients of the pooled rising phase fits
    params_L_F1_set = cell(1, numsets);    % stores passive parameters estimated 
                        %    from the long pulse coefficients of the averaged falling phase fits
    params_S_R1_set = cell(1, numsets);    % stores passive parameters estimated 
                        %    from the short pulse coefficients of the averaged rising phase fits
    cpa_set = cell(1, numsets);    % stores vector of current pulse amplitudes (pA) for each sweep in the set
    pw_set = cell(1, numsets);    % stores vector of pulse width (ms) for each sweep in the set
    dvrec_set = cell(1, numsets);    % stores vector of overall change in membrane potential recorded (mV) 
                    %        for each sweep in the set
    rmse_R_set = cell(1, numsets);    % stores a vector of root-mean-squared errors (mV) in the rising phase for each sweep
                    %        for each sweep in the set
    rmse_F_set = cell(1, numsets);    % stores a vector of root-mean-squared errors (mV) in the falling phase for each sweep
                    %        for each sweep in the set
    parfor iSet = 1:numsets
        tic;
        c = celln_set(iSet);                % current cellID #
        v_val = holdPotential(iSet);            % current Vhold value
        fprintf(['The current cellID # is ', num2str(c), ' with Vhold = ', num2str(v_val), ' mV ...\n']);

        fn_set{iSet} = ['Cell', num2str(c), '_v', num2str(v_val)];    % file base for this set
        dataPathsThis = dataFileNames{iSet};            % files for this set
        nSwps = numel(dataPathsThis);            % total number of sweeps for this set
        fprintf(['The total number of sweeps for this set is ', num2str(nSwps), ' ...\n']);
        if nSwps > 0
            fprintf(['LOADING .mat files for the set ', fn_set{iSet}, ' ...\n']);

            % Load vectors from matfiles
            [tvec0, ivec0s, vvec0s, ivec1s] = load_vectors(dataPathsThis, cprWin);

            % Analyze passive parameters such as input resistance (MOhm)
            fprintf('ANALYZING passive parameters for %s ...\n', fn_set{iSet});
            [passiveParams{iSet}, cpa_set{iSet}, pw_set{iSet}, dvrec_set{iSet}, rmse_R_set{iSet}, rmse_F_set{iSet}, ...
                params_L_F2_set{iSet}, params_S_R2_set{iSet}, params_L_F1_set{iSet}, params_S_R1_set{iSet}] = ...
                find_passive_params (tvec0, ivec0s, vvec0s, ...
                    cpWin, cprWin, cpMid, plotpassiveflag, outFolder, fn_set{iSet}, ivec1s, fitMode);

            [passiveParams{iSet}, fitResults{iSet}, fitObject{iSet}, ...
                goodnessOfFit{iSet}, algorithmInfo{iSet}, ...
                decision{iSet}, allResults{iSet}] = ...
                find_passive_params (tvec0, ivec0s, vvec0s, ...
                                     'PulseWindow', cpWin, ...
                                     'PulseResponseWindow', cprWin, ...
                                     'PlotFlag', plotpassiveflag, ...
                                     'OutFolder', outFolder, ...
                                     'FileBase', fn_set{iSet}, ...
                                     'Ivec1s', ivec1s, ...
                                     'Suffix', suffix, 'TitleMod', titleMod);
        end
        toc;
        fprintf('\n');
    end

    %% Convert cell arrays to two dimensions
    fn_cv = cell(maxCellNumber, nvholds);
    params_cv = cell(maxCellNumber, nvholds);
    params_L_F2_cv = cell(maxCellNumber, nvholds);
    params_S_R2_cv = cell(maxCellNumber, nvholds);
    params_L_F1_cv = cell(maxCellNumber, nvholds);
    params_S_R1_cv = cell(maxCellNumber, nvholds);
    cpa_cv = cell(maxCellNumber, nvholds);
    pw_cv = cell(maxCellNumber, nvholds);
    dvrec_cv = cell(maxCellNumber, nvholds);
    rmse_R_cv = cell(maxCellNumber, nvholds);
    rmse_F_cv = cell(maxCellNumber, nvholds);
    for iSet = 1:numsets
        c = celln_set(iSet);                % current cellID #
        v_val = holdPotential(iSet);            % current Vhold value
        vn = find(VholdBC == v_val);
        fn_cv{c, vn} = fn_set{iSet};
        params_cv{c, vn} = passiveParams{iSet};
        params_L_F2_cv{c, vn} = params_L_F2_set{iSet};
        params_S_R2_cv{c, vn} = params_S_R2_set{iSet};
        params_L_F1_cv{c, vn} = params_L_F1_set{iSet};
        params_S_R1_cv{c, vn} = params_S_R1_set{iSet};
        cpa_cv{c, vn} = cpa_set{iSet};
        pw_cv{c, vn} = pw_set{iSet};
        dvrec_cv{c, vn} = dvrec_set{iSet};
        rmse_R_cv{c, vn} = rmse_R_set{iSet};
        rmse_F_cv{c, vn} = rmse_F_set{iSet};
    end
end

%% Save outputs
% Modify file names with suffix
outSheetNameMod = replace(outSheetName, '.xlsx', [suffix, '.xlsx']);
outMatNameByCellsMod = replace(outMatNameByCells, '.mat', [suffix, '.mat']);
outMatNameBySetsMod = replace(outMatNameBySets, '.mat', [suffix, '.mat']);

% Combine into a single table and Print to Excel file
cellInfoTable = table(cellId, cellName, fileBase, ...
                        dataFileNames, decision, fitObject);
passiveParamsTable = struct2table([passiveParams{:}]);
fitResultsTable = struct2table([fitResults{:}]);
goodnessOfFitTable = struct2table([goodnessOfFit{:}]);
resultsTable = [cellInfoTable, passiveParamsTable, ...
                fitResultsTable, goodnessOfFitTable];
writetable(resultsTable, fullfile(outFolder, outSheetNameMod));
% algorithmInfoTable = struct2table([algorithmInfo{:}]);
% resultsTable = [cellInfoTable, passiveParamsTable, ...
%                 fitResultsTable, goodnessOfFitTable, algorithmInfoTable];

% Save variables for each set/cell
if strcmp(groupmode, 'cell')        % Group the sweeps by cellidrow
    save(fullfile(outFolder, outMatNameByCellsMod), ...
            'cellId', 'cellName', 'fileBase', ...
            'dataFileNames', 'dataFilePaths', ...
            'passiveParams', 'fitResults', 'fitObject', 'goodnessOfFit', ...
            'algorithmInfo', 'decision', 'allResults', ...
            'logheaderSwpInfo', 'logvariablesSwpInfo', ...
            'logheaderParams', 'logvariablesParams', ...
            '-v7.3');
elseif strcmp(groupmode, 'cell_actVhold')           % TODO: Fix everything 
    save(fullfile(outFolder, outMatNameBySetsMod), ...
            'cellId', 'VholdBC', 'cellName', 'fileBase', ...
            'dataFileNames', 'dataFilePaths', ...
            'passiveParams', 'fitResults', 'fitObject', 'goodnessOfFit', ...
            'algorithmInfo', 'decision', 'allResults', ...
            'logheaderSwpInfo', 'logvariablesSwpInfo', ...
            'logheaderParams', 'logvariablesParams', ...
            '-v7.3');
end

% Resave variables for each sweep
resave_variables_for_sweep(outFolder, ...
            logheaderSwpInfo, logvariablesSwpInfo, ...
            logheaderParams, logvariablesParams, ...
            suffix, nTotalSweeps, fnrow, ...
            dataFileNames, cellId, passiveParams, allResults);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resave_variables_for_sweep (outFolder, ...
            logheaderSwpInfo, logvariablesSwpInfo, ...
            logheaderParams, logvariablesParams, ...
            suffix, nTotalSweeps, dataFileName, ...
            dataFileNamesSet, cellIdSet, passiveParams, allResults)
% Resave fitted variables for each sweep from grouped data for each set/cell
%   Note: holdPotential is empty for groupmode == 'cell'

% Construct output file name
outputfilename = ['dclampPassiveLog_bySwps', suffix, '.mat'];

% Count the number of sets
nSets = numel(dataFileNamesSet);

% Count the number of sweeps in each set
nSwpsEachSet = cellfun(@numel, dataFileNamesSet);

% Check if total number of sweeps is correct
nTotalSweeps2 = sum(nSwpsEachSet);
if nTotalSweeps2 ~= nTotalSweeps
    error('Total number of sweeps is incorrect!');
end

% Initialize variables
cellId = zeros(1, nTotalSweeps);         % cell number
setNumber = zeros(1, nTotalSweeps);      % set number
holdPotential = zeros(1, nTotalSweeps);  % rough holding potential (mV)
pulseAmplitude = zeros(1, nTotalSweeps); % current pulse amplitudes (pA)
pulseWidth = zeros(1, nTotalSweeps);     % pulse width (ms)
voltageChange = zeros(1, nTotalSweeps);  % change in membrane potential (mV)
paramNames = fieldnames(passiveParams{1});
nParams = numel(paramNames);
for iParam = 1:nParams
    eval([paramNames{iParam}, ' = zeros(1, nTotalSweeps);']);
end

% Find information for each sweep
for iSet = 1:nSets
    % Get the number of sweeps in this set
    nSwpsThis = nSwpsEachSet(iSet);

    % Get the file names for this set
    dataFileNamesThis = dataFileNamesSet{iSet};

    % Get the cell ID for this set
    cellIdThis = cellIdSet(iSet);

    % Get the parameters for this set
    passiveParamsThis = passiveParams{iSet};

    % Get the other results for this set
    allResultsThis = allResults{iSet};
    pulseAmplitudeThis = allResultsThis.pulseAmplitude;
    pulseWidthThis = allResultsThis.pulseWidth;
    voltageChangeThis = allResultsThis.voltageChange;


    for iSwp = 1:nSwpsThis
        % Get the current file name
        fileNameThis = dataFileNamesThis{iSwp};

        % Find the sweep index to save
        idxThis = find_in_strings(fileNameThis, dataFileName);

        % Save info for this sweep
        cellId(idxThis) = cellIdThis;
        setNumber(idxThis) = iSet;
        pulseAmplitude(idxThis) = pulseAmplitudeThis(iSwp);
        pulseWidth(idxThis) = pulseWidthThis(iSwp);
        voltageChange(idxThis) = voltageChangeThis(iSwp);
        if isstruct(passiveParamsThis)
            for iParam = 1:nParams
                eval([paramNames{iParam}, ...
                        '(idxThis) = passiveParamsThis.(paramNames{iParam});']);
            end
        else
            for iParam = 1:nParams
                eval([paramNames{iParam}, '(idxThis) = NaN;']);
            end
        end
    end
end

passivedatafn = fullfile(outFolder, outputfilename);
command = ['save(passivedatafn, ', ...
    '''logheaderSwpInfo'', ''logvariablesSwpInfo'', ', ...
    '''logheaderParams'', ''logvariablesParams'', ', ...
    '''dataFileName'', ''cellId'', ''setNumber'', ', ...
    '''pulseAmplitude'', ''pulseWidth'', ''voltageChange'', ']; 
for iParam = 1:nParams
    command = [command, sprintf('''%s'', ', paramNames{iParam})];
end
command = [command, '''-v7.3'');'];
eval(command);

% if ~isempty(holdPotential)
%     VholdBCrow(ind) = holdPotential(iSet);
% end
%rmse_R_row = zeros(1, nTotalSweeps);     % root-mean-squared error (mV) in the rising phase
%rmse_F_row = zeros(1, nTotalSweeps);     % root-mean-squared error (mV) in the falling phase
% rmse_R_row(idxThis) = rmse_R_set{iSet}(iSwp);
% rmse_F_row(idxThis) = rmse_F_set{iSet}(iSwp);
%, ...
%    '''rmse_R_row'', ''rmse_F_row'', '];
% if ~isempty(holdPotential)
%     command = [command, '''holdPotential'', '];
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% OLD CODE

numswps = length(fnrow);                % total number of sweeps of interest

%% Set suffices for each fitMode
[suffix, ~] = SpecsForFitmode (fitMode);

%% Create folders for saving files
outfolder_pass = fullfile(outFolder, ['/passive', suffix, '/']);
if exist(outfolder_pass, 'dir') ~= 7
    mkdir(outfolder_pass);
    fprintf('Made directory %s \n', outfolder_pass);
end

% for iSet = 32*3 - 2
% for iSet = 17*3
% for iSet = 25*3 - 2
% for iSet = 1
% for iSet = 3

    c = ceil(iSet/nvholds);                % current cellID #
    v = mod(iSet - 1, nvholds) + 1;            % current vhold #

    if (c - 1) * nvholds + v ~= iSet
        error('code is incorrect!');
    end

    v_val = Vhold(v);                % current Vhold value
    fprintf(['The current Vhold is ', num2str(v_val), ' mV ...\n']);

maxsets = maxCellNumber * nvholds;                % total number of sets (might not have sweeps of interest)

Vhold = [-60; -65; -70];    % Possible Vhold values (LJP-corrected)
vrow = m.vrow;

cprWin = [95, 500];            % Window in which the current pulse response would lie (ms)

if exist('/media/adamX/m3ha/', 'dir') == 7
    homeDirectory = '/media/adamX/m3ha/';
elseif exist('/scratch/al4ng/m3ha/', 'dir') == 7
    homeDirectory = '/scratch/al4ng/m3ha/';
else
    error('Valid homeDirectory does not exist!');
end

[params_cell{iCell}, cpa_cell{iCell}, pw_cell{iCell}, dvrec_cell{iCell}, rmse_R_cell{iCell}, rmse_F_cell{iCell}, ...
    params_L_F2_cell{iCell}, params_S_R2_cell{iCell}, params_L_F1_cell{iCell}, params_S_R1_cell{iCell}] = ...
    find_passive_params (tvec0, ivec0s, vvec0s, ...
        cpWin, cprWin, cpMid, plotpassiveflag, outFolder, fileBase{iCell}, ivec1s, fitMode);

if preallocateflag
    filesCell = cell(1, nCells);            % stores files for each cell
    cellId = zeros(1, nCells);        % stores cell ID number for each cell
end

% Get the maximum cell number recorded 
maxCellNumber = max(cellidrow);

% Count the number of cells to use
ct = 0;
for iCellNumber = 1:maxCellNumber
    % Find all the sweep indices with this cell number
    indThisCellNumber = find(cellidrow == iCellNumber);

    % If there is a sweep, include this cell
    if ~isempty(indThisCellNumber)
        % Increment cell count
        ct = ct + 1;

        % Store the .mat file names for this cell
        filesCell{ct} = transpose(fnrow(indThisCellNumber));

        % Store the cell number for this cell
        cellId(ct) = iCellNumber;
    end
end
nCells = ct;

params_cell = cell(1, nCells);    % stores passive parameters for the cell
params_L_F2_cell = cell(1, nCells);    % stores passive parameters estimated 
                    %    from the long pulse coefficients of the pooled falling phase fits
params_S_R2_cell = cell(1, nCells);    % stores passive parameters estimated 
                    %    from the short pulse coefficients of the pooled rising phase fits
params_L_F1_cell = cell(1, nCells);    % stores passive parameters estimated 
                    %    from the long pulse coefficients of the averaged falling phase fits
params_S_R1_cell = cell(1, nCells);    % stores passive parameters estimated 
                    %    from the short pulse coefficients of the averaged rising phase fits
cpa_cell = cell(1, nCells);    % stores a vector of current pulse amplitudes (pA) for each sweep in the cell
pw_cell = cell(1, nCells);    % stores a vector of pulse widths (ms) for each sweep in the cell
dvrec_cell = cell(1, nCells);    % stores a vector of overall changes in membrane potential recorded (mV) 
                %        for each sweep in the cell
rmse_R_cell = cell(1, nCells);% stores a vector of root-mean-squared errors (mV) in the rising phase for each sweep
                %        for each sweep in the cell
rmse_F_cell = cell(1, nCells);% stores a vector of root-mean-squared errors (mV) in the falling phase for each sweep
                %        for each sweep in the cell

% Preallocate according to the fitMode and groupmode
%%% NOT FINISHED
if preallocateflag
    nCells = 49;
    if fitMode == 0
        numsets = 3*49;            % Total number of cell-Vhold conditions recorded
    elseif fitMode == 1
        numsets = 88;
    elseif fitMode == 2
        numsets = 88;
    end
end

preallocateflag = 0;        % Only use this after total number of sets are certain
if preallocateflag
    dataFileNamesSet = cell(1, numsets);            % stores files for each set of conditions, linear form
    celln_set = zeros(1, numsets);            % stores cell ID number for each set of conditions
    holdPotential = zeros(1, numsets);        % stores Vhold for each set of conditions
end

cpMid = 105;                % approximate midpoint of the current pulse (ms)

[tvec0, ivec0s, vvec0s, ivec1s] = load_vectors(dataPathsThis, cprWin);

function [tvec0, ivec0s, vvec0s, ivec1s] = load_vectors(dataFiles, cprWin);
% Load vectors from matfiles

% Obtain full time vector from .mat data of first sweep
firstfile = dataFiles{1};
if exist(firstfile, 'file') ~= 2
    error(['This mat file: ', firstfile, ' is missing!!']);
end
m = matfile(firstfile);
tvec0 = m.d_orig(:, 1);

% Extract info from time vector
sims = tvec0(2) - tvec0(1);            % sampling interval in ms
cprwin_begin = find(tvec0 >= cprWin(1), 1);    % first index of interest
cprwin_end = find(tvec0 <= cprWin(2), 1, 'last');    % final index of interest
cprwin_ind = cprwin_begin:cprwin_end;        % indices of interest
ndps = length(cprwin_ind);            % total number of data points of interest

% Restrict time vector
tvec0 = tvec0(cprwin_ind);

% Obtain restricted vectors from .mat data of all sweeps
nSwps = numel(dataFiles);            % total number of sweeps for this set
ivec0s = zeros(ndps, nSwps);
vvec0s = zeros(ndps, nSwps);
ivec1s = zeros(ndps, nSwps);
parfor iSwp = 1:nSwps        % FOR each sweep
    thisfile = dataFiles{iSwp};
    if exist(thisfile, 'file') ~= 2
        error(['This mat file: ', thisfile, ' is missing!!']);
    end
    m = matfile(thisfile);
    fileinfo = sprintf('First file == %s\nThis file == %s\n', firstfile, thisfile);

    % Check if time vectors are the same
    tvec0_new = m.d_orig(:, 1);
    sims_new = tvec0_new(2) - tvec0_new(1);
    cprwin_begin_new = find(tvec0_new >= cprWin(1), 1);    % first index of interest
    cprwin_end_new = find(tvec0_new <= cprWin(2), 1, 'last');    % final index of interest
    cprwin_ind_new = cprwin_begin_new:cprwin_end_new;        % indices of interest
    ndps_new = length(cprwin_ind_new);
    if ndps_new ~= ndps
        error('number of data points not the same!\n%s', fileinfo);
    elseif sims_new ~= sims
        error('sampling interval not the same!\n%s', fileinfo);
    end

    % Extract data restricted to cprWin
    ivec0s(:, iSwp) = m.d_orig(cprwin_ind, 3);
    vvec0s(:, iSwp) = m.d_orig(cprwin_ind, 4);
    ivec1s(:, iSwp) = m.d_mf(cprwin_ind, 3);
end

%% Save variables for each set
save(fullfile(outFolder, outputfilename), ...
        'logheaderParams', 'logvariablesParams', ...
        'fn_set', 'dataFileNames', 'celln_set', 'holdPotential', ...
        'passiveParams', 'cpa_set', 'pw_set', 'dvrec_set', ...
        'rmse_R_set', 'rmse_F_set', ...
        'params_L_F2_set', 'params_S_R2_set', 'params_L_F1_set', 'params_S_R1_set', ...
        'fn_cv', 'files_cv', 'params_cv', 'cpa_cv', 'pw_cv', 'dvrec_cv', ...
        'rmse_R_cv', 'rmse_F_cv', ...
        'params_L_F2_cv', 'params_S_R2_cv', 'params_L_F1_cv', 'params_S_R1_cv', ...
        '-v7.3');

%% Resave variables for each sweep
resave_variables_for_sweep(outFolder, logheaderSwpInfo, logvariablesSwpInfo, ...
            logheaderParams, logvariablesParams, ...
            suffix, nTotalSweeps, fnrow, numsets, dataFileNames, celln_set, ...
            cpa_set, pw_set, dvrec_set, rmse_R_set, rmse_F_set, passiveParams, holdPotential)

if strcmp(groupmode, 'cell')        % Group the sweeps by cellidrow
    nSets = nCells;
elseif strcmp(groupmode, 'cell_actVhold')           % TODO: Fix everything 
    nSets = numsets;
end

if ~isempty(holdPotential)
end

    if debugflag
        infolder = fullfile(homeDirectory, '/data_dclamp/take4/debug/');
    else
    end

outputfilename = ['dclampPassiveLog_byCells', suffix, '.mat'];
outputfilename = ['dclampPassiveLog_bySets', suffix, '.mat'];

%}
