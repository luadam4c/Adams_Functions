function [done, outparams, hfig] = m3ha_optimizer_4compgabab (outparams, hfig)
% OPTIMIZER  Passes parameters to NEURON (which runs simulations and saves
% data as .dat files), loads .dat file data, plots data and compares with
% experimental data. Save all this as well (timestamp for now)
%
% realData  Experimental data from OPTIMIZERGUI.m
% 
% outparams   GUI parameters from OPTIMIZERGUI.m
%   outparams.swpWeights : weights of sweeps
%   outparams.swpedges : left and right edges of sweep region to analyze
%   outparams.ltsWeights : weights of LTS parameters
%   outparams.modeselected : mode to run OPTIMIZER.m in (either
%       'modebutton_auto' or 'modebutton_manual')
%   outparams.neuronparams : values of parameters for NEURON 
%   outparams.neuronparamnames : names of parameters for NEURON (e.g. pcabar1,
%       actshift_itGHK, etc)
%   outparams.simplexParams : values of simplex parameters for
%   FMINSEARCH2.m
%   outparams.simplexParamNames : names of simplex parameters for
%   FMINSEARCH2.m
%   outparams.sortedswpnum : sweep numbers of sweeps called into
%   OPTIMIZER.m; use to figure out amplitude of current injected during CCIV
%   protocol
% 
% done (0 or 1)   Output to pass to OPTIMIZERGUI to signal end of
%                     NEURON and FMINSEARCH2 (if in auto mode)
%
% Requires:
%       cd/check_dir.m
%       cd/find_in_strings.m
%       cd/isemptystruct.m
%       cd/locate_functionsdir.m
%       cd/m3ha_compare_neuronparams2.m
%       cd/m3ha_fminsearch3.m
%       cd/m3ha_log_errors_params.m
%       cd/m3ha_neuron_create_new_initial_params.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/set_fields_zero.m
%       cd/restore_fields.m
%       cd/structs2vecs.m
%       ~/Downloaded_Functions/subplotsqueeze.m
%
% Used by:
%       cd/singleneuronfitting42.m and later versions
%       cd/m3ha_optimizergui_4compgabab.m
% 
% By Christine Lee Kyuyoung 2011-01-16
% last modified by CLK 2014-04
% 2016-07-15 - Added MANUAL WITH JITTER
% 2016-07-19 - Changed cols{k} to cols(cgn, :)
% 2016-07-19 - Changed order of subplot(5,1,1) & subplot(5,1,2) in simtraces
% 2016-07-20 - Modified JITTER mode and added AUTO WITH JITTER
% 2016-07-20 - Changed npercg & ncg, simtraces axes
% 2016-10-03 - Added current pulse response
% 2016-10-06 - Reorganized code
% 2016-10-06 - Moved update_sliderposition to optimizergui_4compgabab.m
% 2016-10-06 - xlimits now uses outparams.fitreg instead of outparams.swpedges
% 2016-10-06 - Renamed lots of figures
% 2016-10-06 - Renamed figure handles so that they are now all in a structure 
%               hfig that is passed to and from functions
% 2016-10-07 - current pulse electrode current and responses are now graphed
% 2016-10-07 - Plot flags are now set to zero before optimization and 
%               set back afterwards
% 2017-01-17 - Modified runauto so that it will fit current pulse response
% 2017-01-25 - Randomize initial parameters and run simplex nInitConds times;
%               added generate_IC() & find_best_params()
% 2017-01-26 - Completed generate_IC()
% 2017-04-21 - Now uses m3ha_log_errors_params.m to generate 
%               errors_and_params_log_manual.csv
% 2017-04-21 - Added outFolderName to output file names
% 2017-05-02 - Now uses parfor to run through all initial conditions
% 2017-05-02 - Added 'Iteration #' to compare_params.csv
% 2017-05-02 - Changed simplex count to outparams.simplexNum and 
%               fixed the fact that it wasn't updated during active fitting
% 2017-05-12 - Now saves runauto_results in outfolder
% 2017-05-13 - Now updates outparams.runnumTotal here
% 2017-05-13 - Now returns errCpr & err in outparams
% 2017-05-13 - Now gets outFolderName from outparams
% 2017-05-13 - Now updates outparams.neuronparams to 
%               initCond.neuronparams in generate_IC()
% 2017-05-13 - Now initializes outparams0 within parfor loop
% 2017-05-15 - Now changes outparams.prefix so that '_cpr' is already 
%               incorporated before passive fitting
% 2017-05-15 - Now changes neuronparams_use before calling fminsearch3_4compgabab 
%               to reflect passive or active
% 2017-05-16 - Now runs NEURON with full plots both before and after fitting
% 2017-05-16 - parfor is now conditional on outparams.MaxNumWorkersIC
% 2017-05-17 - update_errorhistoryplot() now plots lts errors if computed
% 2017-05-17 - Added outparams.fitPassiveFlag & outparams.fitActiveFlag
% 2017-05-17 - Moved update_sweeps_figures() to m3ha_neuron_run_and_analyze.m
% 2017-05-19 - Fixed the fact that simplex.initError was used to select  
%               best simplex; simplex.error is now simplex.totalError
% 2017-05-22 - Changed line width and indentation
% 2017-05-22 - Added outparams.fitTogetherFlag
% 2017-05-23 - Removed modeselected from outparams and 
%               replaced with updated outparams.runMode
% 2017-05-23 - Added otherwise to all switch statements
% 2017-06-19 - Fixed runmode under runmanual
% 2017-07-29 - Now saves the best parameters under /pfiles/ if err improves
% 2017-08-10 - Now plots activation/inactivation & I-V curves after optimization
% 2018-01-24 - Added isdeployed
% 2018-03-08 - Changed compare_neuronparams() to compare_neuronparams2()
% 2018-04-24 - Now does not use parfor to run through all initial conditions
%               Some initial conditions take much longer than others
% 2018-08-08 - Now makes sure plotOverlappedFlag is false during fitting
% 2018-08-08 - Now forces constrains all initial conditions for 
%               epas and gpas according to epasEstimate and RinEstimate
% 2018-08-10 - Changed outparams.fitregCpr to use outparams.fitwinCpr
% 2018-08-16 - Now does not use parfor for active fitting
% 2018-12-10 - Updated placement of jitterFlag
% 2018-12-11 - Moved code to m3ha_neuron_create_new_initial_params.m
% 2019-11-13 - Now copies final parameters without renaming
% 2019-11-15 - Now passes in outparams for varargin 
%                   for m3ha_neuron_run_and_analyze.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% If not compiled, add directories to search path for required functions
if ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for check_dir.m, etc.
    addpath(fullfile(functionsDirectory, 'Adams_Functions')); 

    % Add path for subplotsqueeze.m, etc.
    addpath(fullfile(functionsDirectory, 'Downloaded_Functions')); 
end

%% Preparation
% Change fitting flags if necessary
% TODO: Fix
%{
if outparams.fitPassiveFlag && ...
    ~any(neuronparamispas .* outparams.neuronparams_use)
    outparams.fitPassiveFlag = false;
    fprintf(['No passive parameter is fitted, ', ...
            'so passive fitting is turned OFF!\n\n']);
end
if outparams.fitActiveFlag && ...
    ~any(~neuronparamispas .* outparams.neuronparams_use)
    outparams.fitActiveFlag = false;
    fprintf(['No active parameter is fitted, ', ...
            'so active fitting is turned OFF!\n\n']);
end
%}

% Initialize the prefix as the date-time-cell stamp
outparams.prefix = outparams.dateTimeCellStamp;

% If in jitter mode and if parameter is checked, add jitter
if outparams.runMode == 3
    outparams.jitterFlag = true;
else
    outparams.jitterFlag = false;
end

%% Run based on manual or auto mode
switch outparams.runMode
case {1, 3}
    [~, ~, outparams, hfig] = runmanual(outparams, hfig);
case 2
    %% Save old NEURON parameters
    oldNeuronParams = outparams.neuronParamsTable.Value;

    %% Run NEURON with baseline parameters with full plots
    prefixOrig = outparams.prefix;                         % save original prefix
    outparams.prefix = [outparams.prefix, '_bef'];          % change prefix for "before fitting"
    [~, ~, outparams, hfig] = runmanual(outparams, hfig);
    outparams.prefix = prefixOrig;                         % restore original prefix
    drawnow

    %% Optimize parameters
    [~, ~, outparams, hfig] = runauto(outparams, hfig);
    drawnow

    %% Run NEURON with best parameters again with full plots
    prefixOrig = outparams.prefix;                         % save original prefix
    outparams.prefix = [outparams.prefix, '_aft'];          % change prefix for "after fitting"
    [~, ~, outparams, hfig] = runmanual(outparams, hfig);
    outparams.prefix = prefixOrig;                         % restore original prefix
    drawnow

    %% If error improved, copy final parameters to bestParamsDirectory
    if ~outparams.debugFlag && ...
            outparams.err{outparams.runnumTotal}.totalError < ...
                outparams.err{1}.totalError       
        % Get the final parameters file name for this cell
        finalParamsFile = fullfile(outparams.outFolder, ...
                            [outparams.prefix, '_aft.p']);

        % Copy to bestParamsDirectory
        copyfile(finalParamsFile, outparams.bestParamsDirectory);
        fprintf('Best parameters copied for %s!!\n\n', outparams.prefix);
    end

    %% Compare NEURON parameters before and after 
    %   by plotting activation/inactivation & I-V curves
    newNeuronParams = outparams.neuronParamsTable.Value;
    neuronParamNames = outparams.neuronParamsTable.Properties.RowNames;
    bothParamValues = {oldNeuronParams, newNeuronParams};
    bothParamNames = {neuronParamNames, neuronParamNames};
    suffices = {['_', outparams.prefix, '_bef'], ...
                ['_', outparams.prefix, '_aft']};
% TODO: Fix
%     m3ha_compare_neuronparams2(bothParamValues, bothParamNames, suffices, ...
%                         'OutFolder', outparams.outFolder);

case 4                                       %%% TODO: Unfinished
    [~, ~, outparams, hfig] = ...
        runauto_w_jitter(outparams, hfig);
otherwise
    outparams.runMode = 1;
    fprintf('Warning: run mode out of range, changed to 1!\n\n');
end

%% Make all figures visible and update
if outparams.showSweepsFlag
    figs = fieldnames(hfig);
    for k = 1:numel(figs)
        set(hfig.(figs{k}), 'Visible', 'on');
        drawnow
    end
end

%% Used for runbutton_toggle in m3ha_optimizergui_4compgabab.m
done = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errCpr, err, outparams, hfig] = runmanual(outparams, hfig)
%% MANUAL or JITTER modes: Simulate each sweep once and compare with real data

% Update run counts
switch outparams.runMode
case {1, 2}         % manual mode is also run before and after auto mode
    outparams.runnumManual = outparams.runnumManual + 1;
case 3
    outparams.runnumJitter = outparams.runnumJitter + 1;
otherwise
    error('runmode should be 1 or 3!\n');
end
outparams.runnumTotal = outparams.runnumTotal + 1;

% Simulate current pulse response only if outparams.fitPassiveFlag is true
if outparams.fitPassiveFlag         % if fitting passive parameters
    %%%%%% 
    %%%%%%%%%%%%%
    % Set sim mode
    outparams.simMode = 'passive';

    % Save original prefix
    prefixOrig = outparams.prefix;

    % Change prefix for passive-parameters-fitting
    outparams.prefix = [outparams.prefix, '_cpr'];

    % Set grouping vector to be vhold
    outparams.grouping = outparams.vhold;

    % Run all simulations once
    [errCpr, outparams.hfig] = ...
        m3ha_neuron_run_and_analyze(outparams.neuronParamsTable, outparams);

    % Change prefix back
    outparams.prefix = prefixOrig;
    %%%%%%%%%%%%%
    %%%%%% 

    % Log errors and parameters and save parameters as .p file
    if ~isemptystruct(errCpr)
        m3ha_log_errors_params(outparams.logFileName, outparams, errCpr);
    end
else
    if outparams.runnumTotal > 1    % if this is not the first run
        % Use errCpr from last run
        errCpr = outparams.errCpr{outparams.runnumTotal-1};   
    else                            % if this is the first run
        % Use empty structure
        errCpr = struct;
    end
end
outparams.errCpr{outparams.runnumTotal} = errCpr;         % store error structure in outparams

% Simulate GABAB IPSC response only if outparams.fitActiveFlag is true
if outparams.fitActiveFlag
    %%%%%% 
    %%%%%%%%%%%%%
    % Set sim mode
    outparams.simMode = 'active';

    % Run all simulations once
    [err, outparams.hfig] = ...
        m3ha_neuron_run_and_analyze(outparams.neuronParamsTable, outparams);              
    %%%%%%%%%%%%%
    %%%%%% 

    % Log errors and parameters and save parameters
    if ~isemptystruct(err)
        m3ha_log_errors_params(outparams.logFileName, outparams, err);
    end
else
    if outparams.runnumTotal > 1    % if this is not the first run
        % Use errCpr from last run
        err = outparams.err{outparams.runnumTotal-1};   
    else                            % if this is the first run
        % Use empty structure
        err = struct;
    end
end
outparams.err{outparams.runnumTotal} = err;                 % store error structure in outparams

% Update error history plot
% TODO: Fix
% update_errorhistoryplot(hfig, outparams);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errCpr, err, outparams, hfig] = runauto(outparams, hfig)
%% AUTO mode: Find best parameters with m3ha_fminsearch3.m

% Update run counts
outparams.runnumAuto = outparams.runnumAuto + 1;
outparams.runnumTotal = outparams.runnumTotal + 1;

% Extract constants from outparams
cellName = outparams.cellName;
oldOutFolderName = outparams.outFolder;

% Extract number of initial conditions
idxNInitConds = find_in_strings('nInitConds', outparams.autoParamNames);
nInitConds = outparams.autoParams(idxNInitConds);

% Turn off all flags for stats and plots for fminsearch
[outparams] = ...
    set_fields_zero(outparams, ...
        'ltsBurstStatsFlag', 'saveLtsInfoFlag', 'saveLtsStatsFlag', ...
        'saveSimCmdsFlag', 'saveStdOutFlag', 'saveSimOutFlag', ...
        'plotIndividualFlag', 'plotOverlappedFlag', 'plotResidualsFlag', ...
        'plotConductanceFlag', 'plotCurrentFlag', ...
        'plotIpeakFlag', 'plotLtsFlag', 'plotStatisticsFlag', ...
        'plotSwpWeightsFlag');

% Fit parameters
if outparams.fitTogetherFlag        % passive and active fitting done together from single initCond
    % Prepare for fitting
    fprintf('Fitting all parameters for cell %s ... \n', cellName);
    tstartAllFit = tic();          % tracks time for all-parameters-fitting
    fprintf('\n');

    % Fit all parameters
    cprInitialConditionsAll = cell(1, nInitConds);            % stores the initial parameters for each passive simplex run
    cprSimplexOutAll = cell(1, nInitConds);    % stores the simplex outputs for each passive simplex run
    cprExitFlagAll = zeros(1, nInitConds);     % stores the exitflags for each passive simplex run
    initialConditionsAll = cell(1, nInitConds);                % stores the initial parameters for each active simplex run
    simplexOutAll = cell(1, nInitConds);        % stores the simplex outputs for each active simplex run
    exitFlagAll = zeros(1, nInitConds);         % stores the exitflags for each active simplex run
    parfor iInitCond = 1:nInitConds
%    for iInitCond = 1:nInitConds
        % Initialize outparams0 for parfor
        outparams0 = outparams;

        % Store the old simplex count
        oldSimplexCt = outparams.simplexNum;

        % Prepare outparams0 for passive fit
        [outparams0, prefixOrig, neuronParamsUseOrig] = ...
            prepare_outparams_passive(outparams0);

        % Prepare outparams0 for simplex
        [outparams0] = ...
            prepare_outparams_simplex(outparams0, oldOutFolderName, oldSimplexCt, iInitCond);

        % Generate initial parameters for this run
        [cprInitialConditionsAll{iInitCond}, outparams0] = ...
            generate_IC(outparams0, iInitCond);

        % Run fminsearch for current pulse response
        %%%%%%
        %%%%%%%%%%%%%
        [cprSimplexOutAll{iInitCond}, cprExitFlagAll(iInitCond)] = ...
            m3ha_fminsearch3(outparams0);
        %%%%%%%%%%%%%
        %%%%%%

        % Restore outparams0 after passive fitting
        [outparams0] = ...
            restore_outparams_passive(outparams0, prefixOrig, neuronParamsUseOrig);

        % Use the optimized parameters
        outparams0.neuronParamsTable = cprSimplexOutAll{iInitCond}.neuronParamsTable;

        % Prepare outparams0 for active fit
        [outparams0, neuronParamsUseOrig] = ...
            prepare_outparams_active(outparams0);

        % Prepare outparams0 for simplex
        [outparams0] = ...
            prepare_outparams_simplex(outparams0, oldOutFolderName, ...
                                        oldSimplexCt + nInitConds, iInitCond);

        % Generate initial parameters for this run
        [initialConditionsAll{iInitCond}, outparams0] = ...
            generate_IC(outparams0, iInitCond);

        % Run fminsearch for GABAB IPSC response
        %%%%%%
        %%%%%%%%%%%%%
        [simplexOutAll{iInitCond}, exitFlagAll(iInitCond)] = ...
            m3ha_fminsearch3(realDataIpscr, outparams0);
        %%%%%%%%%%%%%
        %%%%%%

        % Restore outparams0 after active fit
        outparams0 = restore_outparams_active(outparams0, neuronParamsUseOrig);
    end

    % Update # of simplex runs and # of simplex steps
    outparams.simplexNum = outparams.simplexNum + nInitConds * 2;
    outparams.simplexIterCount = ...
        outparams.simplexIterCount + ...
        sum(cellfun(@(x) x.ctIterations, cprSimplexOutAll)) + ...
        sum(cellfun(@(x) x.ctIterations, simplexOutAll));

    % Find best of the optimized parameters
    cprCompareparamsfile = fullfile(oldOutFolderName, ...
                                    [outparams.prefix, '_cpr_compare_params_IC_', ...
                                    num2str(outparams.simplexNum - nInitConds * 2 + 1), ...
                                    'to', num2str(outparams.simplexNum - nInitConds), '.csv']);
    [cprSimplexOutBest, errCpr] = ...
        find_best_params(cprSimplexOutAll, cprInitialConditionsAll, cprCompareparamsfile);
    compareparamsfile = fullfile(oldOutFolderName, ...
                                [outparams.prefix, '_compare_params_IC_', ...
                                num2str(outparams.simplexNum - nInitConds + 1), ...
                                'to', num2str(outparams.simplexNum), '.csv']);
    [simplexOutBest, err] = ...
        find_best_params(simplexOutAll, initialConditionsAll, compareparamsfile);

    % Update outparams.neuronParamsTable to the best of the optimized parameters
    outparams.neuronParamsTable = simplexOutBest.neuronParamsTable;

    % Save outputs in matfile
    if outparams.saveMatFileFlag
        save(fullfile(oldOutFolderName, ...
                        [outparams.prefix, '_runauto_results_IC_', ...
                        num2str(outparams.simplexNum - nInitConds + 1), ...
                        'to', num2str(outparams.simplexNum), '.mat']), ...
            'cprInitialConditionsAll', 'cprSimplexOutAll', 'cprExitFlagAll', ...
            'initialConditionsAll', 'simplexOutAll', 'exitFlagAll', ...
            'outparams', '-v7.3');
    end


    % Print time elapsed
    time_taken_allfit = toc(tstartAllFit);
    fprintf('It took %g seconds to run the simplex method on %d initial conditions!!\n', ...
            time_taken_allfit, nInitConds);
    fprintf('\n');

else                                % do passive fitting, find best params, then do active fitting

    % Fit passive parameters
    if outparams.fitPassiveFlag && ~outparams.fitTogetherFlag
        % Prepare for passive-parameters-fitting
        fprintf('Fitting passive parameters for cell %s ... \n', cellName);
        tStartPassiveFit = tic();      % tracks time for passive-parameters-fitting
        fprintf('\n');

        % Prepare outparams for passive fit
        [outparams, prefixOrig, neuronParamsUseOrig] = ...
            prepare_outparams_passive(outparams);
        
        % Optimize passive parameters by fitting to current pulse response
        cprInitialConditionsAll = cell(1, nInitConds);            % stores the initial parameters for each simplex run
        cprSimplexOutAll = cell(1, nInitConds);    % stores the simplex outputs for each simplex run
        cprExitFlagAll = zeros(1, nInitConds);     % stores the exitflags for each simplex run
%        parfor iInitCond = 1:nInitConds
        for iInitCond = 1:nInitConds
            % Initialize outparams0 for parfor
            outparams0 = outparams;

            % Store the old simplex count
            oldSimplexCt = outparams.simplexNum;

            % Prepare outparams0 for simplex
            [outparams0] = ...
                prepare_outparams_simplex(outparams0, oldOutFolderName, oldSimplexCt, iInitCond);

            % Generate initial parameters for this run
            [cprInitialConditionsAll{iInitCond}, outparams0] = generate_IC(outparams0, iInitCond);

            % Run fminsearch for current pulse response
            %%%%%%
            %%%%%%%%%%%%%
            [cprSimplexOutAll{iInitCond}, cprExitFlagAll(iInitCond)] = ...
                m3ha_fminsearch3(outparams0);
            %%%%%%%%%%%%%
            %%%%%%
        end

        % Restore outparams after passive fit
        [outparams] = ...
            restore_outparams_passive(outparams, prefixOrig, neuronParamsUseOrig);

        % Update # of simplex runs and # of simplex steps
        outparams.simplexNum = outparams.simplexNum + nInitConds;    
        outparams.simplexIterCount = outparams.simplexIterCount + ...
                                        sum(cellfun(@(x) x.ctIterations, cprSimplexOutAll));

        % Find best of the optimized parameters
        cprCompareparamsfile = fullfile(oldOutFolderName, ...
                                        [outparams.prefix, '_compare_params_IC_', ...
                                        num2str(outparams.simplexNum - nInitConds + 1), ...
                                        'to', num2str(outparams.simplexNum), '.csv']);
        [cprSimplexOutBest, errCpr] = ...
            find_best_params(cprSimplexOutAll, cprInitialConditionsAll, cprCompareparamsfile);

        % Update outparams.neuronParamsTable to the best of the optimized parameters
        outparams.neuronParamsTable = cprSimplexOutBest.neuronParamsTable;

        % Save outputs in matfile
        if outparams.saveMatFileFlag
            save(fullfile(oldOutFolderName, ...
                            [outparams.prefix, '_runauto_results_IC_', ...
                            num2str(outparams.simplexNum - nInitConds + 1), ...
                            'to', num2str(outparams.simplexNum), '.mat']), ...
                'cprInitialConditionsAll', 'cprSimplexOutAll', ...
                'cprExitFlagAll', 'cprSimplexOutBest', 'errCpr', ...
                'outparams', '-v7.3');
        end

        % Print time elapsed
        timeTakenPassiveFit = toc(tStartPassiveFit);
        fprintf('It took %g seconds to run the simplex method on %d initial conditions!!\n', ...
                timeTakenPassiveFit, nInitConds);
        fprintf('\n');
    else
        if outparams.runnumTotal > 1    % if this is not the first run
            errCpr = outparams.errCpr{outparams.runnumTotal-1};       % use errCpr from last run    
        else                            % if this is the first run
            errCpr = struct;                                           % empty structure
        end
    end

    % Fit active parameters
    if outparams.fitActiveFlag && ~outparams.fitTogetherFlag
        % Prepare for active-parameters-fitting
        fprintf('Performing active-parameters-fitting for cell %s ... \n', cellName);
        tStartActiveFit = tic();                % tracks time for active-parameters-fitting

        % Prepare outparams for active fit
        [outparams, neuronParamsUseOrig] = prepare_outparams_active(outparams);

        % Optimize active parameters by fitting to IPSC response
        initialConditionsAll = cell(1, nInitConds);                % stores the initial parameters for each simplex run
        simplexOutAll = cell(1, nInitConds);        % stores the simplex outputs for each simplex run
        exitFlagAll = zeros(1, nInitConds);         % stores the exitflags for each simplex run
%        parfor iInitCond = 1:nInitConds
        for iInitCond = 1:nInitConds
            % Initialize outparams0 for parfor
            outparams0 = outparams;

            % Store the old simplex count
            oldSimplexCt = outparams.simplexNum;

            % Prepare outparams0 for simplex
            [outparams0] = ...
                prepare_outparams_simplex(outparams0, oldOutFolderName, oldSimplexCt, iInitCond);

            % Generate initial parameters for this run
            [initialConditionsAll{iInitCond}, outparams0] = ...
                generate_IC(outparams0, iInitCond);

            % Run fminsearch for GABAB IPSC response
            %%%%%%
            %%%%%%%%%%%%%
            [simplexOutAll{iInitCond}, exitFlagAll(iInitCond)] = ...
                m3ha_fminsearch3(realDataIpscr, outparams0);
            %%%%%%%%%%%%%
            %%%%%%
        end

        % Restore outparams after active fit
        outparams = restore_outparams_active(outparams, neuronParamsUseOrig);

        % Update # of simplex runs and # of simplex steps
        outparams.simplexNum = outparams.simplexNum + nInitConds;
        outparams.simplexIterCount = ...
            outparams.simplexIterCount + ...
            sum(cellfun(@(x) x.ctIterations, simplexOutAll));

        % Find best of the optimized parameters
        compareparamsfile = fullfile(oldOutFolderName, ...
                                    [outparams.prefix, '_compare_params_IC_', ...
                                    num2str(outparams.simplexNum - nInitConds + 1), ...
                                    'to', num2str(outparams.simplexNum), '.csv']);
        [simplexOutBest, err] = ...
            find_best_params(simplexOutAll, initialConditionsAll, compareparamsfile);

        % Update outparams.neuronParamsTable to the best of the optimized parameters
        outparams.neuronParamsTable = simplexOutBest.neuronParamsTable;

        % Save outputs in matfile
        if outparams.saveMatFileFlag
            save(fullfile(oldOutFolderName, ...
                            [outparams.prefix, '_runauto_results_IC_', ...
                            num2str(outparams.simplexNum - nInitConds + 1), ...
                            'to', num2str(outparams.simplexNum), '.mat']), ...
                'initialConditionsAll', 'simplexOutAll', 'exitFlagAll', ...
                'simplexOutBest', 'err', ...
                'outparams', '-v7.3');
        end

        % Print time elapsed
        timeTakenActiveFit = toc(tStartActiveFit);
        fprintf('It took %g seconds to run the simplex method on %d initial conditions!!\n', ...
                    timeTakenActiveFit, nInitConds);
        fprintf('\n');
    else
        if outparams.runnumTotal > 1    % if this is not the first run
            % Use err from last run
            err = outparams.err{outparams.runnumTotal - 1};
        else                            % if this is the first run
            % Use empty structure
            err = struct;
        end
    end
end

% Store error structure in outparams
outparams.errCpr{outparams.runnumTotal} = errCpr;
outparams.err{outparams.runnumTotal} = err;        

% Restore flags for stats and plots and parameter usage
[outparams] = ...
    restore_fields(outparams, ...
        'ltsBurstStatsFlag', 'saveLtsInfoFlag', 'saveLtsStatsFlag', ...
        'saveSimCmdsFlag', 'saveStdOutFlag', 'saveSimOutFlag', ...
        'plotIndividualFlag', 'plotOverlappedFlag', 'plotResidualsFlag', ...
        'plotConductanceFlag', 'plotCurrentFlag', ...
        'plotIpeakFlag', 'plotLtsFlag', 'plotStatisticsFlag', ...
        'plotSwpWeightsFlag');

% Update error history plot
update_errorhistoryplot(hfig, outparams);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errCpr, err, outparams, hfig] = runauto_w_jitter (outparams, hfig)
%% AUTO WITH JITTER mode: TODO

% Update run counts
outparams.runnumAutoWithJitter = outparams.runnumAutoWithJitter + 1;
outparams.runnumTotal = outparams.runnumTotal + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outparams, prefixOrig, neuronParamsUseOrig] = prepare_outparams_passive (outparams)
%% Prepare outparams for passive fit

% Turn flag for passive-parameters-fitting on
outparams.simMode = 'passive';          

% Change prefix for passive-parameters-fitting
prefixOrig = outparams.prefix;
outparams.prefix = [outparams.prefix, '_cpr'];

% Save original parameter usage
neuronParamsUseOrig = outparams.neuronParamsTable{:, 'InUse'};

% Look for all active parameters
indParamsIsActive = find(~outparams.neuronParamsTable.IsPassive);

% Turn off active parameters
outparams.neuronParamsTable{indParamsIsActive, 'InUse'} = ...
    zeros(length(indParamsIsActive), 1);

% Set simplexParams to the passive ones
outparams.simplexParams = outparams.simplexParamsPassive;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outparams] = restore_outparams_passive (outparams, prefixOrig, neuronParamsUseOrig)
%% Restore outparams after passive fit

% Turn flag for passive-parameters-fitting off
outparams.simMode = 'active';

% Restore prefix
outparams.prefix = prefixOrig;

% Restore original parameter usage
outparams.neuronParamsTable{:, 'InUse'} = neuronParamsUseOrig;

% Reset outparams.simplexParams for safety
outparams.simplexParams = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outparams, neuronParamsUseOrig] = prepare_outparams_active (outparams)
%% Prepare outparams for active fit

% Save original parameter usage
neuronParamsUseOrig = outparams.neuronParamsTable.InUse;

% Look for all passive parameters
indParamsIsPassive = find(outparams.neuronParamsTable.IsPassive);

% Turn off passive parameters
outparams.neuronParamsTable{indParamsIsPassive, 'InUse'} = ...
    zeros(length(indParamsIsPassive), 1);

%Set simplexParams to the active ones
outparams.simplexParams = outparams.simplexParamsActive;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outparams = restore_outparams_active (outparams, neuronParamsUseOrig)
%% Restore outparams after active fit

% Restore original parameter usage
outparams.neuronParamsTable{:, 'InUse'} = neuronParamsUseOrig;

% Reset outparams.simplexParams for safety
outparams.simplexParams = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outparams = prepare_outparams_simplex (outparams, oldOutFolderName, ...
                                                oldSimplexCt, initCondNum)
% Prepare outparams for simplex

% Update simplex number
outparams.simplexNum = oldSimplexCt + initCondNum;

% Create a subfolder for simplex outputs and update outFolderName
outparams.outFolder = ...
    fullfile(oldOutFolderName, ['simplex_', num2str(outparams.simplexNum)]);

% Make sure the directory exists
check_dir(outparams.outFolder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [initCond, outparams] = generate_IC (outparams, initCondNum)
%% Generate randomized initial conditions, except for the first run

% Retrieve from outparams
prevNeuronParamsTable = outparams.neuronParamsTable;
simplexNum = outparams.simplexNum;
% RinEstimate = outparams.RinEstimate;

% Set up initial parameters for this initial condition
if initCondNum == 1
    % Just use the previous set of parameters in the first run
    newNeuronParamsTable = prevNeuronParamsTable;
else
    % Create a unique seed for this set of simplex runs
    rng(simplexNum, 'twister');

    % Randomize parameters for the rest of runs
    newNeuronParamsTable = ...
        m3ha_neuron_create_new_initial_params(prevNeuronParamsTable);
end

%% Generate new parameter names for the initCond structure
% TODO: Might not be necessary
% Determine whether each parameter need to be changed
isInUse = newNeuronParamsTable{:, 'InUse'};

% Get the names of all parameters in use
paramsInUseNames = newNeuronParamsTable.Properties.RowNames(isInUse);

% Modify parameter names for storage
paramsInUseNames = cellfun(@(x) strcat(x, '_0'), paramsInUseNames, ...
                            'UniformOutput', false);

% Get the values of all parameters in use
paramsInUseNewValue = newNeuronParamsTable{isInUse, 'InitValue'};

%% Update outparams structure
% Update outparams to these initial parameters
outparams.neuronParamsTable = newNeuronParamsTable;

%% Store in initCond structure
% Store the seed of the random number generator
initCond.randomSeed = rng;

% Store everything else
initCond.initCondNum = initCondNum;
initCond.simplexNum = simplexNum;
initCond.paramsInUseNames = paramsInUseNames;
initCond.paramsInUseValues = paramsInUseNewValue;
initCond.neuronParamsTable = newNeuronParamsTable;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [simplexOutBest, errBest] = ...
        find_best_params(simplexOutAll, initialConditionsAll, compareparamsfile)
%% Find best parameters

% Count the number of initial conditions
nInitConds = numel(simplexOutAll);                    % # of simplex runs

% Count the number of parameters in use
nInUse = length(simplexOutAll{1}.paramsInUseValues);   % # of parameters in use

% Print parameters for all simplex runs in output csv file
fid = fopen(compareparamsfile, 'w');    % open csv file for writing
fprintf(fid, '%s, %s, %s, ', ...
        'Iteration #', 'Final total error', 'Initial total error');
fprintf(fid, repmat('%s, ', 1, nInUse), simplexOutAll{1}.paramsInUseNames{:});
fprintf(fid, repmat('%s, ', 1, nInUse), initialConditionsAll{1}.paramsInUseNames{:});
fprintf(fid, '%s, %s, %s, %s, %s, %s, %s, ', ...
        'Final error change', 'Error tolerance', ...
        'Final parameter change', 'Parameter change tolerance', ...
        'Total iterations', 'Total function evaluations', ...
        'Last simplex step');
fprintf(fid, '%s, %s, ', ...
        'First error change', 'First parameter change');
fprintf(fid, '%s, %s\n', 'Maximum iterations', 'Maximum function evaluations');
fprintf('\n');
for iInitCond = 1:nInitConds
    fprintf(fid, '%d, %g, %g, ', ...
            iInitCond, simplexOutAll{iInitCond}.totalError, simplexOutAll{iInitCond}.initError);
    fprintf(fid, repmat('%g, ', 1, nInUse), simplexOutAll{iInitCond}.paramsInUseValues);
    fprintf(fid, repmat('%g, ', 1, nInUse), initialConditionsAll{iInitCond}.paramsInUseValues);
    fprintf(fid, '%g, %g, %g, %g, %d, %d, %s, ', ...
            simplexOutAll{iInitCond}.maxErrorChange, simplexOutAll{iInitCond}.relativeErrorTolerance, ...
            simplexOutAll{iInitCond}.maxParamChange, simplexOutAll{iInitCond}.relativeParamTolerance, ...
            simplexOutAll{iInitCond}.ctIterations, simplexOutAll{iInitCond}.ctEvals, ...
            simplexOutAll{iInitCond}.lastHow);
    fprintf(fid, '%g, %g, ', ...
            simplexOutAll{iInitCond}.firstMaxErrorChange, ...
            simplexOutAll{iInitCond}.firstMaxParamChange);
    fprintf(fid, '%d, %d\n', ...
            simplexOutAll{iInitCond}.maxIterations, simplexOutAll{iInitCond}.maxFunctionEvaluations);
end
fclose(fid);                            % close csv file

% Convert simplexOutAll to a cell array of row vectors
[simplexOutVecs, vecNames] = structs2vecs(simplexOutAll); 

% Find the index of the total error vector
idxTotalError = find_in_strings('totalError', vecNames, ...
                                    'Searchmode', 'exact');

% Extract the total error vector
simplexOutTotalErrors = simplexOutVecs{idxTotalError};

% Find index of run with smallest final total error
[~, best] = min(simplexOutTotalErrors);    

% Output the best parameters and the associated simplex outputs
simplexOutBest = simplexOutAll{best};
errBest = simplexOutBest.err;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_errorhistoryplot(hfig, outparams)
%% Shows and updates Error History figure

rn = outparams.runnumTotal;     % current run number
if outparams.fitActiveFlag
    err = outparams.err{rn};        % current error structure
elseif outparams.fitPassiveFlag
    err = outparams.errCpr{rn};
end

% Make the plot the current figure
if outparams.fitActiveFlag
    set(0, 'CurrentFigure', hfig.errorhistory);
elseif outparams.fitPassiveFlag
    set(0, 'CurrentFigure', hfig.cprerrorhistory);
end

% Plot the error
if ~isempty(err)
    if outparams.fitActiveFlag
        % Plot the total error
        subplot(3, 2, 1);
        update_subplot(rn, err.totalError, [], 'total error', 'o', 'b');

        % Plot the average sweep error
        subplot(3, 2, 2);
        update_subplot(rn, err.avgSwpError, [], 'sweep error', 'o', 'b');

        % Plot the average LTS error
        subplot(3, 2, 3);
        update_subplot(rn, err.avgLtsError, [], 'LTS error', 'o', 'b');

        % Plot the average LTS amp error
        subplot(3, 2, 4);
        update_subplot(rn, err.avgLtsAmpError, [], 'amp error', 'o', 'b');

        % Plot the average LTS time error
        subplot(3, 2, 5);
        update_subplot(rn, err.avgLtsDelayError, [], 'time error', 'o', 'b');

        % Plot the average LTS slope error
        subplot(3, 2, 6);
        update_subplot(rn, err.avgLtsSlopeError, [], 'slope error', 'o', 'b');
    else                % if no lts error
        % Just plot the total error
        update_subplot(rn, err.totalError, 'run number', 'total error', 'o', 'b');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_subplot(x, y, x_label, y_label, markerstyle, color)
% Update subplot

% Plot new marker
plot(x, y, 'Marker', markerstyle, 'Color', color, 'MarkerFaceColor', 'auto');

% Adjust y limits
if x == 1
    hold on;
    initymax = y * 1.1;             % initial ymax
    ylim([0 initymax]);
    if ~isempty(x_label)
        xlabel(x_label); 
    end
    if ~isempty(y_label)
        ylabel(y_label);
    end
else
    ylimits = get(gca, 'YLim');
    if ylimits(2) < y               % rescale axes if error is greater than axis limit
        ylim([0, y * 1.1]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%