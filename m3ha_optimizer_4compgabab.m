function [done, outParams, hfig] = m3ha_optimizer_4compgabab (outParams, hfig)
% OPTIMIZER  Passes parameters to NEURON (which runs simulations and saves
% data as .dat files), loads .dat file data, plots data and compares with
% experimental data. Save all this as well (timestamp for now)
%
% realData  Experimental data from OPTIMIZERGUI.m
% 
% outParams   GUI parameters from OPTIMIZERGUI.m
%   outParams.swpWeights : weights of sweeps
%   outParams.swpedges : left and right edges of sweep region to analyze
%   outParams.ltsWeights : weights of LTS parameters
%   outParams.modeselected : mode to run OPTIMIZER.m in (either
%       'modebutton_auto' or 'modebutton_manual')
%   outParams.neuronparams : values of parameters for NEURON 
%   outParams.neuronparamnames : names of parameters for NEURON (e.g. pcabar1,
%       actshift_itGHK, etc)
%   outParams.simplexParams : values of simplex parameters for
%   FMINSEARCH2.m
%   outParams.simplexParamNames : names of simplex parameters for
%   FMINSEARCH2.m
%   outParams.sortedswpnum : sweep numbers of sweeps called into
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
%
% Used by:
%       cd/singleneuronfitting42.m and later versions
%       cd/m3ha_optimizergui_4compgabab.m

% File History:
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
% 2016-10-06 - xlimits now uses outParams.fitreg instead of outParams.swpedges
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
% 2017-05-02 - Changed simplex count to outParams.simplexNum and 
%               fixed the fact that it wasn't updated during active fitting
% 2017-05-12 - Now saves runauto_results in outfolder
% 2017-05-13 - Now updates outParams.runnumTotal here
% 2017-05-13 - Now returns errCpr & err in outParams
% 2017-05-13 - Now gets outFolderName from outParams
% 2017-05-13 - Now updates outParams.neuronparams to 
%               initCond.neuronparams in generate_IC()
% 2017-05-13 - Now initializes outparams0 within parfor loop
% 2017-05-15 - Now changes outParams.prefix so that '_cpr' is already 
%               incorporated before passive fitting
% 2017-05-15 - Now changes neuronparams_use before calling fminsearch3_4compgabab 
%               to reflect passive or active
% 2017-05-16 - Now runs NEURON with full plots both before and after fitting
% 2017-05-16 - parfor is now conditional on outParams.MaxNumWorkersIC
% 2017-05-17 - update_errorhistoryplot() now plots lts errors if computed
% 2017-05-17 - Added outParams.fitCprFlag & outParams.fitIpscrFlag
% 2017-05-17 - Moved update_sweeps_figures() to m3ha_neuron_run_and_analyze.m
% 2017-05-19 - Fixed the fact that simplex.initError was used to select  
%               best simplex; simplex.error is now simplex.totalError
% 2017-05-22 - Changed line width and indentation
% 2017-05-22 - Added outParams.oldFitTogetherFlag
% 2017-05-23 - Removed modeselected from outParams and 
%               replaced with updated outParams.runMode
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
% 2018-08-10 - Changed outParams.fitregCpr to use outParams.fitwinCpr
% 2018-08-16 - Now does not use parfor for active fitting
% 2018-12-10 - Updated placement of jitterFlag
% 2018-12-11 - Moved code to m3ha_neuron_create_new_initial_params.m
% 2019-11-13 - Now copies final parameters without renaming
% 2019-11-15 - Now passes in outParams for varargin 
%                   for m3ha_neuron_run_and_analyze.m
% 2019-11-18 - Fixe fields that are set to zero during fitting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Change fitting flags if necessary
% TODO: Fix
%{
if outParams.fitCprFlag && ...
    ~any(neuronparamispas .* outParams.neuronparams_use)
    outParams.fitCprFlag = false;
    fprintf(['No passive parameter is fitted, ', ...
            'so passive fitting is turned OFF!\n\n']);
end
if outParams.fitIpscrFlag && ...
    ~any(~neuronparamispas .* outParams.neuronparams_use)
    outParams.fitIpscrFlag = false;
    fprintf(['No active parameter is fitted, ', ...
            'so active fitting is turned OFF!\n\n']);
end
%}

% Initialize the prefix as the date-time-cell stamp
if ~isfield(outParams, 'prefix')
    outParams.prefix = outParams.dateTimeCellStamp;
end

% If in jitter mode and if parameter is checked, add jitter
if outParams.runMode == 3
    outParams.jitterFlag = true;
else
    outParams.jitterFlag = false;
end

%% Run based on manual or auto mode
switch outParams.runMode
case {1, 3}
    [~, ~, outParams, hfig] = runmanual(outParams, hfig);
case 2
    %% Save old NEURON parameters
    oldNeuronParams = outParams.neuronParamsTable.Value;

    %% Run NEURON with baseline parameters with full plots
    prefixOrig = outParams.prefix;                         % save original prefix
    outParams.prefix = [outParams.prefix, '_bef'];          % change prefix for "before fitting"
    [~, ~, outParams, hfig] = runmanual(outParams, hfig);
    outParams.prefix = prefixOrig;                         % restore original prefix
    drawnow

    %% Optimize parameters
    [~, ~, outParams, hfig] = runauto(outParams, hfig);
    drawnow

    %% Run NEURON with best parameters again with full plots
    prefixOrig = outParams.prefix;                         % save original prefix
    outParams.prefix = [outParams.prefix, '_aft'];          % change prefix for "after fitting"
    [~, ~, outParams, hfig] = runmanual(outParams, hfig);
    outParams.prefix = prefixOrig;                         % restore original prefix
    drawnow

    %% If error is not worse, copy final parameters to bestParamsDirectory
    if outParams.err{outParams.runnumTotal}.totalError <= ...
                outParams.err{1}.totalError       
        % Get the final parameters file name for this cell
        finalParamsFile = fullfile(outParams.outFolder, ...
                            [outParams.prefix, '_aft_params.csv']);

        % Copy to bestParamsDirectory
        copyfile(finalParamsFile, outParams.bestParamsDirectory);
        fprintf('Best parameters copied for %s!!\n\n', outParams.prefix);
    end

    %% Compare NEURON parameters before and after 
    %   by plotting activation/inactivation & I-V curves
    newNeuronParams = outParams.neuronParamsTable.Value;
    neuronParamNames = outParams.neuronParamsTable.Properties.RowNames;
    bothParamValues = {oldNeuronParams, newNeuronParams};
    bothParamNames = {neuronParamNames, neuronParamNames};
    suffixes = {['_', outParams.prefix, '_bef'], ...
                ['_', outParams.prefix, '_aft']};
% TODO: Fix
%     m3ha_compare_neuronparams2(bothParamValues, bothParamNames, suffixes, ...
%                         'OutFolder', outParams.outFolder);

case 4                                       %%% TODO: Unfinished
    [~, ~, outParams, hfig] = ...
        runauto_w_jitter(outParams, hfig);
otherwise
    outParams.runMode = 1;
    fprintf('Warning: run mode out of range, changed to 1!\n\n');
end

%% Used for runbutton_toggle in m3ha_optimizergui_4compgabab.m
done = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errCpr, err, outParams, hfig] = runmanual(outParams, hfig)
%% MANUAL or JITTER modes: Simulate each sweep once and compare with real data

% Update run counts
switch outParams.runMode
case {1, 2}         % manual mode is also run before and after auto mode
    outParams.runnumManual = outParams.runnumManual + 1;
case 3
    outParams.runnumJitter = outParams.runnumJitter + 1;
otherwise
    error('runmode should be 1 or 3!\n');
end
outParams.runnumTotal = outParams.runnumTotal + 1;

% Simulate current pulse response only if outParams.fitCprFlag is true
if outParams.fitCprFlag         % if fitting passive parameters
    %%%%%% 
    %%%%%%%%%%%%%
    % Set sim mode
    outParams.simMode = 'passive';

    % Save original prefix
    prefixOrig = outParams.prefix;

    % Change prefix for passive-parameters-fitting
    outParams.prefix = [outParams.prefix, '_cpr'];

    % Set grouping vector to be vHold
    outParams.grouping = outParams.vHold;

    % Run all simulations once
    [errCpr, outParams.hfig] = ...
        m3ha_neuron_run_and_analyze(outParams.neuronParamsTable, outParams);

    % Change prefix back
    outParams.prefix = prefixOrig;
    %%%%%%%%%%%%%
    %%%%%% 

    % Log errors and parameters and save parameters as .p file
    if ~isemptystruct(errCpr)
        m3ha_log_errors_params(outParams.logFileName, outParams, errCpr);
    end
else
    if outParams.runnumTotal > 1    % if this is not the first run
        % Use errCpr from last run
        errCpr = outParams.errCpr{outParams.runnumTotal-1};   
    else                            % if this is the first run
        % Use empty structure
        errCpr = struct;
    end
end

% Store error structure in outParams
outParams.errCpr{outParams.runnumTotal} = errCpr;

% Simulate GABAB IPSC response only if outParams.fitIpscrFlag is true
if outParams.fitIpscrFlag
    %%%%%% 
    %%%%%%%%%%%%%
    % Set sim mode
    outParams.simMode = 'active';

    % Run all simulations once
    [err, outParams.hfig] = ...
        m3ha_neuron_run_and_analyze(outParams.neuronParamsTable, outParams);              
    %%%%%%%%%%%%%
    %%%%%% 

    % Log errors and parameters and save parameters
    if ~isemptystruct(err)
        m3ha_log_errors_params(outParams.logFileName, outParams, err);
    end
else
    if outParams.runnumTotal > 1    % if this is not the first run
        % Use errCpr from last run
        err = outParams.err{outParams.runnumTotal-1};   
    else                            % if this is the first run
        % Use empty structure
        err = struct;
    end
end

% Store error structure in outParams
outParams.err{outParams.runnumTotal} = err;

% Update error history plot
% TODO: Fix
% update_errorhistoryplot(hfig, outParams);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errCpr, err, outParams, hfig] = runauto(outParams, hfig)
%% AUTO mode: Find best parameters with m3ha_fminsearch3.m

% Update run counts
outParams.runnumAuto = outParams.runnumAuto + 1;
outParams.runnumTotal = outParams.runnumTotal + 1;

% Extract constants from outParams
cellName = outParams.cellName;
oldOutFolderName = outParams.outFolder;

% Extract number of initial conditions
idxNInitConds = find_in_strings('nInitConds', outParams.autoParamNames);
nInitConds = outParams.autoParams(idxNInitConds);

% Turn off all flags for stats and plots for fminsearch
[outParams] = ...
    set_fields_zero(outParams, ...
        'saveLtsInfoFlag', 'saveLtsStatsFlag', ...
        'saveSimCmdsFlag', 'saveStdOutFlag', 'saveSimOutFlag', ...
        'plotConductanceFlag', 'plotCurrentFlag', ...
        'plotIndividualFlag', 'plotResidualsFlag', 'plotOverlappedFlag', ...
        'plotIpeakFlag', 'plotLtsFlag', 'plotStatisticsFlag', ...
        'plotSwpWeightsFlag');

% Fit parameters
if outParams.fitMode == 5
    % Don't run current pulse response
    if outParams.runnumTotal > 1
        % If this is not the first run, use errCpr from last run
        errCpr = outParams.errCpr{outParams.runnumTotal - 1};
    else
        % If this is the first run, use empty structure
        errCpr = struct;
    end

    % Prepare for active-parameters-fitting
    fprintf('Fitting to IPSC responses for cell %s ... \n', cellName);
    tStartActiveFit = tic();

    % Set simplexParams to the active ones
    outParams.simplexParams = outParams.simplexParamsActive;

    % Optimize active parameters by fitting to IPSC response
    initialConditionsAll = cell(1, nInitConds);
    simplexOutAll = cell(1, nInitConds);
    exitFlagAll = zeros(1, nInitConds);
%    parfor iInitCond = 1:nInitConds
    for iInitCond = 1:nInitConds
        % Initialize outparams0 for parfor
        outparams0 = outParams;

        % Store the old simplex count
        oldSimplexCt = outParams.simplexNum;

        % Prepare outparams0 for simplex
        outparams0 = prepare_outparams_simplex(outparams0, oldOutFolderName, ...
                                                oldSimplexCt, iInitCond);

        % Generate initial parameters for this run
        [initialConditionsAll{iInitCond}, outparams0] = ...
            generate_IC(outparams0, iInitCond);

        % Run fminsearch for GABAB IPSC response
        %%%%%%
        %%%%%%%%%%%%%
        [simplexOutAll{iInitCond}, exitFlagAll(iInitCond)] = ...
            m3ha_fminsearch3(outparams0);
        %%%%%%%%%%%%%
        %%%%%%
    end

    % Reset outParams.simplexParams for safety
    outParams.simplexParams = [];

    % Update # of simplex runs and # of simplex steps
    outParams.simplexNum = outParams.simplexNum + nInitConds;
    outParams.simplexIterCount = ...
        outParams.simplexIterCount + ...
        sum(cellfun(@(x) x.ctIterations, simplexOutAll));

    % Find best of the optimized parameters
    compareparamsfile = fullfile(oldOutFolderName, ...
                            [outParams.prefix, '_compare_params_IC_', ...
                            num2str(outParams.simplexNum - nInitConds + 1), ...
                            'to', num2str(outParams.simplexNum), '.csv']);
    [simplexOutBest, err] = ...
        find_best_params(simplexOutAll, initialConditionsAll, ...
                        compareparamsfile);

    % Update outParams.neuronParamsTable to the best of the optimized parameters
    outParams.neuronParamsTable = simplexOutBest.neuronParamsTable;

    % Save outputs in matfile
    if outParams.saveMatFileFlag
        save(fullfile(oldOutFolderName, ...
                        [outParams.prefix, '_runauto_results_IC_', ...
                        num2str(outParams.simplexNum - nInitConds + 1), ...
                        'to', num2str(outParams.simplexNum), '.mat']), ...
            'initialConditionsAll', 'simplexOutAll', 'exitFlagAll', ...
            'simplexOutBest', 'err', ...
            'outParams', '-v7.3');
    end

    % Print time elapsed
    timeTakenActiveFit = toc(tStartActiveFit);
    fprintf('It took %g seconds to run the simplex method on %d initial conditions!!\n', ...
                timeTakenActiveFit, nInitConds);
    fprintf('\n');
elseif outParams.fitMode == 4
    % Passive and active fitting done together from single initial conditions
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
        outparams0 = outParams;

        % Store the old simplex count
        oldSimplexCt = outParams.simplexNum;

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
            m3ha_fminsearch3(outparams0);
        %%%%%%%%%%%%%
        %%%%%%

        % Restore outparams0 after active fit
        outparams0 = restore_outparams_active(outparams0, neuronParamsUseOrig);
    end

    % Update # of simplex runs and # of simplex steps
    outParams.simplexNum = outParams.simplexNum + nInitConds * 2;
    outParams.simplexIterCount = ...
        outParams.simplexIterCount + ...
        sum(cellfun(@(x) x.ctIterations, cprSimplexOutAll)) + ...
        sum(cellfun(@(x) x.ctIterations, simplexOutAll));

    % Find best of the optimized parameters
    cprCompareparamsfile = fullfile(oldOutFolderName, ...
                                    [outParams.prefix, '_cpr_compare_params_IC_', ...
                                    num2str(outParams.simplexNum - nInitConds * 2 + 1), ...
                                    'to', num2str(outParams.simplexNum - nInitConds), '.csv']);
    [cprSimplexOutBest, errCpr] = ...
        find_best_params(cprSimplexOutAll, cprInitialConditionsAll, cprCompareparamsfile);
    compareparamsfile = fullfile(oldOutFolderName, ...
                                [outParams.prefix, '_compare_params_IC_', ...
                                num2str(outParams.simplexNum - nInitConds + 1), ...
                                'to', num2str(outParams.simplexNum), '.csv']);
    [simplexOutBest, err] = ...
        find_best_params(simplexOutAll, initialConditionsAll, compareparamsfile);

    % Update outParams.neuronParamsTable to the best of the optimized parameters
    outParams.neuronParamsTable = simplexOutBest.neuronParamsTable;

    % Save outputs in matfile
    if outParams.saveMatFileFlag
        save(fullfile(oldOutFolderName, ...
                        [outParams.prefix, '_runauto_results_IC_', ...
                        num2str(outParams.simplexNum - nInitConds + 1), ...
                        'to', num2str(outParams.simplexNum), '.mat']), ...
            'cprInitialConditionsAll', 'cprSimplexOutAll', 'cprExitFlagAll', ...
            'initialConditionsAll', 'simplexOutAll', 'exitFlagAll', ...
            'outParams', '-v7.3');
    end


    % Print time elapsed
    timeTakenAllFit = toc(tstartAllFit);
    fprintf('It took %g seconds to run the simplex method on %d initial conditions!!\n', ...
            timeTakenAllFit, nInitConds);
    fprintf('\n');
else
    % Do passive fitting, find best params, then do active fitting
    % Fit passive parameters
    if outParams.fitMode == 1 || outParams.fitMode == 3
        % Prepare for passive-parameters-fitting
        fprintf('Fitting passive parameters for cell %s ... \n', cellName);
        tStartPassiveFit = tic();      % tracks time for passive-parameters-fitting
        fprintf('\n');

        % Prepare outParams for passive fit
        [outParams, prefixOrig, neuronParamsUseOrig] = ...
            prepare_outparams_passive(outParams);
        
        % Optimize passive parameters by fitting to current pulse response
        cprInitialConditionsAll = cell(1, nInitConds);            % stores the initial parameters for each simplex run
        cprSimplexOutAll = cell(1, nInitConds);    % stores the simplex outputs for each simplex run
        cprExitFlagAll = zeros(1, nInitConds);     % stores the exitflags for each simplex run
%        parfor iInitCond = 1:nInitConds
        for iInitCond = 1:nInitConds
            % Initialize outparams0 for parfor
            outparams0 = outParams;

            % Store the old simplex count
            oldSimplexCt = outParams.simplexNum;

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

        % Restore outParams after passive fit
        [outParams] = ...
            restore_outparams_passive(outParams, prefixOrig, neuronParamsUseOrig);

        % Update # of simplex runs and # of simplex steps
        outParams.simplexNum = outParams.simplexNum + nInitConds;    
        outParams.simplexIterCount = outParams.simplexIterCount + ...
                                        sum(cellfun(@(x) x.ctIterations, cprSimplexOutAll));

        % Find best of the optimized parameters
        cprCompareparamsfile = fullfile(oldOutFolderName, ...
                                        [outParams.prefix, '_compare_params_IC_', ...
                                        num2str(outParams.simplexNum - nInitConds + 1), ...
                                        'to', num2str(outParams.simplexNum), '.csv']);
        [cprSimplexOutBest, errCpr] = ...
            find_best_params(cprSimplexOutAll, cprInitialConditionsAll, cprCompareparamsfile);

        % Update outParams.neuronParamsTable to the best of the optimized parameters
        outParams.neuronParamsTable = cprSimplexOutBest.neuronParamsTable;

        % Save outputs in matfile
        if outParams.saveMatFileFlag
            save(fullfile(oldOutFolderName, ...
                            [outParams.prefix, '_runauto_results_IC_', ...
                            num2str(outParams.simplexNum - nInitConds + 1), ...
                            'to', num2str(outParams.simplexNum), '.mat']), ...
                'cprInitialConditionsAll', 'cprSimplexOutAll', ...
                'cprExitFlagAll', 'cprSimplexOutBest', 'errCpr', ...
                'outParams', '-v7.3');
        end

        % Print time elapsed
        timeTakenPassiveFit = toc(tStartPassiveFit);
        fprintf('It took %g seconds to run the simplex method on %d initial conditions!!\n', ...
                timeTakenPassiveFit, nInitConds);
        fprintf('\n');
    else
        if outParams.runnumTotal > 1    % if this is not the first run
            errCpr = outParams.errCpr{outParams.runnumTotal-1};       % use errCpr from last run
        else                            % if this is the first run
            errCpr = struct;                                          % empty structure
        end
    end

    % Fit active parameters
    if outParams.fitMode == 2 || outParams.fitMode == 3
        % Prepare for active-parameters-fitting
        fprintf('Performing active-parameters-fitting for cell %s ... \n', cellName);
        tStartActiveFit = tic();                % tracks time for active-parameters-fitting

        % Prepare outParams for active fit
        [outParams, neuronParamsUseOrig] = prepare_outparams_active(outParams);

        % Optimize active parameters by fitting to IPSC response
        initialConditionsAll = cell(1, nInitConds);                % stores the initial parameters for each simplex run
        simplexOutAll = cell(1, nInitConds);        % stores the simplex outputs for each simplex run
        exitFlagAll = zeros(1, nInitConds);         % stores the exitflags for each simplex run
%        parfor iInitCond = 1:nInitConds
        for iInitCond = 1:nInitConds
            % Initialize outparams0 for parfor
            outparams0 = outParams;

            % Store the old simplex count
            oldSimplexCt = outParams.simplexNum;

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
                m3ha_fminsearch3(outparams0);
            %%%%%%%%%%%%%
            %%%%%%
        end

        % Restore outParams after active fit
        outParams = restore_outparams_active(outParams, neuronParamsUseOrig);

        % Update # of simplex runs and # of simplex steps
        outParams.simplexNum = outParams.simplexNum + nInitConds;
        outParams.simplexIterCount = ...
            outParams.simplexIterCount + ...
            sum(cellfun(@(x) x.ctIterations, simplexOutAll));

        % Find best of the optimized parameters
        compareparamsfile = fullfile(oldOutFolderName, ...
                                    [outParams.prefix, '_compare_params_IC_', ...
                                    num2str(outParams.simplexNum - nInitConds + 1), ...
                                    'to', num2str(outParams.simplexNum), '.csv']);
        [simplexOutBest, err] = ...
            find_best_params(simplexOutAll, initialConditionsAll, compareparamsfile);

        % Update outParams.neuronParamsTable to the best of the optimized parameters
        outParams.neuronParamsTable = simplexOutBest.neuronParamsTable;

        % Save outputs in matfile
        if outParams.saveMatFileFlag
            save(fullfile(oldOutFolderName, ...
                            [outParams.prefix, '_runauto_results_IC_', ...
                            num2str(outParams.simplexNum - nInitConds + 1), ...
                            'to', num2str(outParams.simplexNum), '.mat']), ...
                'initialConditionsAll', 'simplexOutAll', 'exitFlagAll', ...
                'simplexOutBest', 'err', ...
                'outParams', '-v7.3');
        end

        % Print time elapsed
        timeTakenActiveFit = toc(tStartActiveFit);
        fprintf('It took %g seconds to run the simplex method on %d initial conditions!!\n', ...
                    timeTakenActiveFit, nInitConds);
        fprintf('\n');
    else
        if outParams.runnumTotal > 1    % if this is not the first run
            % Use err from last run
            err = outParams.err{outParams.runnumTotal - 1};
        else                            % if this is the first run
            % Use empty structure
            err = struct;
        end
    end
else
    % Use empty structures
    errCpr = struct;
    err = struct;
end

% Store error structure in outParams
outParams.errCpr{outParams.runnumTotal} = errCpr;
outParams.err{outParams.runnumTotal} = err;        

% Restore flags for stats and plots and parameter usage
[outParams] = ...
    restore_fields(outParams, ...
        'saveLtsInfoFlag', 'saveLtsStatsFlag', ...
        'saveSimCmdsFlag', 'saveStdOutFlag', 'saveSimOutFlag', ...
        'plotConductanceFlag', 'plotCurrentFlag', ...
        'plotIndividualFlag', 'plotResidualsFlag', 'plotOverlappedFlag', ...
        'plotIpeakFlag', 'plotLtsFlag', 'plotStatisticsFlag', ...
        'plotSwpWeightsFlag');

% Update error history plot
update_errorhistoryplot(hfig, outParams);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errCpr, err, outParams, hfig] = runauto_w_jitter (outParams, hfig)
%% AUTO WITH JITTER mode: TODO

% Update run counts
outParams.runnumAutoWithJitter = outParams.runnumAutoWithJitter + 1;
outParams.runnumTotal = outParams.runnumTotal + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outParams, prefixOrig, neuronParamsUseOrig] = ...
                prepare_outparams_passive (outParams)
%% Prepare outParams for passive fit

% Turn flag for passive-parameters-fitting on
outParams.simMode = 'passive';

% Change prefix for passive-parameters-fitting
prefixOrig = outParams.prefix;
outParams.prefix = [outParams.prefix, '_cpr'];

% Save original parameter usage
neuronParamsUseOrig = outParams.neuronParamsTable{:, 'InUse'};

% Look for all active parameters
indParamsIsActive = find(~outParams.neuronParamsTable.IsPassive);

% Turn off active parameters
outParams.neuronParamsTable{indParamsIsActive, 'InUse'} = ...
    zeros(length(indParamsIsActive), 1);

% Set simplexParams to the passive ones
outParams.simplexParams = outParams.simplexParamsPassive;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outParams] = restore_outparams_passive (outParams, prefixOrig, ...
                                                neuronParamsUseOrig)
%% Restore outParams after passive fit

% Turn flag for passive-parameters-fitting off
outParams.simMode = 'active';

% Restore prefix
outParams.prefix = prefixOrig;

% Restore original parameter usage
outParams.neuronParamsTable{:, 'InUse'} = neuronParamsUseOrig;

% Reset outParams.simplexParams for safety
outParams.simplexParams = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outParams, neuronParamsUseOrig] = prepare_outparams_active (outParams)
%% Prepare outParams for active fit

% Save original parameter usage
neuronParamsUseOrig = outParams.neuronParamsTable.InUse;

% Look for all passive parameters
indParamsIsPassive = find(outParams.neuronParamsTable.IsPassive);

% Turn off passive parameters
outParams.neuronParamsTable{indParamsIsPassive, 'InUse'} = ...
    zeros(length(indParamsIsPassive), 1);

% Set simplexParams to the active ones
outParams.simplexParams = outParams.simplexParamsActive;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outParams = restore_outparams_active (outParams, neuronParamsUseOrig)
%% Restore outParams after active fit

% Restore original parameter usage
outParams.neuronParamsTable{:, 'InUse'} = neuronParamsUseOrig;

% Reset outParams.simplexParams for safety
outParams.simplexParams = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outParams = prepare_outparams_simplex (outParams, oldOutFolderName, ...
                                                oldSimplexCt, initCondNum)
% Prepare outParams for simplex

% Update simplex number
outParams.simplexNum = oldSimplexCt + initCondNum;

% Create a subfolder for simplex outputs and update outFolderName
outParams.outFolder = ...
    fullfile(oldOutFolderName, ['simplex_', num2str(outParams.simplexNum)]);

% Make sure the directory exists
check_dir(outParams.outFolder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [initCond, outParams] = generate_IC (outParams, initCondNum)
%% Generate randomized initial conditions, except for the first run

% Retrieve from outParams
prevNeuronParamsTable = outParams.neuronParamsTable;
simplexNum = outParams.simplexNum;
% RinEstimate = outParams.RinEstimate;

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

%% Update outParams structure
% Update outParams to these initial parameters
outParams.neuronParamsTable = newNeuronParamsTable;

%% Store in initCond structure
% Store the seed of the random number generator
initCond.randomSeed = rng;

% Store everything else
initCond.initCondNum = initCondNum;
initCond.simplexNum = simplexNum;
initCond.paramsInUseNames = paramsInUseNames;
initCond.paramsInUseValues = paramsInUseNewValue;
initCond.neuronParamsTable = newNeuronParamsTable;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_errorhistoryplot(hfig, outParams)
%% Shows and updates Error History figure

rn = outParams.runnumTotal;     % current run number
if outParams.fitIpscrFlag
    err = outParams.err{rn};        % current error structure
elseif outParams.fitCprFlag
    err = outParams.errCpr{rn};
end

% Make the plot the current figure
if outParams.fitIpscrFlag
    set(0, 'CurrentFigure', hfig.errorhistory);
elseif outParams.fitCprFlag
    set(0, 'CurrentFigure', hfig.cprerrorhistory);
end

% Plot the error
if ~isempty(err)
    if outParams.fitIpscrFlag
        % Extract error ratios from outParams struct
        errorWeights = outParams.errorWeights;
        missedLtsError = outParams.missedLtsError;
        falseLtsError = outParams.falseLtsError;

        % Create y-axis labels
        totalErrorLabel = 'total error';
        avgSwpErrorLabel = sprintf('sweep error (%d)', errorWeights(1));
        ltsMisMatchErrorLabel = sprintf('match error (+%d, %d)', ...
                                missedLtsError, falseLtsError);
        avgLtsAmpErrorLabel = sprintf('amp error (%d)', errorWeights(2));
        avgLtsDelayErrorLabel = sprintf('time error (%d)', errorWeights(3));
        avgLtsSlopeErrorLabel = sprintf('slope error (%d)', errorWeights(4));

        % Plot the total error
        subplot(3, 2, 1);
        update_subplot(rn, err.totalError, [], totalErrorLabel, 'o', 'b');

        % Plot the average sweep error
        subplot(3, 2, 2);
        update_subplot(rn, err.avgSwpError, [], avgSwpErrorLabel, 'o', 'b');

        % Plot the LTS mismatch error
        subplot(3, 2, 3);
        update_subplot(rn, err.ltsMisMatchError, [], ...
                        ltsMisMatchErrorLabel, 'o', 'b');

        % Plot the average LTS amp error
        subplot(3, 2, 4);
        update_subplot(rn, err.avgLtsAmpError, [], ...
                        avgLtsAmpErrorLabel, 'o', 'b');

        % Plot the average LTS time error
        subplot(3, 2, 5);
        update_subplot(rn, err.avgLtsDelayError, [], ...
                        avgLtsDelayErrorLabel, 'o', 'b');

        % Plot the average LTS slope error
        subplot(3, 2, 6);
        update_subplot(rn, err.avgLtsSlopeError, [], ...
                        avgLtsSlopeErrorLabel, 'o', 'b');
    else                % if no lts error
        % Just plot the total error
        update_subplot(rn, err.totalError, 'run number', 'total error', 'o', 'b');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_subplot(x, y, xLabel, yLabel, markerstyle, color)
% Update subplot

% Plot new marker
plot(x, y, 'Marker', markerstyle, 'Color', color, 'MarkerFaceColor', 'auto');

% Adjust y limits
if x == 1
    hold on;
    % Initialize ymax
    initymax = y * 1.1;
    ylim([0, initymax]);
    if ~isempty(xLabel)
        xlabel(xLabel); 
    end
    if ~isempty(yLabel)
        ylabel(yLabel);
    end
else
    ylimits = get(gca, 'YLim');
    % Rescale axes if error is greater than axis limit
    if ylimits(2) < y
        ylim([0, y * 1.1]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%% Make all figures visible and update
if outParams.showSweepsFlag
    figs = fieldnames(hfig);
    for k = 1:numel(figs)
        set(hfig.(figs{k}), 'Visible', 'on');
        drawnow
    end
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
