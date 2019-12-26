function [errorStruct, hFig, simData] = ...
                m3ha_neuron_run_and_analyze (neuronParamsTableOrFile, varargin)
%% Runs and analyzes "one iteration" of NEURON simulations (once for each of the sweeps)
% Usage: [errorStruct, hFig, simData] = ...
%               m3ha_neuron_run_and_analyze (neuronParamsTableOrFile, varargin)
% Explanation:
%       TODO
% 
% Examples:
%       [err, hFig, simData] = m3ha_neuron_run_and_analyze('20191218T1135_D101310_aft_params');
% 
% Outputs: 
%       TODO
%       simData     - simulated data
%                   specified as a numeric array
%                       or a cell array of numeric arrays
% Arguments:
%       neuronParamsTableOrFile
%                   - table(s) of file(s) of single neuron parameters with 
%                       parameter names as 'RowNames' and with variables:
%                       'Value': value of the parameter
%                       'LowerBound': lower bound of the parameter
%                       'UpperBound': upper bound of the parameter
%                       'JitterPercentage': jitter percentage of the parameter
%                       'IsLog': whether the parameter is 
%                                   to be varied on a log scale
%                       Note: Only the 'Value' column is needed
%                               if jitterFlag is false
%                   must be a 2d table or a cell array of 2d tables
%                       or a character vector or a string 
%                       or a cell array of character vectors
%       varargin    - 'Hfig': handles structure for figures
%                   must be a TODO
%                   default == TODO
%                   - 'BuildMode': TC neuron build mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'passive' - insert leak channels only
%                       'active'  - insert both passive and active channels
%                       or a cell array of them TODO
%                   default == 'active'
%                   - 'SimMode': simulation mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'passive' - simulate a current pulse response
%                       'active'  - simulate an IPSC response
%                   default == 'active'
%                   - 'ColumnMode': column mode
%                   must be recognized by m3ha_select_raw_traces.m
%                   default == set in m3ha_select_raw_traces.m
%                   - 'RowMode': row mode
%                   must be recognized by m3ha_select_raw_traces.m
%                   default == set in m3ha_select_raw_traces.m
%                   - 'AttemptNumber': attempt number
%                   must be recognized by m3ha_select_raw_traces.m
%                   default == set in m3ha_select_raw_traces.m
%                   - 'DataMode': data mode
%                   must be recognized by m3ha_select_raw_traces.m
%                   default == set in m3ha_select_raw_traces.m
%                   - 'NSweeps': number of sweeps
%                   must be a positive integer scalar
%                   default == numel(realData) or 1
%                   - 'FileNames': file names provided
%                   must be a character vector, a string vector 
%                       or a cell array of character vectors
%                   default == none provided
%                   - 'FileBases': base of filenames (without extension) 
%                                   corresponding to each vector
%                   must be a character vector, a string vector 
%                       or a cell array of character vectors
%                   default == read from provided data file names
%                               or set in decide_on_filebases.m
%                   - 'OutFolder': the directory where outputs will be placed
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'Prefix': prefix to prepend to file names
%                   must be a character array
%                   default == extract_common_prefix(fileBases)
%                   - 'DebugFlag': whether debugging
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'CustomHoldCurrentFlag': whether to use a custom 
%                                               holding current
%                   must be numeric 1 or 0
%                   default == false
%                   - 'OnHpcFlag': whether on a high performance 
%                                   computing server
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'GenerateDataFlag': whether generating surrogate data
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'AverageCprFlag': whether to average current pulse 
%                                       responses according to vHold
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'BootstrapCprFlag': whether to bootstrap-average current  
%                                       pulse responses according to vHold
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Normalize2InitErrFlag': whether to normalize errors
%                                               to initial errors
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'SaveImportLogFlag': whether to log imported files
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveParamsFlag': whether to save simulation parameters
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveSimCmdsFlag': whether to save simulation commands
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveStdOutFlag': whether to save standard outputs
%                                           when there are no errors
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'SaveSimOutFlag': whether to save simulation outputs
%                                           when there are no errors
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveLtsInfoFlag': whether to save LTS info
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveLtsStatsFlag': whether to save LTS statistics
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotAllFlag': whether to plot all types of plots
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotIndividualFlag': whether to plot individual traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotResidualsFlag': whether to plot residuals
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotOverlappedFlag': whether to plot overlapped traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotConductanceFlag': whether to plot conductance traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotCurrentFlag': whether to plot current traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotIpeakFlag': whether to current peak analyses
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotLtsFlag': whether to plot vtrace/LTS/burst analyses
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotStatisticsFlag': whether to plot LTS statistics
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotSwpWeightsFlag': whether to show a green 'ON' 
%                                               for sweeps in use
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == set in plot_fitted_traces.m
%                   - 'PlotMarkFlag': whether to plot the way Mark 
%                                       wants plots to look
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'ShowSweepsFlag': whether to show sweep figures
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'JitterFlag': whether to introduce noise in parameters
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Grouping': a grouping vector used to group traces
%                   must be a vector
%                   default == []
%                   - 'CprWindow': current pulse response window in ms
%                   must be a numeric vector with 2 elements
%                   default == [0, 360] + timeToStabilize
%                   - 'IpscTime': start time of IPSC in ms
%                   must be a numeric scalar
%                   default == 1000 + timeToStabilize
%                   - 'IpscPeakWindow': window (ms) to look for IPSC peak
%                   must be a numeric vector
%                   default == [1000, 1300] + timeToStabilize
%                   - 'IpscrWindow': IPSC response window in ms
%                   must be a numeric vector with 2 elements
%                   default == [0, 8000] + timeToStabilize
%                   - 'OutFilePath': path to NEURON output file(s)
%                   must be a characeter vector, a string array 
%                       or a cell array of character arrays
%                   default == 'auto'
%                   - 'Tstop': simulation end time
%                   must be a numeric vector
%                   default == based on simMode, cprWindow(2) & ipscrWindow(2)
%                   - 'RealDataIpscr': recorded data to compare against
%                   must be a numeric array
%                       or a cell array of numeric arrays
%                   default == []
%                   - 'RealDataCpr': recorded data to compare against
%                   must be a numeric array
%                       or a cell array of numeric arrays
%                   default == []
%                   - 'HoldPotentialIpscr': holding potential for (mV)
%                                           IPSC response
%                   must be a numeric vector
%                   default == -70 mV
%                   - 'HoldPotentialCpr': holding potential for (mV)
%                                           current pulse response
%                   must be a numeric vector
%                   default == -70 mV
%                   - 'CurrentPulseAmplitudeIpscr': current pulse amplitude (nA)
%                   must be a numeric vector
%                   default == -0.050 nA
%                   - 'CurrentPulseAmplitudeCpr': current pulse amplitude (nA)
%                   must be a numeric vector
%                   default == -0.050 nA
%                   - 'GababAmpIpscr': GABA-B IPSC amplitude (uS)
%                   must be a numeric vector
%                   default == set of all 12 input waveforms
%                   - 'GababTriseIpscr': GABA-B IPSC rising time constant (ms)
%                   must be a numeric vector
%                   default == set of all 12 input waveforms
%                   - 'GababTfallFastIpscr': GABA-B IPSC falling phase 
%                                           fast component time constant (ms)
%                   must be a numeric vector
%                   default == set of all 12 input waveforms
%                   - 'GababTfallSlowIpscr': GABA-B IPSC falling phase 
%                                           slow component time constant (ms)
%                   must be a numeric vector
%                   default == set of all 12 input waveforms
%                   - 'GababWeightIpscr': GABA-B IPSC falling phase 
%                                           fast component weight
%                   must be a numeric vector
%                   default == set of all 12 input waveforms
%                   - 'HoldCurrentIpscr': custom holding current (nA)
%                                           for IPSC response
%                   must be a numeric vector
%                   default == 0 nA
%                   - 'HoldCurrentCpr': custom holding current (nA)
%                                       for current pulse response
%                   must be a numeric vector
%                   default == 0 nA
%                   - 'HoldCurrentNoiseIpscr': custom holding current noise (nA)
%                                               for IPSC response
%                   must be a numeric vector
%                   default == 0 nA
%                   - 'HoldCurrentNoiseCpr': custom holding current noise (nA)
%                                               for current pulse response
%                   must be a numeric vector
%                   default == 0 nA
%                   - 'RowConditions': row conditions for plotting
%                       Note: each row is assigned a different color
%                   must be a numeric 2D array
%                   default == RowConditionsIpscr or RowConditionsCpr
%                   - 'RowConditionsIpscr': row conditions for plotting
%                                               IPSC response
%                       Note: each row is assigned a different color
%                   must be a numeric 2D array
%                   default == transpose(1:nSweeps)
%                   - 'RowConditionsCpr': row conditions for plotting
%                                               current pulse response
%                       Note: each row is assigned a different color
%                   must be a numeric 2D array
%                   default == transpose(1:nSweeps)
%                   - 'FitWindowCpr': time window to fit (ms)
%                                       for current pulse response
%                   must be a numeric vector with 2 elements
%                   default == [100, 250] + timeToStabilize
%                   - 'FitWindowIpscr': time window to fit (ms)
%                                       for IPSC response
%                   must be a numeric vector with 2 elements
%                   default == [1000, 8000] + timeToStabilize
%                   - 'BaseWindowCpr': baseline window (ms)
%                                       for current pulse response
%                   must be a numeric vector with 2 elements
%                   default == [0, 100] + timeToStabilize
%                   - 'BaseWindowIpscr': baseline window (ms)
%                                       for IPSC response
%                   must be a numeric vector with 2 elements
%                   default == [0, 1000] + timeToStabilize
%                   - 'BaseNoiseCpr': baseline noise (mV)
%                                       for current pulse response
%                   must be a numeric vector
%                   default == set in compute_default_sweep_info.m
%                   - 'BaseNoiseIpscr': baseline noise (mV)
%                                       for IPSC response
%                   must be a numeric vector
%                   default == set in compute_default_sweep_info.m
%                   - 'SweepWeightsCpr': sweep weights 
%                                       for current pulse response
%                   must be a numeric vector
%                   default == set in compute_default_sweep_info.m
%                   - 'SweepWeightsIpscr': sweep weights
%                                       for IPSC response
%                   must be a numeric vector
%                   default == set in compute_default_sweep_info.m
%                   - 'LtsFeatureWeights': LTS feature weights for averaging
%                   must be empty or a numeric vector with length == nSweeps
%                   default == set in compute_lts_errors.m
%                   default == set in compute_lts_errors.m
%                   - 'MissedLtsError': a dimensionless error that penalizes 
%                                       a misprediction of the existence of LTS
%                   must be empty or a numeric vector with length == nSweeps
%                   default == set in compute_lts_errors.m
%                   - 'FalseLtsError': a dimensionless error that penalizes 
%                                       a misprediction of the absence of LTS
%                   must be empty or a numeric vector with length == nSweeps
%                   default == set in compute_lts_errors.m
%                   - 'Match2FeatureErrorRatio': ratio of LTS match error to 
%                                                   LTS feature error
%                   must be empty or a numeric vector with length == nSweeps
%                   default == set in compute_lts_errors.m
%                   - 'Lts2SweepErrorRatio': ratio of LTS error to sweep error
%                   must be empty or a numeric vector with length == nSweeps
%                   default == set in compute_single_neuron_errors.m
%
% Requires:
%       ~/m3ha/optimizer4gabab/singleneuron4compgabab.hoc
%       cd/argfun.m
%       cd/compute_combined_data.m
%       cd/compute_default_sweep_info.m
%       cd/compute_residuals.m
%       cd/compute_sampling_interval.m
%       cd/compute_single_neuron_errors.m
%       cd/create_simulation_output_filenames.m
%       cd/decide_on_colormap.m
%       cd/decide_on_filebases.m
%       cd/extract_columns.m
%       cd/extract_common_prefix.m
%       cd/extract_fileparts.m
%       cd/extract_substrings.m
%       cd/extract_subvectors.m
%       cd/find_window_endpoints.m
%       cd/force_matrix.m
%       cd/load_neuron_outputs.m
%       cd/load_params.m
%       cd/log_arraytext.m
%       cd/m3ha_extract_cell_name.m
%       cd/m3ha_neuron_create_simulation_params.m
%       cd/m3ha_neuron_create_TC_commands.m
%       cd/plot_fitted_traces.m
%       cd/m3ha_select_raw_traces.m
%       cd/parse_ipsc.m
%       cd/parse_lts.m
%       cd/parse_pulse_response.m
%       cd/run_neuron.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/test_var_difference.m
%
% Used by:    
%       cd/m3ha_fminsearch3.m
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_optimizer_4compgabab.m

% File History:
% 2014-04-XX - Created by Christine
% 2016-07-19 - changed ltsWindow
% 2016-07-20 - Added JITTER mode option
% 2016-07-20 - changed di
% 2016-07-20 - Analyze & plot LTS data
% 2016-10-04 - Changed the way NEURON is run: Added parfor and simCommands, etc.
% 2016-10-04 - Changed the method for execution to be from a here document 
%               written in the Matlab code directly
% 2016-10-05 - Moved outparams.lts2SweepErrorRatio to 
%               the calculation of total error (instead of total LTS error)
% 2016-10-05 - Added current pulse response
% 2016-10-05 - Changed the way error calculations are organized; 
%               added compute_and_compare_statistics
% 2016-10-06 - Reorganized code
% 2016-10-06 - Removed err_stats2_sim{bi} & err_stats2_real{bi}
% 2016-10-06 - Changed swpreg to fitreg
% 2016-10-06 - Renamed figure handles so that they are now all in a structure 
%               hFig that is passed to and from functions
% 2016-10-07 - outparams.currpulse(iSwp) is now already in nA
% 2016-10-07 - Added cprflag, findLtsFlag, ltsBurstStatsFlag, ltsErrorFlag
% 2016-10-14 - Updated outputs for m3ha_find_ipsc_peak & m3ha_find_lts
% 2016-10-14 - Fixed outparams.fitreg to fitreg inside parfor loop
% 2016-10-14 - Changed from root mean-squared error to mean-squared error
% 2017-01-14 - Added build() to simCommands
% 2017-01-14 - Added simMode to both build() and sim() of simCommands
% 2017-02-06 - Fixed fprintf of status to be %d from %s
% 2017-04-21 - Added outFolder to output file names
% 2017-04-24 - Renamed simulation_params -> neuronparams
% 2017-05-13 - Now gets outFolder from outparams
% 2017-05-15 - Changed outparams.prefix so that '_cpr' is already incorporated
% 2017-05-15 - Made simMode the last argument of build()
% 2017-05-16 - Added saveParamsFlag, saveStdOutFlag, 
%               saveSimCmdsFlag & saveSimOutFlag
% 2017-05-17 - Added saveLtsInfoFlag & saveLtsStatsFlag
% 2017-05-17 - Updated color groups to reflect pharm-g incr pairs
% 2017-05-17 - Now uses root-mean-squared errors instead of mean-squared errors
% 2017-05-17 - Weighting of errors is now normalized
% 2017-05-17 - Made 'a' column vectors
% 2017-05-17 - Now computes average sweep error instead of total sweep error
% 2017-05-17 - Made all errors dimensionless by normalizing by the real data value
% 2017-05-17 - Normalize by initial error instead, rmse now has units of mV again
% 2017-05-17 - Added outparams.isLtsError
% 2017-05-17 - Fixed the case where simulation correctly predicted 
%               the non-existence of LTS (should be 0 error)
% 2017-05-19 - Moved update_sweeps_figures() here from m3ha_optimizer_4compgabab.m
% 2017-05-19 - Separated update_sweeps_figures() into 
%               decide_on_xlimits(), plot_overlapped_traces(),
%               plot_conductance_traces(), plot_current_traces() 
%               && plot_fitted_traces()
% 2017-05-22 - Changed line width and indentation
% 2017-05-23 - Removed modeselected from outparams and replaced with updated outparams.runMode
% 2017-07-27 - Now extracts all NEURON parameters to workspace
% 2017-07-28 - Added commands for calling adjust_globalpas, adjust_leak, 
%               adjust_IT, adjust_Ih, adjust_IA, adjust_IKir, adjust_INaP
% 2017-07-28 - Changed precision of parameters from %6.3f to %g
% 2017-07-28 - Made normalization to starting error an option (outparams.normalize2InitErrFlag)
% 2017-07-28 - Made sweep error dimensionless by dividing by the holding potential
% 2017-07-28 - Made sweep error dimensionless by dividing by the the maximum noise of sweep
% 2017-07-29 - Now normalizes LTS amp error by maximum noise instead
% 2017-07-29 - Now normalizes LTS time error by ioffset instead
% 2017-07-29 - Now normalizes LTS slope error by 
%               slope*(2*maximum noise/peakprom + 2*ioffset/peakwidth) instead
% 2017-07-29 - Now normalizes LTS time error by peakwidth instead
% 2017-08-11 - Added outparams.colmode and load across trials NEURON parameters
%               from outparams.neuronparamsBest for each cell
% 2017-08-12 - Now loads across cells NEURON parameters from 
%               outparams.neuronparamsAcrossCells when fitting across trials
% 2017-08-21 - Fixed ncg = size(colorMap, 2) -> size(colorMap, 1);
%               colorMap{cgn} -> colorMap(cgn, :)
% 2017-08-29 - Removed dend0
% 2017-11-09 - Replaced saveas with save_all_figtypes
% 2017-11-09 - Now uses plot_bar.m instead of errorbar
% 2017-12-21 - Added some more output ~s to m3ha_find_lts.m 
% 2018-01-24 - Added isdeployed
% 2018-03-02 - Added runNeuronCommand, moduleLoadCommands, moduleLoadCommandsHpc
% 2018-03-02 - Added onHpcFlag as an optional argument
% 2018-03-07 - Added outparams.showStatisticsFlag
% 2018-03-09 - Fixed the case of plotting scatter plots when the standard
%               deviation is zero
% 2018-05-21 - Changed voltage uncertainty from outparams.maxNoise to 
%               outparams.baseNoise or outparams.baseNoiseCpr
% 2018-05-21 - Now uses nanmean to compute rmse() to account for 
%               NaN in cpr response values
% 2018-05-21 - Changed color groups for current pulse response to 3 groups
% 2018-07-09 - Now uses compute_rms_error
% 2018-07-12 - BT - Added plotting residuals
% 2018-08-09 - Now does not divide by baseNoise, but use baseNoise
%               to modify sweep weights inst
% 2018-08-09 - Now plots sweep weights
% 2018-08-10 - Changed fitreg to fitWindow
% 2018-08-26 - Now links x and y axes for subplots
% 2018-10-01 - Added holdcurrentNoise
% 2018-10-15 - Updated usage of TC3
% 2018-10-16 - Now uses save_params.m
% 2018-10-19 - Made 'RealData' as an optional parameter
% 2018-11-16 - Separated 'RealDataCpr' from 'RealDataIpscr'
% 2019-01-08 - Reorganized code so that one can run a single simulation easily
%                   from a set of parameters
% 2019-01-09 - Added 'GenerateDataFlag' as an optional parameter
% 2019-01-14 - Added 'BootstrapCprFlag' as an optional parameter
% 2019-05-08 - Updated usage of plot_bar.m
% 2019-10-13 - Updated simulated data column numbers
% 2019-11-15 - Added 'IpscTime' as an optional parameter
% 2019-11-15 - Added 'IpscPeakWindow' as an optional parameter
% 2019-11-15 - Added 'FileBases' as an optional parameter
% 2019-11-17 - Added 'PlotConductanceFlag' as an optional parameter
% 2019-11-17 - Added 'PlotCurrentFlag' as an optional parameter
% 2019-12-03 - Added 'FileNames' as an optional parameter
% 2019-12-04 - If 'FileNames' is passed in, use the same GABA-B IPSC parameters
%               as those used in dynamic clamp
% 2019-12-17 - Fixed GABAB parameter bug
% 2019-12-18 - Added 'Match2FeatureErrorRatio' as an optional parameter
% 2019-12-19 - Now distinguishes buildMode from simMode
% 2019-12-19 - Removed y axis link for current pulse response plots
% 2019-12-21 - Added 'PlotAllFlag' as an optional argument and change
%                   default to not plot anything
% 2019-12-23 - Added 'SaveImportLogFlag' as an optional parameter


%% Hard-coded parameters
validBuildModes = {'active', 'passive'};
validSimModes = {'active', 'passive'};
hocFile = 'singleneuron4compgabab.hoc';
maxRowsWithOneOnly = 8;
verbose = false;
lineWidthParallel = 1;
figTypes = 'png';
importedSuffix = 'imported_files';

% The following must be consistent with both dclampDataExtractor.m & ...
%   singleneuron4compgabab.hoc
cprWinOrig = [0, 360];          % current pulse response window (ms), original
ipscTimeOrig = 1000;            % time of IPSC application (ms), original
ipscrWinOrig = [0, 8000];       % IPSC response window (ms), original
ipscpWinOrig = [1000, 1300];    % window (ms) in which IPSC reaches peak 
                                %   amplitude , original
                                %   Based on observation, IPSCs are 
                                %   not influenced by LTSs before 1300 ms

% The following must be consistent with singleneuron4compgabab.hoc
timeToStabilize = 2000;         % padded time (ms) to make sure initial value 
                                %   of simulations are stabilized

% Default time windows to fit
%   Note: must be consistent with singleneuronfitting.m
baseWinCprOrig      = [0, 100];     % window in which the current pulse response base line woule lie (ms), original
fitWinCprOrig    = [100, 250];      % time window to fit for current pulse response (ms), original
baseWinIpscrOrig    = [0, 1000];    % window in which the IPSC response base line woule lie (ms), original
fitWinIpscrOrig  = [1000, 8000];    % time window to fit for IPSC response (ms), original

% For plotting
figNumberConductance = 102;
figNumberCurrent = 103;
figNumberIndividual = 104;
figNumberResiduals = 105;
figNumberOverlapped = 106;

%% Column numbers for recorded data
%   Note: Must be consistent with m3ha_resave_sweeps.m
TIME_COL_REC = 1;
VOLT_COL_REC = 2;
CURR_COL_REC = 3;
COND_COL_REC = 4;

%% Column numbers for simulated data
%   Note: Must be consistent with singleneuron4compgabab.hoc
TIME_COL_SIM = 1;
VOLT_COL_SIM = 2;
DEND1_COL_SIM = 3;
DEND2_COL_SIM = 4;
IDCLAMP_COL_SIM = 5;
GGABAB_COL_SIM = 6;
iCP_COL_SIM = 7;
IEXT_COL_SIM = 8;
ICA_COL_SIM = 9;
ITM_COL_SIM = 10;
ITMINF_COL_SIM = 11;
ITH_COL_SIM = 12;
ITHINF_COL_SIM = 13;
IH_COL_SIM = 14;
IHM_COL_SIM = 15;
IKA_COL_SIM = 16;
IAM1_COL_SIM = 17;
IAH1_COL_SIM = 18;
IAM2_COL_SIM = 19;
IAH2_COL_SIM = 20;
IKKIR_COL_SIM = 21;
IKIRM_COL_SIM = 22;
INAPNA_COL_SIM = 23;
INAPM_COL_SIM = 24;
INAPH_COL_SIM = 25;

%% Default values for optional arguments
hFigDefault = '';               % no prior hFig structure by default
buildModeDefault = 'active';    % insert active channels by default
simModeDefault = 'active';      % simulate active responses by default
columnModeDefault = [];         % set in m3ha_select_raw_traces.m
rowModeDefault = [];            % set in m3ha_select_raw_traces.m
attemptNumberDefault = [];      % set in m3ha_select_raw_traces.m
dataModeDefault = [];           % set in m3ha_select_raw_traces.m
nSweepsDefault = [];            % set later
fileNamesDefault = {};          % none provided by default
fileNamesCprDefault = {};       % none provided by default
fileNamesIpscrDefault = {};     % none provided by default
fileBasesDefault = {};          % set later
outFolderDefault = '' ;         % set later
prefixDefault = '';             % set later
debugFlagDefault = false;       % not in debug mode by default
customHoldCurrentFlagDefault = 0; % don't use custom hold current by default
onHpcFlagDefault = false;       % not on a high performance computing
                                %   server by default
generateDataFlagDefault = false;% not generating surrogate data by default
averageCprFlagDefault = false;  % don't average current pulse responses 
                                %   according to vHold by default
bootstrapCprFlagDefault = false;% don't bootstrap-average current pulse  
                                %   responses according to vHold by default
normalize2InitErrFlagDefault = false;
saveImportLogFlagDefault = true;    % save imported file log by default
saveParamsFlagDefault = true;   % save simulation parameters by default
saveSimCmdsFlagDefault = true;  % save simulation commands by default
saveStdOutFlagDefault = false;  % save standard outputs only if error by default
saveSimOutFlagDefault = true;   % save simulation outputs by default
saveLtsInfoFlagDefault = true;  % save LTS info by default
saveLtsStatsFlagDefault = true; % save LTS statistics by default
plotAllFlagDefault = false;         % plot nothing by default
plotIndividualFlagDefault = [];     % set later
plotResidualsFlagDefault = [];      % set later
plotOverlappedFlagDefault = [];     % set later
plotConductanceFlagDefault = [];    % set later
plotCurrentFlagDefault = [];        % set later
plotIpeakFlagDefault = [];          % set later
plotLtsFlagDefault = [];            % set later
plotStatisticsFlagDefault = [];     % set later
plotSwpWeightsFlagDefault = 'auto'; % set in m3ha_plot_inidividual_traces.m
plotMarkFlagDefault = true;         % the way Mark wants plots to look
showSweepsFlagDefault = true;       % whether to show sweep figures
jitterFlagDefault = false;          % no jitter by default
groupingDefault = [];               % no grouping by default
cprWindowDefault = cprWinOrig + timeToStabilize;
ipscTimeDefault = ipscTimeOrig + timeToStabilize;
ipscPeakWindowDefault = ipscpWinOrig + timeToStabilize;
ipscrWindowDefault = ipscrWinOrig + timeToStabilize;
outFilePathDefault = 'auto';    % set later
tstopDefault = [];              % set later
realDataIpscrDefault = [];      % no data to compare against by default
realDataCprDefault = [];        % no data to compare against by default
holdPotentialIpscrDefault = -70;% (mV)
holdPotentialCprDefault = -70;  % (mV)
currentPulseAmplitudeIpscrDefault = -0.050;  % (nA)
currentPulseAmplitudeCprDefault = -0.050;  % (nA)
gababAmpIpscrDefault = [];           % set later
gababTriseIpscrDefault = [];         % set later
gababTfallFastIpscrDefault = [];     % set later
gababTfallSlowIpscrDefault = [];     % set later
gababWeightIpscrDefault = [];        % set later
holdCurrentIpscrDefault = 0;    % (nA)
holdCurrentCprDefault = 0;      % (nA)
holdCurrentNoiseIpscrDefault = 0;% (nA)
holdCurrentNoiseCprDefault = 0; % (nA)
rowConditionsDefault = [];          % set later
rowConditionsIpscrDefault = [];     % set later
rowConditionsCprDefault = [];       % set later
fitWindowCprDefault = fitWinCprOrig + timeToStabilize;
fitWindowIpscrDefault = fitWinIpscrOrig + timeToStabilize;
baseWindowCprDefault = baseWinCprOrig + timeToStabilize;
baseWindowIpscrDefault = baseWinIpscrOrig + timeToStabilize;
baseNoiseCprDefault = [];       % set later
baseNoiseIpscrDefault = [];     % set later
sweepWeightsCprDefault = [];    % set later
sweepWeightsIpscrDefault = [];  % set later
ltsFeatureWeightsDefault = [];  % set in compute_lts_errors.m
missedLtsErrorDefault = [];     % set in compute_lts_errors.m
falseLtsErrorDefault = [];      % set in compute_lts_errors.m
match2FeatureErrorRatioDefault = []; % set in compute_lts_errors.m
lts2SweepErrorRatioDefault = [];% set in compute_single_neuron_errors.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'neuronParamsTableOrFile', ...
    @(x) validateattributes(x, {'table', 'char', 'string', 'cell'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'HFig', hFigDefault);
addParameter(iP, 'BuildMode', buildModeDefault, ...
    @(x) any(validatestring(x, validBuildModes)));
addParameter(iP, 'SimMode', simModeDefault, ...
    @(x) any(validatestring(x, validSimModes)));
addParameter(iP, 'ColumnMode', columnModeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer'}));
addParameter(iP, 'RowMode', rowModeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer'}));
addParameter(iP, 'AttemptNumber', attemptNumberDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer'}));
addParameter(iP, 'DataMode', dataModeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'integer'}));
addParameter(iP, 'NSweeps', nSweepsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'FileNames', fileNamesDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['FileNames must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'FileNamesCpr', fileNamesCprDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['FileNamesCpr must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'FileNamesIpscr', fileNamesIpscrDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['FileNamesIpscr must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'FileBases', fileBasesDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['FileBases must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'DebugFlag', debugFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'CustomHoldCurrentFlag', customHoldCurrentFlagDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'binary'}));
addParameter(iP, 'OnHpcFlag', onHpcFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'GenerateDataFlag', generateDataFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'AverageCprFlag', averageCprFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'BootstrapCprFlag', bootstrapCprFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Normalize2InitErrFlag', normalize2InitErrFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveImportLogFlag', saveImportLogFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveParamsFlag', saveParamsFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveSimCmdsFlag', saveSimCmdsFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveStdOutFlag', saveStdOutFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveSimOutFlag', saveSimOutFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveLtsInfoFlag', saveLtsInfoFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveLtsStatsFlag', saveLtsStatsFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotAllFlag', plotAllFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotIndividualFlag', plotIndividualFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotResidualsFlag', plotResidualsFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotOverlappedFlag', plotOverlappedFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotConductanceFlag', plotConductanceFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotCurrentFlag', plotCurrentFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotIpeakFlag', plotIpeakFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotLtsFlag', plotLtsFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotStatisticsFlag', plotStatisticsFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotSwpWeightsFlag', plotSwpWeightsFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotMarkFlag', plotMarkFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ShowSweepsFlag', showSweepsFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'JitterFlag', jitterFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Grouping', groupingDefault);
addParameter(iP, 'CprWindow', cprWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
addParameter(iP, 'IpscTime', ipscTimeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'IpscPeakWindow', ipscPeakWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'IpscrWindow', ipscrWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
addParameter(iP, 'OutFilePath', outFilePathDefault, ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));
addParameter(iP, 'Tstop', tstopDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'RealDataIpscr', realDataIpscrDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['RealDataIpscr must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'RealDataCpr', realDataCprDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['RealDataCpr must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'HoldPotentialIpscr', holdPotentialIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'HoldPotentialCpr', holdPotentialCprDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'CurrentPulseAmplitudeIpscr', ...
                currentPulseAmplitudeIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'CurrentPulseAmplitudeCpr', ...
                currentPulseAmplitudeCprDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'GababAmpIpscr', gababAmpIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'GababTriseIpscr', gababTriseIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'GababTfallFastIpscr', gababTfallFastIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'GababTfallSlowIpscr', gababTfallSlowIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'GababWeightIpscr', gababWeightIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'HoldCurrentIpscr', holdCurrentIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'HoldCurrentCpr', holdCurrentCprDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'HoldCurrentNoiseIpscr', holdCurrentNoiseIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'HoldCurrentNoiseCpr', holdCurrentNoiseCprDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'RowConditions', rowConditionsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'RowConditionsIpscr', rowConditionsIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'RowConditionsCpr', rowConditionsCprDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'FitWindowCpr', fitWindowCprDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
addParameter(iP, 'FitWindowIpscr', fitWindowIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
addParameter(iP, 'BaseWindowCpr', baseWindowCprDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
addParameter(iP, 'BaseWindowIpscr', baseWindowIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
addParameter(iP, 'BaseNoiseCpr', baseNoiseCprDefault, ...
    @(x) assert(isnumericvector(x), 'BaseNoiseCpr must be a numeric vector!'));
addParameter(iP, 'BaseNoiseIpscr', baseNoiseIpscrDefault, ...
    @(x) assert(isnumericvector(x), 'BaseNoiseIpscr must be a numeric vector!'));
addParameter(iP, 'SweepWeightsCpr', sweepWeightsCprDefault, ...
    @(x) assert(isnumericvector(x), 'SweepWeightsCpr must be a numeric vector!'));
addParameter(iP, 'SweepWeightsIpscr', sweepWeightsIpscrDefault, ...
    @(x) assert(isnumericvector(x), 'SweepWeightsIpscr must be a numeric vector!'));
addParameter(iP, 'LtsFeatureWeights', ltsFeatureWeightsDefault, ...
    @(x) assert(isnumericvector(x), 'LtsFeatureWeights must be a numeric vector!'));
addParameter(iP, 'MissedLtsError', missedLtsErrorDefault, ...
    @(x) assert(isnumericvector(x), 'MissedLtsError must be a numeric vector!'));
addParameter(iP, 'FalseLtsError', falseLtsErrorDefault, ...
    @(x) assert(isnumericvector(x), 'FalseLtsError must be a numeric vector!'));
addParameter(iP, 'Match2FeatureErrorRatio', match2FeatureErrorRatioDefault, ...
    @(x) assert(isnumericvector(x), 'Match2FeatureErrorRatio must be a numeric vector!'));
addParameter(iP, 'Lts2SweepErrorRatio', lts2SweepErrorRatioDefault, ...
    @(x) assert(isnumericvector(x), 'Lts2SweepErrorRatio must be a numeric vector!'));

% Read from the Input Parser
parse(iP, neuronParamsTableOrFile, varargin{:});
hFig = iP.Results.HFig;
buildMode = validatestring(iP.Results.BuildMode, validBuildModes);
simMode = validatestring(iP.Results.SimMode, validSimModes);
columnMode = iP.Results.ColumnMode;
rowMode = iP.Results.RowMode;
attemptNumber = iP.Results.AttemptNumber;
dataMode = iP.Results.DataMode;
nSweepsUser = iP.Results.NSweeps;
fileNames = iP.Results.FileNames;
fileNamesCpr = iP.Results.FileNamesCpr;
fileNamesIpscr = iP.Results.FileNamesIpscr;
fileBases = iP.Results.FileBases;
outFolder = iP.Results.OutFolder;
prefix = iP.Results.Prefix;
debugFlag = iP.Results.DebugFlag;
customHoldCurrentFlag = iP.Results.CustomHoldCurrentFlag;
onHpcFlag = iP.Results.OnHpcFlag;
generateDataFlag = iP.Results.GenerateDataFlag;
averageCprFlag = iP.Results.AverageCprFlag;
bootstrapCprFlag = iP.Results.BootstrapCprFlag;
normalize2InitErrFlag = iP.Results.Normalize2InitErrFlag;
saveImportLogFlag = iP.Results.SaveImportLogFlag;
saveParamsFlag = iP.Results.SaveParamsFlag;
saveSimCmdsFlag = iP.Results.SaveSimCmdsFlag;
saveStdOutFlag = iP.Results.SaveStdOutFlag;
saveSimOutFlag = iP.Results.SaveSimOutFlag;
saveLtsInfoFlag = iP.Results.SaveLtsInfoFlag;
saveLtsStatsFlag = iP.Results.SaveLtsStatsFlag;
plotAllFlag = iP.Results.PlotAllFlag;
plotIndividualFlag = iP.Results.PlotIndividualFlag;
plotResidualsFlag = iP.Results.PlotResidualsFlag;
plotOverlappedFlag = iP.Results.PlotOverlappedFlag;
plotConductanceFlag = iP.Results.PlotConductanceFlag;
plotCurrentFlag = iP.Results.PlotCurrentFlag;
plotIpeakFlag = iP.Results.PlotIpeakFlag;
plotLtsFlag = iP.Results.PlotLtsFlag;
plotStatisticsFlag = iP.Results.PlotStatisticsFlag;
plotSwpWeightsFlag = iP.Results.PlotSwpWeightsFlag;
plotMarkFlag = iP.Results.PlotMarkFlag;
showSweepsFlag = iP.Results.ShowSweepsFlag;
jitterFlag = iP.Results.JitterFlag;
grouping = iP.Results.Grouping;
cprWindow = iP.Results.CprWindow;
ipscTime = iP.Results.IpscTime;
ipscPeakWindow = iP.Results.IpscPeakWindow;
ipscrWindow = iP.Results.IpscrWindow;
outFilePath = iP.Results.OutFilePath;
tstop = iP.Results.Tstop;
realDataIpscr = iP.Results.RealDataIpscr;
realDataCpr = iP.Results.RealDataCpr;
holdPotentialIpscr = iP.Results.HoldPotentialIpscr;
holdPotentialCpr = iP.Results.HoldPotentialCpr;
currentPulseAmplitudeIpscr = iP.Results.CurrentPulseAmplitudeIpscr;
currentPulseAmplitudeCpr = iP.Results.CurrentPulseAmplitudeCpr;
gababAmpIpscr = iP.Results.GababAmpIpscr;
gababTriseIpscr = iP.Results.GababTriseIpscr;
gababTfallFastIpscr = iP.Results.GababTfallFastIpscr;
gababTfallSlowIpscr = iP.Results.GababTfallSlowIpscr;
gababWeightIpscr = iP.Results.GababWeightIpscr;
holdCurrentIpscr = iP.Results.HoldCurrentIpscr;
holdCurrentCpr = iP.Results.HoldCurrentCpr;
holdCurrentNoiseIpscr = iP.Results.HoldCurrentNoiseIpscr;
holdCurrentNoiseCpr = iP.Results.HoldCurrentNoiseCpr;
rowConditions = iP.Results.RowConditions;
rowConditionsIpscr = iP.Results.RowConditionsIpscr;
rowConditionsCpr = iP.Results.RowConditionsCpr;
fitWindowCpr = iP.Results.FitWindowCpr;
fitWindowIpscr = iP.Results.FitWindowIpscr;
baseWindowCpr = iP.Results.BaseWindowCpr;
baseWindowIpscr = iP.Results.BaseWindowIpscr;
baseNoiseCpr = iP.Results.BaseNoiseCpr;
baseNoiseIpscr = iP.Results.BaseNoiseIpscr;
sweepWeightsCpr = iP.Results.SweepWeightsCpr;
sweepWeightsIpscr = iP.Results.SweepWeightsIpscr;
ltsFeatureWeights = iP.Results.LtsFeatureWeights;
missedLtsError = iP.Results.MissedLtsError;
falseLtsError = iP.Results.FalseLtsError;
match2FeatureErrorRatio = iP.Results.Match2FeatureErrorRatio;
lts2SweepErrorRatio = iP.Results.Lts2SweepErrorRatio;

%% Preparation
% Initialize outputs
errorStruct = struct;
hFig = struct;
simData = [];

% Create an experiment identifier
expStr = prefix;

% Parse first argument
if istext(neuronParamsTableOrFile)
    % Read in the table
    neuronParamsTable = load_params(neuronParamsTableOrFile);

    % Look for a cell name
    cellName = m3ha_extract_cell_name(neuronParamsTableOrFile, ...
                                        'ForceSingleOutput', true);

    % Set a default expStr
    if isempty(expStr)
        % Try extracting a common prefix
        paramFileBases = extract_fileparts(neuronParamsTableOrFile, 'base');
        commonPrefix = extract_common_prefix(paramFileBases);

        % Decide on the prefix
        if ~isempty(commonPrefix)
            expStrPrefix = commonPrefix;
        elseif ~isempty(cellName)
            expStrPrefix = cellName;
        else
            expStrPrefix = 'unnamed';
        end

        % Decide on the suffix
        if strcmpi(simMode, 'passive')
            expStr = strcat(expStrPrefix, '_cpr');
        else
            expStr = strcat(expStrPrefix, '_ipscr');
        end
    end
else
    % The first argument is(are) the parameter table(s)
    neuronParamsTable = neuronParamsTableOrFile;
    
    % Cell name is not provided yet at this stage
    cellName = '';
end

% Choose data files to compare against if cell name provided
if isempty(fileNames)
    if strcmpi(simMode, 'passive') && ~isempty(fileNamesCpr)
        fileNames = fileNamesCpr;
    elseif strcmpi(simMode, 'active') && ~isempty(fileNamesIpscr)
        fileNames = fileNamesIpscr;
    elseif ~isempty(cellName)
        % Select data files
        [fileNames, rowConditions] = ...
            m3ha_select_raw_traces(cellName, 'ColumnMode', columnMode, ...
                            'RowMode', rowMode, 'AttemptNumber', attemptNumber, ...
                            'DataMode', dataMode);
    end
elseif isempty(cellName)
    % Look for a cell name
    cellName = m3ha_extract_cell_name(fileNames, 'ForceSingleOutput', true);
end

% Create a default output folder
if isempty(outFolder)
    if ~isempty(expStr)
        outFolder = expStr;
    else
        outFolder = pwd;
    end
end

% Make sure the output folder exists
check_dir(outFolder);

% Create output paths
importedPath = fullfile(outFolder, [expStr, '_', importedSuffix, '.txt']);

% Update flags
if verbose
    fprintf('Setting default flags for %s ...\n', fileBase);
end
[plotIndividualFlag, plotResidualsFlag, plotOverlappedFlag, ...
plotConductanceFlag, plotCurrentFlag, plotIpeakFlag, ...
plotLtsFlag, plotStatisticsFlag] = ...
    argfun(@(x) set_default_flag(x, plotAllFlag), ...
                plotIndividualFlag, plotResidualsFlag, plotOverlappedFlag, ...
                plotConductanceFlag, plotCurrentFlag, plotIpeakFlag, ...
                plotLtsFlag, plotStatisticsFlag);

% Decide on simulation-mode-dependent variables
if strcmpi(simMode, 'passive')
    % Import data to compare if not provided but file names provided
    if isempty(realDataCpr) && ~isempty(fileNames)
        % Log names of files imported
        if saveImportLogFlag
            log_arraytext(importedPath, fileNames);
        end

        % Import data
        [realDataCpr, sweepInfoCpr] = ...
            m3ha_import_raw_traces(fileNames, 'ImportMode', 'passive', ...
                        'Verbose', verbose, 'OutFolder', outFolder);

        % Update file names
        fileNames = sweepInfoCpr.fileNames;

        % Update row conditions if needed
        if numel(fileNames) <= 3
            rowConditions = transpose(1:3);
        end

        % If not generating data, don't set hold current noise
        if ~generateDataFlag
            nSweepsCpr = height(sweepInfoCpr);
            sweepInfoCpr{:, 'holdCurrentNoise'} = zeros(nSweepsCpr, 1);
        end

        % Append suffices to variable names
        sweepInfoCpr.Properties.VariableNames = ...
            strcat(sweepInfoCpr.Properties.VariableNames, 'Cpr');

        % Convert sweep info tables to scalar structures
        sweepInfoCprStruct = table2struct(sweepInfoCpr, 'ToScalar', true);

        % Extract from structure
        currentPulseAmplitudeCpr = sweepInfoCprStruct.currentPulseAmplitudeCpr;
        holdPotentialCpr = sweepInfoCprStruct.holdPotentialCpr;
        holdCurrentCpr = sweepInfoCprStruct.holdCurrentCpr;
        holdCurrentNoiseCpr = sweepInfoCprStruct.holdCurrentNoiseCpr;
        baseNoiseCpr = sweepInfoCprStruct.baseNoiseCpr;
        sweepWeightsCpr = sweepInfoCprStruct.sweepWeightsCpr;
    end

    % Read passive fit data and params
    realData = realDataCpr;
    currentPulseAmplitude = currentPulseAmplitudeCpr;
    holdPotential = holdPotentialCpr;
    holdCurrent = holdCurrentCpr;
    holdCurrentNoise = holdCurrentNoiseCpr;
    fitWindow = fitWindowCpr;
    baseWindow = baseWindowCpr;
    baseNoise = baseNoiseCpr;
    sweepWeights = sweepWeightsCpr;
    errorMode = 'SweepOnly';
elseif strcmpi(simMode, 'active')
    % Import data to compare if not provided but file names provided
    if isempty(realDataIpscr) && ~isempty(fileNames)
        % Log names of files imported
        log_arraytext(importedPath, fileNames);

        % Import data
        [realDataIpscr, sweepInfoIpscr] = ...
            m3ha_import_raw_traces(fileNames, 'ImportMode', 'active', ...
                        'Verbose', verbose, 'OutFolder', outFolder);

        % Update file names
        fileNames = sweepInfoIpscr.fileNames;

        % If not generating data, don't set hold current noise
        if ~generateDataFlag
            nSweepsIpscr = height(sweepInfoIpscr);
            sweepInfoIpscr{:, 'holdCurrentNoise'} = zeros(nSweepsIpscr, 1);
        end

        % Append suffices to variable names
        sweepInfoIpscr.Properties.VariableNames = ...
            strcat(sweepInfoIpscr.Properties.VariableNames, 'Ipscr');

        % Convert sweep info tables to scalar structures
        sweepInfoIpscrStruct = table2struct(sweepInfoIpscr, 'ToScalar', true);

        % Extract from structure
        currentPulseAmplitudeIpscr = ...
            sweepInfoIpscrStruct.currentPulseAmplitudeIpscr;
        holdPotentialIpscr = sweepInfoIpscrStruct.holdPotentialIpscr;
        holdCurrentIpscr = sweepInfoIpscrStruct.holdCurrentIpscr;
        holdCurrentNoiseIpscr = sweepInfoIpscrStruct.holdCurrentNoiseIpscr;
        baseNoiseIpscr = sweepInfoIpscrStruct.baseNoiseIpscr;
        sweepWeightsIpscr = sweepInfoIpscrStruct.sweepWeightsIpscr;
        gababAmpIpscr = sweepInfoIpscrStruct.gababAmpIpscr;
        gababTriseIpscr = sweepInfoIpscrStruct.gababTriseIpscr;
        gababTfallFastIpscr = sweepInfoIpscrStruct.gababTfallFastIpscr;
        gababTfallSlowIpscr = sweepInfoIpscrStruct.gababTfallSlowIpscr;
        gababWeightIpscr = sweepInfoIpscrStruct.gababWeightIpscr;
    end

    % Read active fit data and params
    realData = realDataIpscr;
    currentPulseAmplitude = currentPulseAmplitudeIpscr;
    holdPotential = holdPotentialIpscr;
    holdCurrent = holdCurrentIpscr;
    holdCurrentNoise = holdCurrentNoiseIpscr;
    fitWindow = fitWindowIpscr;
    baseWindow = baseWindowIpscr;
    baseNoise = baseNoiseIpscr;
    sweepWeights = sweepWeightsIpscr;
    errorMode = 'Sweep&LTS';
end

% Decide on row conditions
if isempty(rowConditions)
    if strcmpi(simMode, 'passive')
        rowConditions = rowConditionsCpr;
    elseif strcmpi(simMode, 'active')
        rowConditions = rowConditionsIpscr;
    end
end

% Decide whether to plot anything
if plotOverlappedFlag || plotConductanceFlag || ...
        plotCurrentFlag || plotIndividualFlag || plotResidualsFlag
    plotFlag = true;
else
    plotFlag = false;
end

% Decide on the number of sweeps to run and compare
nSweeps = decide_on_nSweeps(realData, nSweepsUser);

% Create file bases if not provided
if isempty(fileBases)
    fileBases = decide_on_filebases(fileNames, nSweeps);
end

% Decide on expStr if not provided
if isempty(expStr)
    % Try extracting the common prefix from the file bases
    expStr = extract_common_prefix(fileBases);

    % If still empty, use the output folder base name
    if isempty(expStr)
        expStr = extract_fileparts(outFolder, 'dirbase');
    end
end

% Decide on x-axis limits for plotting
if plotFlag
    xLimits = decide_on_xlimits(fitWindow, baseWindow, simMode, plotMarkFlag);
end

% Set figure visibility status
if showSweepsFlag
    visibleStatus = 'on';
else
    visibleStatus = 'off';
end

% Decide on rowConditions and nRows
[rowConditions, nRows] = ...
    decide_on_row_conditions(rowConditions, nSweeps, maxRowsWithOneOnly);

% Decide on the color map for individual and residual plots
colorMapIndividual = decide_on_colormap('r', nRows);

% Decide on the colors for parallel plots
colorMapParallel = decide_on_colormap([], 4);
if nSweeps > nRows
    nColumns = ceil(nSweeps / nRows);
    nSlots = nColumns * nRows;
    colorMapParallel = reshape(repmat(reshape(colorMapParallel, 1, []), ...
                        nColumns, 1), nSlots, 3);
end

% Print to standard output
if verbose
    fprintf('Preparing to run simulations for %s ...\n', expStr);
end

% Create output file paths if not provided
if strcmpi(outFilePath, 'auto')
    outFilePath = create_simulation_output_filenames(nSweeps, ...
                        'OutFolder', outFolder, 'Prefix', expStr, ...
                        'Suffix', fileBases);
end

% Create a table of simulation parameters
simParamsTable = m3ha_neuron_create_simulation_params(neuronParamsTable, ...
                        'Prefix', expStr, 'OutFolder', outFolder, ...
                        'SaveParamsFlag', saveParamsFlag, ...
                        'JitterFlag', jitterFlag, 'NSims', nSweeps, ...
                        'CprWindow', cprWindow, 'IpscrWindow', ipscrWindow, ...
                        'BuildMode', buildMode, 'SimMode', simMode, ...
                        'OutFilePath', outFilePath, 'Tstop', tstop, ...
                        'HoldPotential', holdPotential, ...
                        'CurrentPulseAmplitude', currentPulseAmplitude, ...
                        'GababAmp', gababAmpIpscr, ...
                        'GababTrise', gababTriseIpscr, ...
                        'GababTfallFast', gababTfallFastIpscr, ...
                        'GababTfallSlow', gababTfallSlowIpscr, ...
                        'GababWeight', gababWeightIpscr, ...
                        'CustomHoldCurrentFlag', customHoldCurrentFlag, ...
                        'HoldCurrent', holdCurrent, ...
                        'HoldCurrentNoise', holdCurrentNoise);

% Create simulation commands to be read by NEURON
simCommands = m3ha_neuron_create_TC_commands(simParamsTable, ...
                        'Prefix', expStr, 'OutFolder', outFolder, ...
                        'SaveSimCmdsFlag', saveSimCmdsFlag);

% Extract vectors from recorded data
%   Note: these will be empty if realData not provided
[tVecs, vVecsRec, iVecsRec, gVecsRec] = ...
    extract_columns(realData, [TIME_COL_REC, VOLT_COL_REC, ...
                                CURR_COL_REC, COND_COL_REC]);

% Compute baseline noise and sweep weights if not provided
[~, ~, baseNoise, sweepWeights] = ...
    compute_default_sweep_info(tVecs, vVecsRec, ...
            'BaseWindow', baseWindow, 'BaseNoise', baseNoise, ...
            'SweepWeights', sweepWeights);

%% Run NEURON
% Print to standard output
if verbose
    fprintf('Running simulations for %s ...\n', expStr);
end

% Run NEURON with the hocfile and attached simulation commands
output = run_neuron(hocFile, 'SimCommands', simCommands, ...
                    'Prefix', expStr, 'OutFolder', outFolder, ...
                    'DebugFlag', debugFlag, 'OnHpcFlag', onHpcFlag, ...
                    'SaveStdOutFlag', saveStdOutFlag);

% Check if there are errors
if any(output.hasError)
    fprintf('Simulations ran into error!\n');
    return
end
                
%% TODO: Break the function here into two

%% Analyze results
% Print to standard output
if verbose
    fprintf('Extracting simulation results for %s ...\n', expStr);
end

% Create an experiment identifier for title
expStrForTitle = strrep(expStr, '_', '\_');

% Load .out files created by NEURON
% If recorded data provided (tVecs not empty at this point),
%   interpolate simulated data to match the time points of recorded data
% Note: This is necessary because CVODE (variable time step method) 
%       is applied in NEURON
simData = load_neuron_outputs('FileNames', outFilePath, 'tVecs', tVecs, ...
                                'RemoveAfterLoad', ~saveSimOutFlag);

% Extract vectors from simulated data
%   Note: these are arrays with 25 columns
if strcmpi(buildMode, 'passive')
    [tVecs, vVecsSim, gVecsSim, iVecsSim, vVecsDend1, vVecsDend2] = ...
        extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                        GGABAB_COL_SIM, IEXT_COL_SIM, ...
                        DEND1_COL_SIM, DEND2_COL_SIM]);
elseif strcmpi(buildMode, 'active')
    [tVecs, vVecsSim, gVecsSim, iVecsSim, icaVecsSim, ...
            itmVecsSim, itminfVecsSim, ithVecsSim, ithinfVecsSim, ...
            ihVecsSim, ihmVecsSim, ikaVecsSim, iam1VecsSim, iah1VecsSim, ...
            iam2VecsSim, iah2VecsSim, ikkirVecsSim, ikirmVecsSim, ...
            inapnaVecsSim, inapmVecsSim, inaphVecsSim] = ...
        extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                        GGABAB_COL_SIM, IEXT_COL_SIM, ...
                        ICA_COL_SIM, ITM_COL_SIM, ITMINF_COL_SIM, ...
                        ITH_COL_SIM, ITHINF_COL_SIM, ...
                        IH_COL_SIM, IHM_COL_SIM, ...
                        IKA_COL_SIM, IAM1_COL_SIM, IAH1_COL_SIM, ...
                        IAM2_COL_SIM, IAH2_COL_SIM, ...
                        IKKIR_COL_SIM, IKIRM_COL_SIM, ...
                        INAPNA_COL_SIM, INAPM_COL_SIM, INAPH_COL_SIM]);
end

% If requested, bootstrap both recorded and simulated responses 
%   according to a grouping condition
if bootstrapCprFlag && ~isempty(realData) && strcmpi(simMode, 'passive')
    % Print to standard output
    if verbose
        fprintf('Bootstrap-averaging results for %s ...\n', expStr);
    end

    % Decide on the combination method
    if bootstrapCprFlag
        method = 'bootmean';
    else
        error('Code logic error!');
    end

    % Combine both recorded and simulated responses 
    realData = compute_combined_data(realData, method, 'Grouping', grouping, ...
                                    'ColNum', [VOLT_COL_REC, CURR_COL_REC]);
    simData = compute_combined_data(simData, method, 'Grouping', grouping, ...
                                    'ColNum', [VOLT_COL_SIM, IEXT_COL_SIM]);

    % Re-extract columns
    [vVecsRec, iVecsRec, gVecsRec] = ...
        extract_columns(realData, [VOLT_COL_REC, CURR_COL_REC, COND_COL_REC]);
    if strcmpi(buildMode, 'passive')
        [tVecs, vVecsSim, gVecsSim, iVecsSim, vVecsDend1, vVecsDend2] = ...
            extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                            GGABAB_COL_SIM, IEXT_COL_SIM, ...
                            DEND1_COL_SIM, DEND2_COL_SIM]);
    elseif strcmpi(buildMode, 'active')
        [tVecs, vVecsSim, gVecsSim, iVecsSim, icaVecsSim, ...
                itmVecsSim, itminfVecsSim, ithVecsSim, ithinfVecsSim, ...
                ihVecsSim, ihmVecsSim, ikaVecsSim, iam1VecsSim, iah1VecsSim, ...
                iam2VecsSim, iah2VecsSim, ikkirVecsSim, ikirmVecsSim, ...
                inapnaVecsSim, inapmVecsSim, inaphVecsSim] = ...
            extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                            GGABAB_COL_SIM, IEXT_COL_SIM, ...
                            ICA_COL_SIM, ITM_COL_SIM, ITMINF_COL_SIM, ...
                            ITH_COL_SIM, ITHINF_COL_SIM, ...
                            IH_COL_SIM, IHM_COL_SIM, ...
                            IKA_COL_SIM, IAM1_COL_SIM, IAH1_COL_SIM, ...
                            IAM2_COL_SIM, IAH2_COL_SIM, ...
                            IKKIR_COL_SIM, IKIRM_COL_SIM, ...
                            INAPNA_COL_SIM, INAPM_COL_SIM, INAPH_COL_SIM]);
    end

    % Re-compute number of sweeps
    nSweeps = numel(realData);

    % Re-compute baseline noise and sweep weights
    [~, ~, baseNoise, sweepWeights] = ...
        compute_default_sweep_info(tVecs, vVecsRec, ...
                                    'BaseWindow', baseWindow);
end

% Analyze the responses and compare
if generateDataFlag
    % Print to standard output
    if verbose
        fprintf('Analyzing responses for %s ...\n', expStr);
    end

    % Compute the sampling interval
    siMs = compute_sampling_interval(tVecs);

    % Decide on spreadsheet names
    featuresFile = fullfile(outFolder, [expStr, '_features.csv']);
    testResultsFile = fullfile(outFolder, [expStr, '_test_results.csv']);

    % Parse the simulated responses
    featuresSim = analyze_response(vVecsSim, iVecsSim, siMs, simMode, ...
                                    'simulated', ipscTime, ipscpWindow);

    % Parse the recorded responses
    if ~isempty(realData)
        featuresRec = analyze_response(vVecsRec, iVecsRec, siMs, simMode, ...
                                        'recorded', ipscTime, ipscpWindow);

        % Combine the features tables
        featuresTable = vertcat(featuresSim, featuresRec);
    else
        featuresTable = featuresSim;
    end

    % Save the results
    writetable(featuresTable, featuresFile);

    % Print to standard output
    if verbose
        fprintf('Comparing simulated versus recorded responses for %s ...\n', ...
                expStr);
    end

    % Do the appropriate comparison test
    if strcmpi(simMode, 'passive')
        featuresToCompare = {'steadyAmplitude', 'tauSlow', 'tauFast'};
    elseif strcmpi(simMode, 'active')
        featuresToCompare = {'peakProm', 'peakWidth', ...
                                'peakTime', 'spikesPerPeak'};
    end

    % Test the difference of features between dataType
    testResults = test_var_difference(featuresTable, featuresToCompare, ...
                                'dataType', 'SheetName', testResultsFile, ...
                                'Prefix', expStr, 'OutFolder', outFolder);
end

% If requested, combine both recorded and simulated responses 
%   according to a grouping condition
if averageCprFlag && ~isempty(realData) && strcmpi(simMode, 'passive')
    % Decide on the combination method
    if averageCprFlag
        method = 'mean';
    else
        error('Code logic error!');
    end

    % Combine both recorded and simulated responses 
    realData = compute_combined_data(realData, method, 'Grouping', grouping, ...
                                    'ColNum', [VOLT_COL_REC, CURR_COL_REC]);
    simData = compute_combined_data(simData, method, 'Grouping', grouping, ...
                                    'ColNum', [VOLT_COL_SIM, IEXT_COL_SIM]);

    % Re-extract columns
    [vVecsRec, iVecsRec, gVecsRec] = ...
        extract_columns(realData, [VOLT_COL_REC, CURR_COL_REC, COND_COL_REC]);
    if strcmpi(buildMode, 'passive')
        [tVecs, vVecsSim, gVecsSim, iVecsSim, vVecsDend1, vVecsDend2] = ...
            extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                            GGABAB_COL_SIM, IEXT_COL_SIM, ...
                            DEND1_COL_SIM, DEND2_COL_SIM]);
    elseif strcmpi(buildMode, 'active')
        [tVecs, vVecsSim, gVecsSim, iVecsSim, icaVecsSim, ...
                itmVecsSim, itminfVecsSim, ithVecsSim, ithinfVecsSim, ...
                ihVecsSim, ihmVecsSim, ikaVecsSim, iam1VecsSim, iah1VecsSim, ...
                iam2VecsSim, iah2VecsSim, ikkirVecsSim, ikirmVecsSim, ...
                inapnaVecsSim, inapmVecsSim, inaphVecsSim] = ...
            extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                            GGABAB_COL_SIM, IEXT_COL_SIM, ...
                            ICA_COL_SIM, ITM_COL_SIM, ITMINF_COL_SIM, ...
                            ITH_COL_SIM, ITHINF_COL_SIM, ...
                            IH_COL_SIM, IHM_COL_SIM, ...
                            IKA_COL_SIM, IAM1_COL_SIM, IAH1_COL_SIM, ...
                            IAM2_COL_SIM, IAH2_COL_SIM, ...
                            IKKIR_COL_SIM, IKIRM_COL_SIM, ...
                            INAPNA_COL_SIM, INAPM_COL_SIM, INAPH_COL_SIM]);
    end

    % Re-compute number of sweeps
    nSweeps = numel(realData);

    % Re-compute baseline noise and sweep weights
    [~, ~, baseNoise, sweepWeights] = ...
        compute_default_sweep_info(tVecs, vVecsRec, ...
                                    'BaseWindow', baseWindow);
end

% Compare with recorded data
if ~isempty(realData)
    % Display message
    if verbose
        fprintf('Computing errors for %s ...\n', expStr);
    end

    % Calculate voltage residuals (simulated - recorded)
    residuals = compute_residuals(vVecsSim, vVecsRec);

    % TODO: Fix this when normalize2InitErrFlag is true
    %       initSwpError needs to be read from the previous errorStruct somehow
    %       initLtsError needs to be read from the previous errorStruct somehow
    initSwpError = NaN;
    initLtsError = NaN;
    if normalize2InitErrFlag
        error('initSwpError needs to be fixed first!')
    end

    % Calculate errors (sweep errors, LTS errors, etc.)
    errorStruct = compute_single_neuron_errors(vVecsSim, vVecsRec, ...
                    'ErrorMode', errorMode, 'TimeVecs', tVecs, ...
                    'IvecsSim', iVecsSim, 'IvecsRec', iVecsRec, ...
                    'FitWindow', fitWindow, 'BaseWindow', baseWindow, ...
                    'BaseNoise', baseNoise, 'SweepWeights', sweepWeights, ...
                    'LtsFeatureWeights', ltsFeatureWeights, ...
                    'MissedLtsError', missedLtsError, ...
                    'FalseLtsError', falseLtsError, ...
                    'Lts2SweepErrorRatio', lts2SweepErrorRatio, ...
                    'Match2FeatureErrorRatio', match2FeatureErrorRatio, ...
                    'NormalizeError', normalize2InitErrFlag, ...
                    'InitSwpError', initSwpError, ...
                    'InitLtsError', initLtsError, ...
                    'IpscTime', ipscTime, 'IpscPeakWindow', ipscPeakWindow, ...
                    'FileBase', fileBases, ...
                    'OutFolder', outFolder, 'Prefix', expStr, ...
                    'SaveLtsInfoFlag', saveLtsInfoFlag, ...
                    'SaveLtsStatsFlag', saveLtsStatsFlag, ...
                    'PlotIpeakFlag', plotIpeakFlag, ...
                    'PlotLtsFlag', plotLtsFlag, ...
                    'PlotStatisticsFlag', plotStatisticsFlag);


    % Extract just the sweep errors
    swpErrors = errorStruct.swpErrors;
else
    residuals = [];
    errorStruct = struct;
    swpErrors = [];
end

%% Plot figures
% Prepare for plotting
if plotFlag
    % Print message
    if strcmpi(simMode, 'passive')
        fprintf('UPDATING current pulse response figures for %s ...\n', expStr);
    elseif strcmpi(simMode, 'active')
        fprintf('UPDATING GABAB IPSC response figures for %s ...\n', expStr);
    end

    % Find the indices of the x-axis limit endpoints
    endPointsForPlots = find_window_endpoints(xLimits, tVecs);

    % Prepare vectors for plotting
    if strcmpi(buildMode, 'passive')
        [tVecs, residuals, vVecsRec, iVecsRec, gVecsRec, ...
            vVecsSim, iVecsSim, gVecsSim, vVecsDend1, vVecsDend2] = ...
            argfun(@(x) prepare_for_plotting(x, endPointsForPlots), ...
                    tVecs, residuals, vVecsRec, iVecsRec, gVecsRec, ...
                    vVecsSim, iVecsSim, gVecsSim, vVecsDend1, vVecsDend2);
    elseif strcmpi(buildMode, 'active')
        [tVecs, residuals, vVecsRec, iVecsRec, gVecsRec, ...
            vVecsSim, iVecsSim, gVecsSim, ...
            icaVecsSim, itmVecsSim, itminfVecsSim, ...
            ithVecsSim, ithinfVecsSim, ihVecsSim, ihmVecsSim, ...
            ikaVecsSim, iam1VecsSim, iah1VecsSim, ...
            iam2VecsSim, iah2VecsSim, ikkirVecsSim, ikirmVecsSim, ...
            inapnaVecsSim, inapmVecsSim, inaphVecsSim] = ...
            argfun(@(x) prepare_for_plotting(x, endPointsForPlots), ...
                    tVecs, residuals, vVecsRec, iVecsRec, gVecsRec, ...
                    vVecsSim, iVecsSim, gVecsSim, ...
                    icaVecsSim, itmVecsSim, itminfVecsSim, ...
                    ithVecsSim, ithinfVecsSim, ihVecsSim, ihmVecsSim, ...
                    ikaVecsSim, iam1VecsSim, iah1VecsSim, ...
                    iam2VecsSim, iah2VecsSim, ikkirVecsSim, ikirmVecsSim, ...
                    inapnaVecsSim, inapmVecsSim, inaphVecsSim);
    end
end

% Plot a comparison of simulated and recorded conductance traces
if plotConductanceFlag
    % Print to standard output
    fprintf('Plotting figure of conductance trace comparison for %s ...\n', ...
            expStr);

    % Select data to plot
    dataForConductanceComparison = {gVecsRec; gVecsSim};

    % Construct matching time vectors
    tVecsForConductanceComparison = repmat({tVecs}, [2, 1]);

    % Construct matching y labels
    yLabelsConductanceComparison = {'Conductance (uS)'; 'Conductance (uS)'};

    % Decide on figure name
    figName = fullfile(outFolder, [expStr, '_conductance_comparison.png']);

    % Decide on figure titles for each subplot
    if strcmpi(simMode, 'passive')
        figSubTitles = {'Recorded conductance during current pulse', ...
                        'Simulated conductance during current pulse'};
    else
        figSubTitles = {'Recorded GABA_B IPSC conductances', ...
                        'Simulated GABA_B IPSC conductances'};
    end
    strcat(figSubTitles, 'for', ' ', expStrForTitle);

    % Plot the conductance traces
    figHandle = set_figure_properties('ClearFigure', true, ...
                    'Visible', visibleStatus, ...
                    'FigNumber', figNumberConductance, ...
                    'Name', 'Conductance traces');
    hFig.conductanceComparison = ...
        plot_traces(tVecsForConductanceComparison, ...
                    dataForConductanceComparison, ...
                    'Verbose', false, 'PlotMode', 'parallel', ...
                    'SubplotOrder', 'list', ...
                    'ColorMode', 'byTraceInPlot', ...
                    'LegendLocation', 'suppress', ...
                    'ColorMap', colorMapParallel, ...
                    'XLimits', xLimits, ...
                    'LinkAxesOption', 'x', 'XUnits', 'ms', ...
                    'YLabel', yLabelsConductanceComparison, ...
                    'FigSubTitles', figSubTitles, 'FigHandle', figHandle, ...
                    'FigName', figName, 'FigTypes', figTypes, ...
                    'LineWidth', lineWidthParallel);

    % TODO: Add this?
    % if nSweeps == 4
    %     legend('Control', 'GAT1 Block', 'GAT3 Block', 'Dual Block')
    % end
end

% Plot a comparison of simulated and recorded current traces
if plotCurrentFlag
    % Print to standard output
    fprintf('Plotting figure of current trace comparison for %s ...\n', ...
            expStr);

    % Select data to plot
    dataForCurrentComparison = {iVecsRec; iVecsSim};

    % Construct matching time vectors
    tVecsForCurrentComparison = repmat({tVecs}, [2, 1]);

    % Construct matching y labels
    yLabelsCurrentComparison = {'Current (nA)'; 'Current (nA)'};

    % Decide on figure name
    figName = fullfile(outFolder, [expStr, '_current_comparison.png']);

    % Decide on figure titles for each subplot
    if strcmpi(simMode, 'passive')
        figSubTitles = {'Recorded current pulses', ...
                        'Simulated current pulses'};
    else
        figSubTitles = {'Recorded GABA_B IPSC currents', ...
                        'Simulated GABA_B IPSC currents'};
    end
    strcat(figSubTitles, 'for', ' ', expStrForTitle);

    % Plot the current traces
    figHandle = set_figure_properties('ClearFigure', true, ...
                    'Visible', visibleStatus, ...
                    'FigNumber', figNumberCurrent, ...
                    'Name', 'Current traces');
    hFig.currentComparison = ...
        plot_traces(tVecsForCurrentComparison, dataForCurrentComparison, ...
                    'Verbose', false, 'PlotMode', 'parallel', ...
                    'SubplotOrder', 'list', ...
                    'ColorMode', 'byTraceInPlot', ...
                    'LegendLocation', 'suppress', ...
                    'ColorMap', colorMapParallel, ...
                    'XLimits', xLimits, ...
                    'LinkAxesOption', 'x', 'XUnits', 'ms', ...
                    'YLabel', yLabelsCurrentComparison, ...
                    'FigSubTitles', figSubTitles, 'FigHandle', figHandle, ...
                    'FigName', figName, 'FigTypes', figTypes, ...
                    'LineWidth', lineWidthParallel);

    % TODO: Add this?
    % if nSweeps == 4
    %     legend('Control', 'GAT1 Block', 'GAT3 Block', 'Dual Block')
    % end
end

% Plot individual simulated traces against recorded traces
if plotIndividualFlag
    % Print to standard output
    fprintf('Plotting figure of individual voltage traces for %s ...\n', ...
            expStr);

    % Decide on figure title and figure name
    figTitle = sprintf('All traces for Experiment %s', expStrForTitle);
    figName = fullfile(outFolder, [expStr, '_individual.png']);

    % Decide on the axes to be linked
    if strcmp(simMode, 'passive')
        linkAxesOption = 'x';
    else
        linkAxesOption = 'xy';
    end

    % Plot the individual traces
    figHandle = set_figure_properties('ClearFigure', true, ...
                    'Visible', visibleStatus, ...
                    'FigNumber', figNumberIndividual, ...
                    'Name', 'All traces');
    hFig.individual = ...
        plot_fitted_traces(tVecs, vVecsSim, 'MinimalLabels', true, ...
                'DataToCompare', vVecsRec, 'PlotMode', 'parallel', ...
                'SubplotOrder', 'bycolor', 'ColorMode', 'byRow', ...
                'ColorMap', colorMapIndividual, ...
                'XLimits', xLimits, 'LinkAxesOption', linkAxesOption, ...
                'FitWindow', fitWindow, 'BaseWindow', baseWindow, ...
                'BaseNoise', baseNoise, 'SweepErrors', swpErrors, ...
                'FigTitle', figTitle, 'FigName', figName, ...
                'FigHandle', figHandle, 'FigTypes', figTypes, ...
                'PlotSwpWeightsFlag', plotSwpWeightsFlag);
end
    
% Plot residuals
if plotResidualsFlag
    % Print to standard output
    fprintf('Plotting figure of residual traces for %s ...\n', ...
            expStr);

    % Decide on figure title and file name
    figTitle = sprintf('Residuals for Experiment %s', expStrForTitle);
    figName = fullfile(outFolder, [expStr, '_residuals.png']);

    % Plot the individual traces
    figHandle = set_figure_properties('ClearFigure', true, ...
                    'Visible', visibleStatus, ...
                    'FigNumber', figNumberResiduals, ...
                    'Name', 'Residual traces');
    hFig.residuals = ...
        plot_fitted_traces(tVecs, residuals, 'MinimalLabels', true, ...
                'PlotMode', 'residuals', ...
                'SubplotOrder', 'bycolor', 'ColorMode', 'byRow', ...
                'ColorMap', colorMapIndividual, ...
                'XLimits', xLimits, 'LinkAxesOption', 'xy', ...
                'FitWindow', fitWindow, 'BaseWindow', baseWindow, ...
                'BaseNoise', baseNoise, 'SweepErrors', swpErrors, ...
                'FigTitle', figTitle, 'FigName', figName, ...
                'FigHandle', figHandle, 'FigTypes', figTypes, ...
                'PlotSwpWeightsFlag', plotSwpWeightsFlag);
end    

% Plot different types of traces with different conditions overlapped
if plotOverlappedFlag
    % Print to standard output
    fprintf('Plotting figure of overlapped traces for %s ...\n', ...
            expStr);

    % Compute processed data
    if strcmpi(buildMode, 'active')
        itm2hVecsSim = (itmVecsSim .^ 2) .* ithVecsSim;
        itminf2hinfVecsSim = (itminfVecsSim .^ 2) .* ithinfVecsSim;
    end

    % Select data to plot
    if strcmpi(buildMode, 'passive')
        dataForOverlapped = {vVecsSim; vVecsDend1; vVecsDend2; iVecsSim};
    elseif strcmpi(buildMode, 'active')
        dataForOverlapped = {vVecsSim; gVecsSim; iVecsSim; ...
                icaVecsSim; itm2hVecsSim; itminf2hinfVecsSim; ...
                itmVecsSim; itminfVecsSim; ithVecsSim; ithinfVecsSim; ...
                ihVecsSim; ihmVecsSim; ...
                ikaVecsSim; iam1VecsSim; iah1VecsSim; ...
                iam2VecsSim; iah2VecsSim; ikkirVecsSim; ikirmVecsSim; ...
                inapnaVecsSim; inapmVecsSim; inaphVecsSim};
    end

    % Construct matching y labels
    if strcmpi(buildMode, 'passive')
        yLabelsOverlapped = {'V_{soma} (mV)'; 'V_{dend1} (mV)'; ...
                            'V_{dend2} (mV)'; 'I_{stim} (nA)'};
    elseif strcmpi(buildMode, 'active')
        yLabelsOverlapped = {'V_{soma} (mV)'; 'g_{GABA_B} (uS)'; ...
                'I_{stim} (nA)'; 'I_{Ca} (mA/cm^2)'; ...
                'm^2h_{T}'; 'm_{\infty}^2h_{\infty,T}'; ...
                'm_{T}'; 'm_{\infty,T}'; 'h_{T}'; 'h_{\infty,T}'; ...
                'I_{h} (mA/cm^2)'; 'm_{h}'; 'I_{A} (mA/cm^2)'; ...
                'm_{1,A}'; 'h_{1,A}'; 'm_{2,A}'; 'h_{2,A}'; ...
                'I_{Kir} (mA/cm^2)'; 'm_{\infty,Kir}'; ...
                'I_{NaP} (mA/cm^2)'; 'm_{\infty,NaP}'; 'h_{NaP}'};
    end

    % Add recorded voltage on the top if exists
    if ~isempty(vVecsRec)
        dataForOverlapped = [{vVecsRec}; dataForOverlapped];
        yLabelsOverlapped = [{'V_{rec} (mV)'}; yLabelsOverlapped];
    end
        
    % Construct matching time vectors
    tVecsForOverlapped = repmat({tVecs}, size(dataForOverlapped));

    % Decide on figure title and file name
    figTitle = sprintf('Overlapped traces for Experiment %s', ...
                        expStrForTitle);
    figName = fullfile(outFolder, [expStr, '_overlapped.png']);

    % Plot overlapped traces
    % TODO: Integrate into m3ha_plot_simulated_traces.m
    nSubPlots = numel(yLabelsOverlapped);
    figHandle = set_figure_properties('AlwaysNew', true, ...
                    'Visible', visibleStatus, ...
                    'FigNumber', figNumberOverlapped, ...
                    'FigExpansion', [1, nSubPlots/16], ...
                    'Name', 'All simulated traces');
    hFig.overlapped = ...
        plot_traces(tVecsForOverlapped, dataForOverlapped, ...
                    'Verbose', false, 'PlotMode', 'parallel', ...
                    'SubplotOrder', 'list', ...
                    'ColorMode', 'byTraceInPlot', ...
                    'LegendLocation', 'suppress', ...
                    'ColorMap', colorMapParallel, ...
                    'XLimits', xLimits, ...
                    'LinkAxesOption', 'x', 'XUnits', 'ms', ...
                    'YLabel', yLabelsOverlapped, ...
                    'FigTitle', figTitle, 'FigHandle', figHandle, ...
                    'FigName', figName, 'FigTypes', figTypes, ...
                    'LineWidth', lineWidthParallel);

    % handles = m3ha_plot_simulated_traces('Directory', outFolder, ...
    %                                     'ColorMap', colorMapParallel);
end

% Print an empty line
if verbose
    fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nSweeps = decide_on_nSweeps (realData, nSweepsUser)
%% Returns the number of sweeps to run and compare

if ~isempty(realData)
    % Simulate as many sweeps as recorded data
    nSweeps = numel(realData);

    % Check whether the user provided a different nSweeps
    if ~isempty(nSweepsUser) && nSweepsUser ~= nSweeps
        fprintf('realData provided, so nSweepsUser ignored!\n\n');
    end
elseif ~isempty(nSweepsUser)
    nSweeps = nSweepsUser;
else
    nSweeps = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rowConditions, nRows] = ...
                decide_on_row_conditions (rowConditions, nSweeps, ...
                                            maxRowsWithOneOnly)
%% Decide on rowConditions and nRows

% Place each sweep on its own row if rowConditions not provided
if isempty(rowConditions)
    % Decide on the number of rows
    if nSweeps <= maxRowsWithOneOnly
        nRows = nSweeps;
    else
        nRows = ceil(sqrt(nSweeps));
    end

    % Label the rows 1, 2, ..., nRows
    rowConditions = transpose(1:nRows);
else
    % Get the number of rows for plotting
    nRows = size(rowConditions, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xLimits = decide_on_xlimits (fitWindow, baseWindow, ...
                                        simMode, plotMarkFlag)
%% Decide on x-axis limits

% Put all window endpoints together
allEndpoints = [baseWindow, fitWindow];
allEndpoints = allEndpoints(:);

if plotMarkFlag && strcmpi(simMode, 'active')
    xLimits = [2800, 4500];
else
    xLimits = [min(allEndpoints), max(allEndpoints)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vecs = prepare_for_plotting(vecs, endPointsForPlots)
%% Prepare vectors for plotting

% Restrict vectors to xLimits to save time on plotting
vecs = extract_subvectors(vecs, 'Endpoints', endPointsForPlots);

% Combine vectors into matrices
vecs = force_matrix(vecs, 'AlignMethod', 'leftAdjustPad');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function featuresTable = analyze_response (vVecs, iVecs, siMs, simMode, ...
                                            dataType, ipscTime, ipscPeakWindow)

% Hard-coded parameters
meanVoltageWindow = 0.5;    % width in ms for calculating mean voltage 
                            %   for input resistance calculations
verbose = false;

% Parse the response (generate statistics)
if strcmpi(simMode, 'passive')
    % Parse the pulse response
    featuresTable = parse_pulse_response(vVecs, siMs, 'PulseVectors', iVecs, ...
                                'SameAsPulse', true, 'DetectPeak', false, ...
                                'FitResponse', true, ...
                                'MeanValueWindowMs', meanVoltageWindow);
elseif strcmpi(simMode, 'active')
    % Parse the IPSC current
    ipscTable = parse_ipsc(iVecs, siMs, 'StimStartMs', ipscTime, ...
                            'PeakWindowMs', ipscPeakWindow, ...
                            'Verbose', verbose);

    % Extract peak delay
    ipscDelay = ipscTable.peakDelayMs;

    % Parse the LTS response
    featuresTable = parse_lts(vVecs, siMs, 'StimStartMs', ipscTime, ...
                                'MinPeakDelayMs', ipscDelay, ...
                                'Verbose', verbose);
end

% Count the number of rows
nRows = height(featuresTable);

% Add a column for dataType ('simulated' or 'recorded')
featuresTable.dataType = repmat({dataType}, nRows, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
