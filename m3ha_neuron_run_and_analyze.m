function [errorStruct, hFig, simData] = ...
                m3ha_neuron_run_and_analyze (neuronParamsTable, varargin)
%% Runs "one iteration" of NEURON (once for each of the sweeps)
% Usage: [errorStruct, hFig, simData] = ...
%               m3ha_neuron_run_and_analyze (neuronParamsTable, varargin)
% Explanation:
%       TODO
% 
% Outputs: 
%       TODO
%       simData     - simulated data
%                   specified as a numeric array
%                       or a cell array of numeric arrays
% Arguments:
%       neuronParamsTable   
%                   - table(s) of single neuron parameters with 
%                       parameter names as 'RowNames' and with variables:
%                       'Value': value of the parameter
%                       'LowerBound': lower bound of the parameter
%                       'UpperBound': upper bound of the parameter
%                       'JitterPercentage': jitter percentage of the parameter
%                       'IsLog': whether the parameter is 
%                                   to be varied on a log scale
%                   must be a 2d table or a cell array of 2d tables
%       varargin    - 'Hfig': handles structure for figures
%                   must be a TODO
%                   default == TODO
%                   - 'SimMode': simulation mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'passive' - simulate a current pulse response
%                       'active'  - simulate an IPSC response
%                   default == 'passive'
%                   - 'NSweeps': number of sweeps
%                   must be a positive integer scalar
%                   default == numel(realData) or 1
%                   - 'Prefix': prefix to prepend to file names
%                   must be a character array
%                   default == ''
%                   - 'OutFolder': the directory where outputs will be placed
%                   must be a string scalar or a character vector
%                   default == pwd
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
%                   - 'Normalize2InitErrFlag': whether to normalize errors
%                                               to initial errors
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'SaveParamsFlag': whether to save simulation parameters
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveSimCmdsFlag': whether to save simulation commands
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveStdOutFlag': whether to save standard outputs
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveSimOutFlag': whether to save simulation outputs
%                                           when there are no errors
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotIndividualFlag': whether to plot individual traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotResidualsFlag': whether to plot residuals
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotOverlappedFlag': whether to plot overlapped traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotConductanceFlag': whether to plot conductance traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotCurrentFlag': whether to plot current traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotIpeakFlag': whether to current peak analyses
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotLtsFlag': whether to plot vtrace/LTS/burst analyses
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotStatisticsFlag': whether to plot LTS statistics
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotSwpWeightsFlag': whether to show a green 'ON' 
%                                       for sweeps in use
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
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
%                   - 'CprWindow': current pulse response window in ms
%                   must be a numeric vector with 2 elements
%                   default == [2000, 2360]
%                   - 'IpscrWindow': IPSC response window in ms
%                   must be a numeric vector with 2 elements
%                   default == [2000, 10000]
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
%                   - 'GababAmp': GABA-B IPSC amplitude (uS)
%                   must be a numeric vector
%                   default == set of all 12 input waveforms
%                   - 'GababTrise': GABA-B IPSC rising time constant (ms)
%                   must be a numeric vector
%                   default == set of all 12 input waveforms
%                   - 'GababTfallFast': GABA-B IPSC falling phase 
%                                           fast component time constant (ms)
%                   must be a numeric vector
%                   default == set of all 12 input waveforms
%                   - 'GababTfallSlow': GABA-B IPSC falling phase 
%                                           slow component time constant (ms)
%                   must be a numeric vector
%                   default == set of all 12 input waveforms
%                   - 'GababWeight': GABA-B IPSC falling phase 
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
%                   default == [2100, 2250]
%                   - 'FitWindowIpscr': time window to fit (ms)
%                                       for IPSC response
%                   must be a numeric vector with 2 elements
%                   default == [2980, 10000]
%                   - 'BaseWindowCpr': baseline window (ms)
%                                       for current pulse response
%                   must be a numeric vector with 2 elements
%                   default == [2000, 2100]
%                   - 'BaseWindowIpscr': baseline window (ms)
%                                       for IPSC response
%                   must be a numeric vector with 2 elements
%                   default == [2000, 2980]
%                   - 'BaseNoiseCpr': baseline noise (mV)
%                                       for current pulse response
%                   must be a numeric vector
%                   default == 1 mV
%                   - 'BaseNoiseIpscr': baseline noise (mV)
%                                       for IPSC response
%                   must be a numeric vector
%                   default == 1 mV
%                   - 'SweepWeightsCpr': sweep weights 
%                                       for current pulse response
%                   must be a numeric vector
%                   default == 1
%                   - 'SweepWeightsIpscr': sweep weights
%                                       for IPSC response
%                   must be a numeric vector
%                   default == 1
%
% Requires:
%       ~/m3ha/optimizer4gabab/singleneuron4compgabab.hoc
%       cd/argfun.m
%       cd/bar_w_CI.m
%       cd/compute_default_sweep_info.m
%       cd/compute_residuals.m
%       cd/compute_rms_error.m
%       cd/compute_sampling_interval.m
%       cd/create_colormap.m
%       cd/extract_columns.m
%       cd/extract_subvectors.m
%       cd/find_IPSC_peak.m
%       cd/find_LTS.m
%       cd/find_window_endpoints.m
%       cd/force_matrix.m
%       cd/load_neuron_outputs.m
%       cd/match_time_points.m
%       cd/compute_single_neuron_errors.m
%       cd/m3ha_neuron_create_simulation_params.m
%       cd/m3ha_neuron_create_TC_commands.m
%       cd/m3ha_plot_individual_traces.m
%       cd/parse_pulse_response.m
%       cd/run_neuron.m
%       cd/save_all_figtypes.m
%       cd/test_difference.m
%
%       cd/save_params.m TODO
%
% Used by:    
%       cd/m3ha_fminsearch3.m
%       ~/m3ha/optimizer4gabab/optimizer_4compgabab.m

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
% 2016-10-14 - Updated outputs for find_IPSC_peak & find_LTS
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
% 2017-05-19 - Moved update_sweeps_figures() here from optimizer_4compgabab.m
% 2017-05-19 - Separated update_sweeps_figures() into 
%               decide_on_xlimits(), plot_overlapped_traces(),
%               plot_conductance_traces(), plot_current_traces() 
%               && m3ha_plot_individual_traces()
% 2017-05-22 - Changed line width and indentation
% 2017-05-23 - Removed modeselected from outparams and replaced with updated outparams.runmode
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
% 2017-11-09 - Now uses bar_w_CI.m instead of errorbar
% 2017-12-21 - Added some more output ~s to find_LTS.m 
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

%% Hard-coded parameters
validSimModes = {'active', 'passive'};
hocFile = 'singleneuron4compgabab.hoc';
timeToStabilize = 2000;
maxRowsWithOneOnly = 8;

%% Column numbers for recorded data
%   Note: Must be consistent with ResaveSweeps.m
TIME_COL_REC = 1;
VOLT_COL_REC = 2;
CURR_COL_REC = 3;

%% Column numbers for simulated data
%   Note: Must be consistent with singleneuron4compgabab.hoc
TIME_COL_SIM = 1;
VOLT_COL_SIM = 2;
DEND1_COL_SIM = 3;
DEND2_COL_SIM = 4;
ITM_COL_SIM = 5;
ITH_COL_SIM = 6;
ICA_COL_SIM = 7;
IHM_COL_SIM = 8;
IDCLAMP_COL_SIM = 9;
GGABAB_COL_SIM = 10;
iCP_COL_SIM = 11;
IEXT_COL_SIM = 12;

%% Default values for optional arguments
hFigDefault = '';               % no prior hFig structure by default
simModeDefault = 'active'; %'passive';     % simulate a current pulse response by default
nSweepsDefault = [];            % set later
prefixDefault = '';             % prepend nothing to file names by default
outFolderDefault = pwd;         % use the present working directory for outputs
                                %   by default
debugFlagDefault = false;       % not in debug mode by default
customHoldCurrentFlagDefault = 0; % don't use custom hold current by default
onHpcFlagDefault = false;       % not on a high performance computing
                                %   server by default
generateDataFlagDefault = false;% not generating surrogate data by default
averageCprFlagDefault = false;  % don't average current pulse responses 
                                %   according to vHold by default
normalize2InitErrFlagDefault = false;
saveParamsFlagDefault = true;   % save simulation parameters by default
saveSimCmdsFlagDefault = true;  % save simulation commands by default
saveStdOutFlagDefault = true;   % save standard outputs by default
saveSimOutFlagDefault = true;   % save simulation outputs by default
plotIndividualFlagDefault = true;   % all zoomed traces plotted by default
plotResidualsFlagDefault = true;    % all residuals plotted by default
plotOverlappedFlagDefault = true;   % all overlapped traces plotted by default
plotConductanceFlagDefault = true;  % all conductance traces plotted by default
plotCurrentFlagDefault = true;      % all current traces plotted by default
plotIpeakFlagDefault = true;        % current peak analyses plotted by default
plotLtsFlagDefault = true;          % LTS/burst analyses plotted by default
plotStatisticsFlagDefault = true;   % LTS & burst statistics plotted by default
plotSwpWeightsFlagDefault = 'auto'; % set in m3ha_plot_inidividual_traces.m
plotMarkFlagDefault = true;         % the way Mark wants plots to look
showSweepsFlagDefault = true;       % whether to show sweep figures
jitterFlagDefault = false;      % no jitter by default
cprWindowDefault = [0, 360] + timeToStabilize;      % (ms)
ipscrWindowDefault = [0, 8000] + timeToStabilize;   % (ms)
outFilePathDefault = 'auto';    % set later
tstopDefault = [];              % set later
realDataIpscrDefault = [];      % no data to compare against by default
realDataCprDefault = [];        % no data to compare against by default
holdPotentialIpscrDefault = -70;% (mV)
holdPotentialCprDefault = -70;  % (mV)
currentPulseAmplitudeIpscrDefault = -0.050;  % (nA)
currentPulseAmplitudeCprDefault = -0.050;  % (nA)
gababAmpDefault = [];           % set later
gababTriseDefault = [];         % set later
gababTfallFastDefault = [];     % set later
gababTfallSlowDefault = [];     % set later
gababWeightDefault = [];        % set later
holdCurrentIpscrDefault = 0;    % (nA)
holdCurrentCprDefault = 0;      % (nA)
holdCurrentNoiseIpscrDefault = 0;% (nA)
holdCurrentNoiseCprDefault = 0; % (nA)
rowConditionsIpscrDefault = []; % set later
rowConditionsCprDefault = [];   % set later
fitWindowCprDefault = [100, 250] + timeToStabilize;     % (ms)
fitWindowIpscrDefault = [980, 8000] + timeToStabilize;  % (ms)
baseWindowCprDefault = [0, 100] + timeToStabilize;      % (ms)
baseWindowIpscrDefault = [0, 980] + timeToStabilize;    % (ms)
baseNoiseCprDefault = [];        % set later
baseNoiseIpscrDefault = [];      % set later
sweepWeightsCprDefault = [];     % set later
sweepWeightsIpscrDefault = [];   % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'neuronParamsTable', ...
    @(x) validateattributes(x, {'table', 'cell'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'HFig', hFigDefault);
addParameter(iP, 'SimMode', simModeDefault, ...
    @(x) any(validatestring(x, validSimModes)));
addParameter(iP, 'NSweeps', nSweepsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
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
addParameter(iP, 'Normalize2InitErrFlag', normalize2InitErrFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveParamsFlag', saveParamsFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveSimCmdsFlag', saveSimCmdsFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveStdOutFlag', saveStdOutFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveSimOutFlag', saveSimOutFlagDefault, ...
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
addParameter(iP, 'CprWindow', cprWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
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
addParameter(iP, 'CurrentPulseAmplitudeIpscr', currentPulseAmplitudeIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'CurrentPulseAmplitudeCpr', currentPulseAmplitudeCprDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'GababAmp', gababAmpDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'GababTrise', gababTriseDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'GababTfallFast', gababTfallFastDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'GababTfallSlow', gababTfallSlowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'GababWeight', gababWeightDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'HoldCurrentIpscr', holdCurrentIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'HoldCurrentCpr', holdCurrentCprDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'HoldCurrentNoiseIpscr', holdCurrentNoiseIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'HoldCurrentNoiseCpr', holdCurrentNoiseCprDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
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
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'BaseNoiseIpscr', baseNoiseIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'SweepWeightsCpr', sweepWeightsCprDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'SweepWeightsIpscr', sweepWeightsIpscrDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Read from the Input Parser
parse(iP, neuronParamsTable, varargin{:});
hFig = iP.Results.HFig;
simMode = validatestring(iP.Results.SimMode, validSimModes);
nSweepsUser = iP.Results.NSweeps;
prefix = iP.Results.Prefix;
outFolder = iP.Results.OutFolder;
debugFlag = iP.Results.DebugFlag;
customHoldCurrentFlag = iP.Results.CustomHoldCurrentFlag;
onHpcFlag = iP.Results.OnHpcFlag;
generateDataFlag = iP.Results.GenerateDataFlag;
averageCprFlag = iP.Results.AverageCprFlag;
normalize2InitErrFlag = iP.Results.Normalize2InitErrFlag;
saveParamsFlag = iP.Results.SaveParamsFlag;
saveSimCmdsFlag = iP.Results.SaveSimCmdsFlag;
saveStdOutFlag = iP.Results.SaveStdOutFlag;
saveSimOutFlag = iP.Results.SaveSimOutFlag;
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
cprWindow = iP.Results.CprWindow;
ipscrWindow = iP.Results.IpscrWindow;
outFilePath = iP.Results.OutFilePath;
tstop = iP.Results.Tstop;
realDataIpscr = iP.Results.RealDataIpscr;
realDataCpr = iP.Results.RealDataCpr;
holdPotentialIpscr = iP.Results.HoldPotentialIpscr;
holdPotentialCpr = iP.Results.HoldPotentialCpr;
currentPulseAmplitudeIpscr = iP.Results.CurrentPulseAmplitudeIpscr;
currentPulseAmplitudeCpr = iP.Results.CurrentPulseAmplitudeCpr;
gababAmp = iP.Results.GababAmp;
gababTrise = iP.Results.GababTrise;
gababTfallFast = iP.Results.GababTfallFast;
gababTfallSlow = iP.Results.GababTfallSlow;
gababWeight = iP.Results.GababWeight;
holdCurrentIpscr = iP.Results.HoldCurrentIpscr;
holdCurrentCpr = iP.Results.HoldCurrentCpr;
holdCurrentNoiseIpscr = iP.Results.HoldCurrentNoiseIpscr;
holdCurrentNoiseCpr = iP.Results.HoldCurrentNoiseCpr;
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

% TODO
% vHold
% cpStartWindowOrig
% cprWindow

%% Preparation
% Decide on simulation-mode-dependent variables
CURR_COL_SIM = IEXT_COL_SIM;
if strcmpi(simMode, 'passive')
    realData = realDataCpr;
    currentPulseAmplitude = currentPulseAmplitudeCpr;
    holdPotential = holdPotentialCpr;
    holdCurrent = holdCurrentCpr;
    holdCurrentNoise = holdCurrentNoiseCpr;
    rowConditions = rowConditionsCpr;
    fitWindow = fitWindowCpr;
    baseWindow = baseWindowCpr;
    baseNoise = baseNoiseCpr;
    sweepWeights = sweepWeightsCpr;
    errorMode = 'SweepOnly';
elseif strcmpi(simMode, 'active')
    realData = realDataIpscr;
    currentPulseAmplitude = currentPulseAmplitudeIpscr;
    holdPotential = holdPotentialIpscr;
    holdCurrent = holdCurrentIpscr;
    holdCurrentNoise = holdCurrentNoiseIpscr;
    rowConditions = rowConditionsIpscr;
    fitWindow = fitWindowIpscr;
    baseWindow = baseWindowIpscr;
    baseNoise = baseNoiseIpscr;
    sweepWeights = sweepWeightsIpscr;
    errorMode = 'Sweep&LTS';
end

% Decide whether to plot anything
if plotOverlappedFlag || plotConductanceFlag || ...
    plotCurrentFlag || plotIndividualFlag || plotResidualsFlag ...
    plotFlag = true;
else
    plotFlag = false;
end

% Decide on the number of sweeps to run and compare
nSweeps = decide_on_nSweeps(realData, nSweepsUser);

% Decide on x-axis limits for plotting
if plotFlag
    xLimits = decide_on_xlimits(fitWindow, baseWindow, simMode, plotMarkFlag);
end

% Decide on figure numbers
if showSweepsFlag
    figNumberIndividual = 104;
else
    figNumberIndividual = [];
end

% Decide on rowConditions and nRows
[rowConditions, nRows] = ...
    decide_on_row_conditions(rowConditions, nSweeps, maxRowsWithOneOnly);

% Decide on the colors for each row in the plots
colorMap = create_colormap(nRows);

% Create output file paths if not provided
if strcmpi(outFilePath, 'auto')
    outFilePath = create_simulation_output_filenames(nSweeps, ...
                            'OutFolder', outFolder, 'Prefix', prefix);
end

% Create a table of simulation parameters
simParamsTable = m3ha_neuron_create_simulation_params(neuronParamsTable, ...
                        'Prefix', prefix, 'OutFolder', outFolder, ...
                        'SaveParamsFlag', saveParamsFlag, ...
                        'JitterFlag', jitterFlag, ...
                        'CprWindow', cprWindow, 'IpscrWindow', ipscrWindow, ...
                        'NSims', nSweeps, 'SimMode', simMode, ...
                        'OutFilePath', outFilePath, 'Tstop', tstop, ...
                        'HoldPotential', holdPotential, ...
                        'CurrentPulseAmplitude', currentPulseAmplitude, ...
                        'GababAmp', gababAmp, 'GababTrise', gababTrise, ...
                        'GababTfallFast', gababTfallFast, ...
                        'GababTfallSlow', gababTfallSlow, ...
                        'GababWeight', gababWeight, ...
                        'CustomHoldCurrentFlag', customHoldCurrentFlag, ...
                        'HoldCurrent', holdCurrent, ...
                        'HoldCurrentNoise', holdCurrentNoise);

% Create simulation commands to be read by NEURON
simCommands = m3ha_neuron_create_TC_commands(simParamsTable, ...
                        'Prefix', prefix, 'OutFolder', outFolder, ...
                        'SaveSimCmdsFlag', saveSimCmdsFlag);

% Extract vectors from recorded data
%   Note: these will be empty if realData not provided
[tVecs, vVecsRec, iVecsRec] = ...
    extract_columns(realData, [TIME_COL_REC, VOLT_COL_REC, CURR_COL_REC]);

% Compute baseline noise and sweep weights if not provided
[~, ~, baseNoise, sweepWeights] = ...
    compute_default_sweep_info(tVecs, vVecsRec, ...
            'BaseWindow', baseWindow, 'BaseNoise', baseNoise, ...
            'SweepWeights', sweepWeights);

%% Run NEURON
% Run NEURON with the hocfile and attached simulation commands
run_neuron(hocFile, 'SimCommands', simCommands, ...
                    'Prefix', prefix, 'OutFolder', outFolder, ...
                    'DebugFlag', debugFlag, 'OnHpcFlag', onHpcFlag, ...
                    'SaveStdOutFlag', saveStdOutFlag);

%% TODO: Break the function here into two

%% Analyze results
% Create an experiment identifier
expStr = prefix;

% Create an experiment identifier for title
expStrForTitle = strrep(expStr, '_', '\_');

% Load .out files created by NEURON
simDataOrig = load_neuron_outputs('FileNames', outFilePath, ...
                                'RemoveAfterLoad', ~saveSimOutFlag);

% If recorded data provided (tVecs not empty at this point), 
%   interpolate simulated data to match the time points of recorded data
% Note: This is necessary because CVODE (variable time step method) 
%       is applied in NEURON
if ~isempty(tVecs)
    simData = cellfun(@(x, y) match_time_points(x, y), ...
                        simDataOrig, tVecs, ...
                        'UniformOutput', false);
else
    simData = simDataOrig;
end

% Extract vectors from simulated data
%   Note: these are arrays with 12 columns
if strcmpi(simMode, 'passive')
    [tVecs, vVecsSim, iVecsSim, vVecsDend1, vVecsDend2] = ...
        extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                        CURR_COL_SIM, DEND1_COL_SIM, DEND2_COL_SIM]);
elseif strcmpi(simMode, 'active')
    [tVecs, vVecsSim, iVecsSim, itmVecsSim, ithVecsSim, ...
            icaVecsSim, ihmVecsSim, gVecsSim] = ...
        extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                        CURR_COL_SIM, ITM_COL_SIM, ITH_COL_SIM, ...
                        ICA_COL_SIM, IHM_COL_SIM, GGABAB_COL_SIM]);
end

% Analyze the responses and compare
if generateDataFlag
    % Compute the sampling interval
    siMs = compute_sampling_interval(tVecs);

    % Decide on spreadsheet names
    featuresFile = fullfile(outFolder, [expStr, '_features.csv']);
    testResultsFile = fullfile(outFolder, [expStr, '_test_results.csv']);

    % Parse the simulated responses
    featuresSim = analyze_response(vVecsSim, iVecsSim, siMs, ...
                                    simMode, 'simulated');

    % Parse the recorded responses
    if ~isempty(realData)
        featuresRec = analyze_response(vVecsRec, iVecsRec, siMs, ...
                                        simMode, 'recorded');

        % Combine the features tables
        featuresTable = vertcat(featuresSim, featuresRec);
    else
        featuresTable = featuresSim;
    end

    % Save the results
    writetable(featuresTable, featuresFile);

    % Do the appropriate comparison test
    if strcmpi(simMode, 'passive')
        featuresToCompare = {'steadyAmplitude', 'tauSlow', 'tauFast'};
    elseif strcmpi(simMode, 'active')
        featuresToCompare = {'TODO', 'TODO'};
    end

    % Test the difference of features between dataType
    testResults = test_difference(featuresTable, featuresToCompare, 'dataType', ...
                                    'SheetName', testResultsFile, ...
                                    'OutFolder', outFolder);
end

% If requested, average both recorded and simulated responses 
%   according to the vHold condition recorded by Christine
% TODO: Use m3ha_average_by_vhold from m3ha_import_raw_traces.m
if averageCprFlag && ~isempty(realData) && strcmpi(simMode, 'passive')
    % Average both recorded and simulated responses 
    [realData, simData] = ...
        m3ha_average_by_vhold(realData, simData, vHold, ...
                            cpStartWindowOrig, cprWindow);

    % Re-extract columns
    [vVecsRec, iVecsRec] = ...
        extract_columns(realData, [VOLT_COL_REC, CURR_COL_REC]);
    if strcmpi(simMode, 'passive')
        [tVecs, vVecsSim, iVecsSim, vVecsDend1, vVecsDend2] = ...
            extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                            CURR_COL_SIM, DEND1_COL_SIM, DEND2_COL_SIM]);
    elseif strcmpi(simMode, 'active')
        [tVecs, vVecsSim, iVecsSim, itmVecsSim, ithVecsSim, ...
                icaVecsSim, ihmVecsSim, gVecsSim] = ...
            extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                            CURR_COL_SIM, ITM_COL_SIM, ITH_COL_SIM, ...
                            ICA_COL_SIM, IHM_COL_SIM, GGABAB_COL_SIM]);
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
    % Calculate voltage residuals (simulated - recorded)
    residuals = compute_residuals(vVecsSim, vVecsRec);

    % TODO: Fix this when normalize2InitErrFlag is true
    initSwpError = 1;
    if normalize2InitErrFlag
        error('initSwpError needs to be fixed first!')
    end

    % Calculate errors (sweep errors, LTS errors, etc.)
    errorStruct = compute_single_neuron_errors(vVecsSim, vVecsRec, ...
                    'ErrorMode', errorMode, 'TimeVecs', tVecs, ...
                    'IvecsSim', iVecsSim, 'IvecsReal', iVecsRec, ...
                    'FitWindow', fitWindow, 'SweepWeights', sweepWeights, ...
                    'BaseWindow', baseWindow, 'BaseNoise', baseNoise, ...
                    'NormalizeError', normalize2InitErrFlag, ...
                    'InitSwpError', initSwpError);

    % Extract just the sweep errors
    swpErrors = errorStruct.swpErrors;
else
    residuals = [];
    errorStruct = struct;
    swpErrors = [];
end

%% Plot figures
if plotFlag
    % Print message
    if strcmpi(simMode, 'passive')
        fprintf('UPDATING current pulse response figures for %s ...\n', prefix);
    elseif strcmpi(simMode, 'active')
        fprintf('UPDATING GABAB IPSC response figures for %s ...\n', prefix);
    end

    % Find the indices of the x-axis limit endpoints
    endPointsForPlots = find_window_endpoints(xLimits, tVecs);

    % Prepare vectors for plotting
    if strcmpi(simMode, 'passive')
        [tVecs, vVecsSim, vVecsRec, residuals, ...
            iVecsSim, vVecsDend1, vVecsDend2] = ...
            argfun(@(x) prepare_for_plotting(x, endPointsForPlots), ...
                    tVecs, vVecsSim, vVecsRec, residuals, ...
                    iVecsSim, vVecsDend1, vVecsDend2);
    elseif strcmpi(simMode, 'active')
        [tVecs, vVecsSim, vVecsRec, residuals, ...
            iVecsSim, itmVecsSim, ithVecsSim, ...
            icaVecsSim, ihmVecsSim, gVecsSim] = ...
            argfun(@(x) prepare_for_plotting(x, endPointsForPlots), ...
                    tVecs, vVecsSim, vVecsRec, residuals, ...
                    iVecsSim, itmVecsSim, ithVecsSim, ...
                    icaVecsSim, ihmVecsSim, gVecsSim);
    end

    % Plot individual simulated traces against recorded traces
    if plotIndividualFlag
        % Print to standard output
        fprintf('Plotting figure of individual voltage traces for %s ...\n', ...
                expStr);

        % Decide on figure title and file name
        figTitle = sprintf('All traces for Experiment %s', expStrForTitle);
        figName = fullfile(outFolder, [expStr, '_individual.png']);

        % Plot the individual traces
        hFig.individual = ...
            m3ha_plot_individual_traces(tVecs, vVecsSim, ...
                    'DataToCompare', vVecsRec, 'PlotMode', 'parallel', ...
                    'SubplotOrder', 'bycolor', 'ColorMode', 'byRow', ...
                    'ColorMap', colorMap, ...
                    'XLimits', xLimits, 'LinkAxesOption', 'x', ...
                    'FitWindow', fitWindow, 'BaseWindow', baseWindow, ...
                    'BaseNoise', baseNoise, 'SweepErrors', swpErrors, ...
                    'FigTitle', figTitle, 'FigName', figName, ...
                    'FigNumber', figNumberIndividual, ...
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
        hFig.residuals = ...
            m3ha_plot_individual_traces(tVecs, residuals, ...
                    'PlotMode', 'residuals', ...
                    'SubplotOrder', 'bycolor', 'ColorMode', 'byRow', ...
                    'ColorMap', colorMap, ...
                    'XLimits', xLimits, 'LinkAxesOption', 'xy', ...
                    'FitWindow', fitWindow, 'BaseWindow', baseWindow, ...
                    'BaseNoise', baseNoise, 'SweepErrors', swpErrors, ...
                    'FigTitle', figTitle, 'FigName', figName, ...
                    'FigNumber', figNumberIndividual, ...
                    'PlotSwpWeightsFlag', plotSwpWeightsFlag);
    end    

    % Plot different types of traces with different conditions overlapped
    if plotOverlappedFlag
        % Print to standard output
        fprintf('Plotting figure of overlapped traces for %s ...\n', ...
                expStr);

        % Select data to plot
        if strcmpi(simMode, 'passive')
            dataForOverlapped = {vVecsSim; vVecsDend1; vVecsDend2; iVecsSim};
        elseif strcmpi(simMode, 'active')
            dataForOverlapped = {vVecsSim; iVecsSim; itmVecsSim; ithVecsSim; ...
                                    ihmVecsSim; icaVecsSim; gVecsSim};
        end

        % Construct matching y labels
        if strcmpi(simMode, 'passive')
            yLabelsOverlapped = {'V_{soma} (mV)'; 'V_{dend1} (mV)'; ...
                                'V_{dend2} (mV)'; 'I_{stim} (nA)'};
        elseif strcmpi(simMode, 'active')
            yLabelsOverlapped = {'V_{soma} (mV)'; 'I_{stim} (nA)'; ...
                                'm_{I_T}'; 'h_{I_T}'; 'm_{I_h}'; ...
                                'I_{Ca} (mA/cm^2)'; 'g_{GABA_B} (uS)'};
        end

        % Add recorded voltage on the top if exists
        if ~isempty(vVecsRec)
            dataForOverlapped = [{vVecsRec}; dataForOverlapped];
            yLabelsOverlapped = [{'V_{rec} (mV)'}; yLabelsOverlapped];
        end
            
        % Construct matching time vectors
        tVecsForOverlapped = repmat({tVecs}, size(dataForOverlapped));

        % Expand the colormap if necessary
        if nSweeps > nRows
            nColumns = ceil(nSweeps / nRows);
            nSlots = nColumns * nRows;
            colorMap = reshape(repmat(reshape(colorMap, 1, []), ...
                                nColumns, 1), nSlots, 3);
        end

        % Decide on figure title and file name
        figTitle = sprintf('Overlapped traces for Experiment %s', ...
                            expStrForTitle);
        figName = fullfile(outFolder, [expStr, '_overlapped.png']);

        % Plot overlapped traces
        hFig.overlapped = ...
            plot_traces(tVecsForOverlapped, dataForOverlapped, ...
                        'Verbose', false, 'PlotMode', 'parallel', ...
                        'SubplotOrder', 'list', 'ColorMode', 'byTraceInPlot', ...
                        'LegendLocation', 'suppress', ...
                        'ColorMap', colorMap, 'XLimits', xLimits, ...
                        'LinkAxesOption', 'x', 'XUnits', 'ms', ...
                        'YLabel', yLabelsOverlapped, ...
                        'FigTitle', figTitle, 'FigName', figName);
    end

    %% TODO TODO
    % if plotConductanceFlag
    %     hFig.conductance = ...
    %         plot_conductance_traces(realData, simData, outparams, hFig, nSweeps, colorMap, ncg, npercg, xlimitsMax);
    % end
    % if plotCurrentFlag
    %     hFig.current = ...
    %         plot_current_traces(realData, simData, outparams, hFig, nSweeps, colorMap, ncg, npercg, xlimitsMax);
    % end

    % Print an empty line
    fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        nRows = floor(sqrt(nSweeps));
    end

    % Label the rows 1, 2, ..., nRows
    rowConditions = transpose(1:nRows);
else
    % Get the number of rows for plotting
    nRows = size(rowConditions, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xLimits = decide_on_xlimits (fitWindow, baseWindow, simMode, plotMarkFlag)
%% Decide on x-axis limits

% Put all window endpoints together
allEndpoints = [baseWindow, fitWindow];
allEndpoints = allEndpoints(:);

if plotMarkFlag && strcmpi(simMode, 'active')
    xLimits = [2800, 4500];
else
    xLimits = [min(allEndpoints), max(allEndpoints)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vecs = prepare_for_plotting(vecs, endPointsForPlots)
%% Prepare vectors for plotting

% Restrict vectors to xLimits to save time on plotting
vecs = extract_subvectors(vecs, 'Endpoints', endPointsForPlots);

% Combine vectors into matrices
vecs = force_matrix(vecs, 'AlignMethod', 'leftAdjustPad');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function featuresTable = analyze_response (vVecs, iVecs, siMs, simMode, dataType)

% Hard-coded parameters
meanVoltageWindow = 0.5;    % width in ms for calculating mean voltage 
                            %   for input resistance calculations

% Parse the response (generate statistics)
if strcmpi(simMode, 'passive')
    % Parse the pulse response
    featuresTable = parse_pulse_response(vVecs, siMs, 'PulseVectors', iVecs, ...
                                'SameAsPulse', true, 'DetectPeak', false, ...
                                'FitResponse', true, ...
                                'MeanValueWindowMs', meanVoltageWindow);
elseif strcmpi(simMode, 'passive')
    % Parse the IPSC current
    % TODO
    ipscTable = parse_ipsc(iVecs, siMs);

    % Parse the IPSC response
    % TODO
    featuresTable = parse_lts(vVecs, siMs, ...
                            'MeanValueWindowMs', meanVoltageWindow);
end

% Count the number of rows
nRows = height(featuresTable);

% Add a column for dataType ('simulated' or 'recorded')
featuresTable.dataType = repmat({dataType}, nRows, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [realData, simData] = m3ha_average_by_vhold(realData, simData, vHold)
%% Average both recorded and simulated responses according to vHold
%% TODO: don't hard-code column numbers
%% TODO: pull out code into functions

%% Average both recorded and simulated responses
% Get the number of data points for the current pulse responses
ndpReal = size(realData{1}, 1);
ndpSim = size(simData{1}, 1);

% Find unique vHold values
vUnique = sort(unique(vHold));
nVhold = length(vUnique);

% Group the traces by unique vHold values, then average the grouped traces
realDataGrouped = cell(nVhold, 1);
simDataGrouped = cell(nVhold, 1);
realDataAveraged = cell(nVhold, 1);
simDataAveraged = cell(nVhold, 1);
for iVhold = 1:nVhold
    % Get the current vHold value
    vnow = vUnique(iVhold);

    % Collect all cpr traces with this vHold value
    realDataGroupedThisVhold = realData(vHold == vnow);
    simDataGroupedThisVhold = simData(vHold == vnow);

    % Preallocate
    realDataAveragedThisVhold = zeros(ndpReal, 4);
    simDataAveragedThisVhold = zeros(ndpSim, 11);

    % Take the time, current and conductance traces from the first trace
    for iCol = [1, 3, 4]
        realDataAveragedThisVhold(:, iCol) = ...
            realDataGroupedThisVhold{1}(:, iCol);
    end
    for iCol = [1, 3:11]
        simDataAveragedThisVhold(:, iCol) = ...
            simDataGroupedThisVhold{1}(:, iCol);
    end

    % Average the recorded voltage traces
    temp1 = cellfun(@(x) x(:, 2), realDataGroupedThisVhold, ...
                    'UniformOutput', false);
    vRealGroupedThisVhold = cell2mat(temp1');
    vRealAveragedThis = nanmean(vRealGroupedThisVhold, 2);
    realDataAveragedThisVhold(:, 2) = vRealAveragedThis;

    % Average the simulated voltage traces
    temp2 = cellfun(@(x) x(:, 2), simDataGroupedThisVhold, ...
                    'UniformOutput', false);
    vSimGroupedThisVhold = cell2mat(temp2');
    vSimAveragedThis = nanmean(vSimGroupedThisVhold, 2);
    simDataAveragedThisVhold(:, 2) = vSimAveragedThis;

    % Save in arrays
    realDataGrouped{iVhold} = realDataGroupedThisVhold;
    realDataAveraged{iVhold} = realDataAveragedThisVhold;
    simDataGrouped{iVhold} = simDataGroupedThisVhold;
    simDataAveraged{iVhold} = simDataAveragedThisVhold;
end

% Replace with averaged data
realData = realDataAveraged;
simData = simDataAveraged;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [err, outparams] = errorcalc(realData, simData, outparams, nSweeps, colorMap, ncg, npercg)
%% Error calculator
%% TODO: Make functions out of this

%% Extract info from outparams
% Needed no matter what
outFolder = outparams.outFolder;
prefix = outparams.prefix;
cprflag = outparams.cprflag;
findLtsFlag = outparams.findLtsFlag;
ltsBurstStatsFlag = outparams.ltsBurstStatsFlag;
ltsErrorFlag = outparams.ltsErrorFlag;
plotIpeakFlag = outparams.plotIpeakFlag;
plotLtsFlag = outparams.plotLtsFlag;
plotStatisticsFlag = outparams.plotStatisticsFlag;
isLtsError = outparams.isLtsError;

if cprflag
    sweepWeights = outparams.sweepWeightsCpr;
else
    sweepWeights = outparams.sweepWeightsIpscr;
end

% Needed no matter and convert to columns
ioffset = outparams.ioffset;
peakprom = outparams.peakprom;
peakwidth = outparams.peakwidth;
ioffset = ioffset(:);
peakprom = peakprom(:);
peakwidth = peakwidth(:);

% Needed no matter what because of parfor
filenames = outparams.filenames;            
timeIPSC = outparams.ipscTime;
ipscpwin = outparams.ipscpWindow;
ltsWindow = outparams.ltsWindow;
if cprflag
    baseNoise = outparams.baseNoiseCpr;
else
    baseNoise = outparams.baseNoiseIpscr;
end
if ~cprflag && findLtsFlag && ltsBurstStatsFlag
    cellID = outparams.cellID;
end

if ltsErrorFlag
    ltsWeights = outparams.ltsWeights;
    lts2SweepErrorRatio = outparams.lts2SweepErrorRatio;
end

%% Calculate constants for efficiency
totalSweepWeights = sum(sweepWeights);
totalltsw = sum(ltsWeights);

%% Calculate maximum noise
if ~cprflag
    maxNoise = 4 * baseNoise;
end

%% Find regions of sweep trace for fitting (ms)
if cprflag
    fitWindow = outparams.fitWindowCpr;
    baseWindow = outparams.baseWindowCpr;
else
    fitWindow = outparams.fitWindowIpscr;
    baseWindow = [];
end

%% Set up vectors for parfor
a = zeros(nSweeps, 1);                  % column vectors

% For sweep error
rmse = a;                               % root-mean-squared error (mV) for each sweep

% For finding LTS
real_ipeakt = a;
sim_ipeakt = a;

% For LTS error & statistics
real_ltsv = a;
real_ltst = a;
real_ltsdvdtv = a;
sim_ltsv = a;
sim_ltst = a;
sim_ltsdvdtv = a;

% For completeness
real_ltsdvdtt = a;
sim_ltsdvdtt = a;

% For LTS statistics
real_pkprom = a;
real_pkwidth = a;
real_pkclass = a;
real_np2der = a;
real_spp = a;
real_btime = a;
real_spb = a;
sim_pkprom = a;
sim_pkwidth = a;
sim_pkclass = a;
sim_np2der = a;
sim_spp = a;
sim_btime = a;
sim_spb = a;

tFit = cell(nSweeps, 1);
vsim = cell(nSweeps, 1);
vreal = cell(nSweeps, 1);
isim = cell(nSweeps, 1);
ireal = cell(nSweeps, 1);

%% Crop and obtain required vectors, compute sum of squares error,
%   find LTS, extract features and compare statistics

% Get the fit region lower and upper bounds
fitWindowLowerBound = fitWindow(:, 1);
fitWindowUpperBound = fitWindow(:, 2);

% Extract the fit regions and compute root-mean-squared error
parfor iSwp = 1:nSweeps
    % Get the time vector for this sweep
    timeVec = realData{iSwp}(:, 1);

    % Find the indices of the time vector for fitting
    indFitWin = find(timeVec >= fitWindowLowerBound(iSwp) & ...
                     timeVec <= fitWindowUpperBound(iSwp));

    % Extract the vectors for fitting
    tFit{iSwp} = timeVec(indFitWin);
    vreal{iSwp} = realData{iSwp}(indFitWin, 2);
    vsim{iSwp} = simData{iSwp}(indFitWin, 2);
    ireal{iSwp} = realData{iSwp}(indFitWin, 3);
    isim{iSwp} = simData{iSwp}(indFitWin, 9);

    % Compute root-mean-squared error (mV) over the fit window
    rmse(iSwp) = compute_rms_error(vreal{iSwp}, vsim{iSwp});
end

if ~cprflag && findLtsFlag
    % For active fitting, also find low-threshold spikes in both real and simulated data
    % Note: The functions find_IPSC_peak.m & find_LTS.m are under /home/Matlab/Adams_Functions
    parfor iSwp = 1:nSweeps
        % Extract file base from matfile name
        [~, fileBase, ~] = fileparts(filenames{iSwp});
        fileBaseReal = strcat(fileBase, '_real');
        fileBaseSim = strcat(fileBase, '_sim');
        
        % Find IPSC peaks for real data
        [real_ipeakt(iSwp), ~, ~, ~] = ...
            find_IPSC_peak(tFit{iSwp}, ireal{iSwp}, timeIPSC, ipscpwin, ...
                plotIpeakFlag, outFolder, fileBaseReal);

        % Find LTSs (if any) for real data
        [~, ~, ~, real_np2der(iSwp), real_pkprom(iSwp), real_pkwidth(iSwp), real_pkclass(iSwp), ...
            real_spp(iSwp), real_ltst(iSwp), real_ltsv(iSwp), ...
            real_ltsdvdtt(iSwp), real_ltsdvdtv(iSwp), real_btime(iSwp), real_spb(iSwp), ~, ~, ~, ~, ~, ~, ~, ~] = ...
            find_LTS(tFit{iSwp}, vreal{iSwp}, timeIPSC, real_ipeakt(iSwp), ...
                maxNoise(iSwp), ltsWindow, plotLtsFlag, outFolder, fileBaseReal);

        % Find IPSC peaks for simulated data
        [sim_ipeakt(iSwp), ~, ~, ~] = ...
            find_IPSC_peak(tFit{iSwp}, isim{iSwp}, timeIPSC, ipscpwin, ...
                plotIpeakFlag, outFolder, fileBaseSim);

        % Find LTSs (if any) for simulated data
        [~, ~, ~, sim_np2der(iSwp), sim_pkprom(iSwp), sim_pkwidth(iSwp), sim_pkclass(iSwp), ...
            sim_spp(iSwp), sim_ltst(iSwp), sim_ltsv(iSwp), ...
            sim_ltsdvdtt(iSwp), sim_ltsdvdtv(iSwp), sim_btime(iSwp), sim_spb(iSwp), ~, ~, ~, ~, ~, ~, ~, ~] = ...
            find_LTS(tFit{iSwp}, vsim{iSwp}, timeIPSC, sim_ipeakt(iSwp), ...
                maxNoise(iSwp), ltsWindow, plotLtsFlag, outFolder, fileBaseSim);

    end
end

%% Compare LTS and burst statistics if GABAB IPSC response is simulated
if ~cprflag && findLtsFlag && ltsBurstStatsFlag
    % Save statistics first
    if outparams.saveLtsInfoFlag
        save(fullfile(outFolder, [prefix, '_LTS_info.mat']), ...
            'real_ipeakt', 'real_ltsv', 'real_ltst', ...
            'real_ltsdvdtv', 'real_ltsdvdtt', ...
            'real_pkprom', 'real_pkwidth', 'real_pkclass', ...
            'real_np2der', 'real_spp', 'real_btime', 'real_spb', ...
            'sim_ipeakt', 'sim_ltsv', 'sim_ltst', ...
            'sim_ltsdvdtv', 'sim_ltsdvdtt', ...
            'sim_pkprom', 'sim_pkwidth', 'sim_pkclass', ...
            'sim_np2der', 'sim_spp', 'sim_btime', 'sim_spb', ...
            '-v7.3');
    end
    compute_and_compare_statistics(nSweeps, colorMap, ncg, npercg, ...
        cellID, outparams, plotStatisticsFlag, ...
        real_ipeakt, real_ltsv, real_ltst, real_ltsdvdtv, real_ltsdvdtt, ...
        real_pkprom, real_pkwidth, real_pkclass, real_np2der, real_spp, real_btime, real_spb, ...
        sim_ipeakt, sim_ltsv, sim_ltst, sim_ltsdvdtv, sim_ltsdvdtt, ...
        sim_pkprom, sim_pkwidth, sim_pkclass, sim_np2der, sim_spp, sim_btime, sim_spb);
end

%% Compute average sweep errors:
% Get the root-mean-square error
err.swperr = rmse;

% Weight by baseline noise level, then make dimensionless by 
%       dividing by the first error if normalize2InitErrFlag is true
if cprflag
    fieldnameIniterr = 'cpr_avgswperr0';
else
    fieldnameIniterr = 'avgswperr0';
end
[err.avgSwpError, outparams] = ...
    compute_avgerror (fieldnameIniterr, err.swperr, outparams, ...
                        sweepWeights, totalSweepWeights);

%% Compute LTS errors and total errors
if ~cprflag && findLtsFlag && ltsErrorFlag

    % Store LTS features in error structure
    err.real_ltsv = real_ltsv;
    err.sim_ltsv = sim_ltsv;
    err.real_ltst = real_ltst;
    err.sim_ltst = sim_ltst;
    err.real_ltsdvdtv = real_ltsdvdtv;
    err.sim_ltsdvdtv = sim_ltsdvdtv;

    % Compute dimensionless LTS errors for each sweep
    [err.ltsAmpErrors, outparams.ltsSwpw] = ...
        compute_ltsdata_error (real_ltsv, sim_ltsv, baseNoise, isLtsError, sweepWeights);
    [err.ltsDelayErrors, ~] = ...
        compute_ltsdata_error (real_ltst, sim_ltst, peakwidth, isLtsError, sweepWeights);
    slopeUncertainty = (2*baseNoise ./ peakprom + 2*ioffset ./ peakwidth) ...
                            .* real_ltsdvdtv;
    [err.ltsSlopeErrors, ~] = ...
        compute_ltsdata_error (real_ltsdvdtv, sim_ltsdvdtv, ...
                                slopeUncertainty, isLtsError, sweepWeights);

    % Compute weighted-root-mean-squared-averaged LTS errors (dimensionless)
    ltsSwpw = outparams.ltsSwpw;
    ltsTotalSweepWeights = sum(outparams.ltsSwpw);
    [err.avgLtsAmpError, outparams] = ...
        compute_avgerror ('avgltsverr0', err.ltsAmpErrors, outparams, ltsSwpw, ltsTotalSweepWeights);
    [err.avgLtsDelayError, outparams] = ...
        compute_avgerror ('avgltsterr0', err.ltsDelayErrors, outparams, ltsSwpw, ltsTotalSweepWeights);
    [err.avgLtsSlopeError, outparams] = ...
        compute_avgerror ('avgltsdvdtverr0', err.ltsSlopeErrors, outparams, ltsSwpw, ltsTotalSweepWeights);
    
    % Average LTS error (dimensionless) is the weighted average of sweep error and lts error, 
    %        weighted by outparams.ltsWeights (set by user)
    err.avgLtsError = (ltsWeights' * [err.avgLtsAmpError; err.avgLtsDelayError; err.avgLtsSlopeError]) / totalltsw;

    % Total error (dimensionless) is the weighted average of sweep error and lts error, 
    %        weighted by outparams.lts2SweepErrorRatio (set by user)
    err.toterr = (err.avgSwpError + lts2SweepErrorRatio * err.avgLtsError) / (1 + lts2SweepErrorRatio);

else
    % Total error (dimensionless) is just the weighted-root-mean-squared average of sweep error
    err.toterr = err.avgSwpError;            

    % Set other errors to -999 for debug
    err.avgLtsAmpError = -999;
    err.avgLtsDelayError = -999;
    err.avgLtsSlopeError = -999;
    err.avgLtsError = -999;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ltsswperr, ltsSwpw] = compute_ltsdata_error (real_ltsdata, sim_ltsdata, uncertainty, isLtsError, sweepWeights)
% Compute LTS error for each sweep based on LTS existence
%% TODO: Make functions out of this

ltsswperr = zeros(size(real_ltsdata));
ltsSwpw = sweepWeights;                        % initialize lts sweep weights to user-defined sweep weights
for iSwp = 1:length(ltsswperr)
    if isnan(real_ltsdata(iSwp)) && isnan(sim_ltsdata(iSwp))          % LTS does not exist in both cases
        % Do not weigh this sweep in computing average error
        ltsSwpw(iSwp) = 0;
    elseif ~isnan(real_ltsdata(iSwp)) && ~isnan(sim_ltsdata(iSwp))    % LTS does exist in both cases
        % Compute dimensionless LTS error by normalizing to uncertainty in measurement
        ltsswperr(iSwp) = (real_ltsdata(iSwp) - sim_ltsdata(iSwp))/uncertainty(iSwp);
    else                                % LTS existence mismatch
        % Set to a predefined error (e.g., 10)
        ltsswperr(iSwp) = isLtsError;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [avgerror, outparams] = compute_avgerror (fieldnameIniterr, allErrors, outparams, allWeights, totalWeights)
% Compute averaged error over all sweeps and normalized to initial error if normalize2InitErrFlag is true
%% TODO: Make functions out of this

% Compute weighted-root-mean-squared-averaged error over all sweeps
avgerror = sqrt(sum(allWeights .* allErrors .^ 2) / totalWeights);

% Record initial error and normalize if normalize2InitErrFlag is true
if outparams.runnumtotal == 1
    % Store as starting error
    outparams.(fieldnameIniterr) = avgerror;

    % Starting normalized error is 1
    if outparams.normalize2InitErrFlag
        avgerror = 1;
    end
else
    % Divide by starting error
    if outparams.normalize2InitErrFlag
        avgerror = avgerror / outparams.(fieldnameIniterr);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ltst_all, ltst_mean_real, ltst_mean_sim, ltst_mean_all, ...
    ltst_std_real, ltst_std_sim, ltst_cv_real, ltst_cv_sim, ltst_cv_all, ltst_max_dev_all, ...
    ltsdvdtv_all, ltsdvdtv_mean_real, ltsdvdtv_mean_sim, ltsdvdtv_mean_all, ...
    ltsdvdtv_std_real, ltsdvdtv_std_sim, ltsdvdtv_cv_real, ltsdvdtv_cv_sim, ...
    ltsdvdtv_cv_all, ltsdvdtv_max_dev_all] ...
        = compute_and_compare_statistics(nSweeps, colorMap, ncg, npercg, ...
            cellID, outparams, plotStatisticsFlag, ...
            real_ipeakt, real_ltsv, real_ltst, real_ltsdvdtv, real_ltsdvdtt, ...
            real_pkprom, real_pkwidth, real_pkclass, real_np2der, real_spp, real_btime, real_spb, ...
            sim_ipeakt, sim_ltsv, sim_ltst, sim_ltsdvdtv, sim_ltsdvdtt, ...
            sim_pkprom, sim_pkwidth, sim_pkclass, sim_np2der, sim_spp, sim_btime, sim_spb)
%% Calculate statistics and plot those of simulated data against real data
%% TODO: Make functions out of this

%% Extract from outparams
outFolder = outparams.outFolder;
prefix = outparams.prefix;

%% Compute LTS statistics
ltsp = zeros(ncg, 2);            % LTS probability
ltsv_real = cell(ncg, 1);
ltst_real = cell(ncg, 1);
ltsdvdtv_real = cell(ncg, 1);
ltsdvdtt_real = cell(ncg, 1);
ltsv_sim = cell(ncg, 1);
ltst_sim = cell(ncg, 1);
ltsdvdtv_sim = cell(ncg, 1);
ltsdvdtt_sim = cell(ncg, 1);
both_has_lts_ct = zeros(ncg, 1);
ltst_all = cell(ncg, 1);
ltst_mean_real = zeros(ncg, 1);
ltst_mean_sim = zeros(ncg, 1);
ltst_mean_all = zeros(ncg, 1);
ltst_std_real = zeros(ncg, 1);
ltst_std_sim = zeros(ncg, 1);
ltst_std_all = zeros(ncg, 1);
ltst_cv_real = zeros(ncg, 1);
ltst_cv_sim = zeros(ncg, 1);
ltst_cv_all = zeros(ncg, 1);
ltst_max_dev_all = zeros(ncg, 1);
ltsdvdtv_all = cell(ncg, 1);
ltsdvdtv_mean_real = zeros(ncg, 1);
ltsdvdtv_mean_sim = zeros(ncg, 1);
ltsdvdtv_mean_all = zeros(ncg, 1);
ltsdvdtv_std_real = zeros(ncg, 1);
ltsdvdtv_std_sim = zeros(ncg, 1);
ltsdvdtv_std_all = zeros(ncg, 1);
ltsdvdtv_cv_real = zeros(ncg, 1);
ltsdvdtv_cv_sim = zeros(ncg, 1);
ltsdvdtv_cv_all = zeros(ncg, 1);
ltsdvdtv_max_dev_all = zeros(ncg, 1);
for cgn = 1:ncg            % color group number
    % LTS probability
    ct_real = 0;    % counts sweeps with LTS
    ct_sim = 0;    % counts sweeps with LTS
    for iSwp = 1:nSweeps
        if ceil(iSwp/npercg) == cgn && ~isnan(real_ltst(iSwp))
            ct_real = ct_real + 1;
        end
        if ceil(iSwp/npercg) == cgn && ~isnan(sim_ltst(iSwp))
            ct_sim = ct_sim + 1;
        end
    end
    ltsp(cgn, 1) = ct_real/npercg;
    ltsp(cgn, 2) = ct_sim/npercg;

    % Data for those sweeps that have LTSs both in real and simulated conditions
    ct = 0;         % counts sweeps with LTSs
    for iSwp = 1:nSweeps
        if ceil(iSwp/npercg) == cgn ...
            && ~isnan(real_ltst(iSwp)) && ~isnan(sim_ltst(iSwp))  
            ct = ct + 1;
            ltsv_real{cgn}(1, ct) = real_ltsv(iSwp);
            ltst_real{cgn}(1, ct) = real_ltst(iSwp);
            ltsdvdtv_real{cgn}(1, ct) = real_ltsdvdtv(iSwp);
            ltsdvdtt_real{cgn}(1, ct) = real_ltsdvdtt(iSwp);
            ltsv_sim{cgn}(1, ct) = sim_ltsv(iSwp);
            ltst_sim{cgn}(1, ct) = sim_ltst(iSwp);
            ltsdvdtv_sim{cgn}(1, ct) = sim_ltsdvdtv(iSwp);
            ltsdvdtt_sim{cgn}(1, ct) = sim_ltsdvdtt(iSwp);
        end
    end
    both_has_lts_ct(cgn) = ct;
    if both_has_lts_ct(cgn) > 0
        % LTS peak time data
        ltst_all{cgn} = [ltst_real{cgn} ltst_sim{cgn}];
        ltst_mean_real(cgn) = mean(ltst_real{cgn});
        ltst_mean_sim(cgn) = mean(ltst_sim{cgn});
        ltst_mean_all(cgn) = mean(ltst_all{cgn});
        ltst_std_real(cgn) = std(ltst_real{cgn});
        ltst_std_sim(cgn) = std(ltst_sim{cgn});
        ltst_std_all(cgn) = std([ltst_real{cgn} ltst_sim{cgn}]);
        ltst_cv_real(cgn) = ltst_std_real(cgn)/ltst_mean_real(cgn);
        ltst_cv_sim(cgn) = ltst_std_sim(cgn)/ltst_mean_sim(cgn);
        ltst_cv_all(cgn) = ltst_std_all(cgn)/ltst_mean_all(cgn);
        ltst_max_dev_all(cgn) = max(abs(ltst_all{cgn} - ltst_mean_all(cgn))/ltst_mean_all(cgn));

        % LTS max slope data
        ltsdvdtv_all{cgn} = [ltsdvdtv_real{cgn} ltsdvdtv_sim{cgn}];
        ltsdvdtv_mean_real(cgn) = mean(ltsdvdtv_real{cgn});
        ltsdvdtv_mean_sim(cgn) = mean(ltsdvdtv_sim{cgn});
        ltsdvdtv_mean_all(cgn) = mean(ltsdvdtv_all{cgn});
        ltsdvdtv_std_real(cgn) = std(ltsdvdtv_real{cgn});
        ltsdvdtv_std_sim(cgn) = std(ltsdvdtv_sim{cgn});
        ltsdvdtv_std_all(cgn) = std([ltsdvdtv_real{cgn} ltsdvdtv_sim{cgn}]);
        ltsdvdtv_cv_real(cgn) = ltsdvdtv_std_real(cgn)/ltsdvdtv_mean_real(cgn);
        ltsdvdtv_cv_sim(cgn) = ltsdvdtv_std_sim(cgn)/ltsdvdtv_mean_sim(cgn);
        ltsdvdtv_cv_all(cgn) = ltsdvdtv_std_all(cgn)/ltsdvdtv_mean_all(cgn);
        ltsdvdtv_max_dev_all(cgn) = max(abs(ltsdvdtv_all{cgn} - ltsdvdtv_mean_all(cgn))/ltsdvdtv_mean_all(cgn));
    end
end
if outparams.saveLtsStatsFlag
    save(fullfile(outFolder, [prefix, '_LTS_statistics.mat']), ...
        'ltst_all', 'ltst_mean_real', 'ltst_mean_sim', 'ltst_mean_all', ...
        'ltst_std_real', 'ltst_std_sim', 'ltst_cv_real', 'ltst_cv_sim', 'ltst_cv_all', 'ltst_max_dev_all', ...
        'ltsdvdtv_all', 'ltsdvdtv_mean_real', 'ltsdvdtv_mean_sim', 'ltsdvdtv_mean_all', ...
        'ltsdvdtv_std_real', 'ltsdvdtv_std_sim', 'ltsdvdtv_cv_real', 'ltsdvdtv_cv_sim', ...
        'ltsdvdtv_cv_all', 'ltsdvdtv_max_dev_all', '-v7.3');
end

%% Compute LTS/burst statistics

% Items to compute
statstitle = {'LTS onset time (ms)', 'LTS time jitter (ms)', 'LTS probability', 'Spikes per LTS', ...
        'Burst onset time (ms)', 'Burst time jitter (ms)', 'Burst probability', 'Spikes per burst'};
statsfilename = {'lts_onset_time', 'lts_time_jitter', 'lts_probability', 'spikes_per_lts', ...
                'burst_onset_time', 'burst_time_jitter', 'burst_probability', 'spikes_per_burst'};
pplabel2 = {'Con', 'GAT1', 'GAT3', 'Dual'};

% Initialize stats vectors
all_stats_real = cell(1, length(statstitle));
mean_stats_real = cell(1, length(statstitle));
std_stats_real = cell(1, length(statstitle));
ct_stats_real = cell(1, length(statstitle));
err_stats_real = cell(1, length(statstitle));
highbar_stats_real = cell(1, length(statstitle));
lowbar_stats_real = cell(1, length(statstitle));
all_stats_sim = cell(1, length(statstitle));
mean_stats_sim = cell(1, length(statstitle));
std_stats_sim = cell(1, length(statstitle));
ct_stats_sim = cell(1, length(statstitle));
err_stats_sim = cell(1, length(statstitle));
highbar_stats_sim = cell(1, length(statstitle));
lowbar_stats_sim = cell(1, length(statstitle));
for bi = 1:length(statstitle)
    all_stats_real{bi} = cell(ncg, 1);
    mean_stats_real{bi} = zeros(ncg, 1);
    std_stats_real{bi} = zeros(ncg, 1);
    ct_stats_real{bi} = zeros(ncg, 1);
    err_stats_real{bi} = zeros(ncg, 1);
    highbar_stats_real{bi} = zeros(ncg, 1);
    lowbar_stats_real{bi} = zeros(ncg, 1);

    all_stats_sim{bi} = cell(ncg, 1);
    mean_stats_sim{bi} = zeros(ncg, 1);
    std_stats_sim{bi} = zeros(ncg, 1);
    ct_stats_sim{bi} = zeros(ncg, 1);
    err_stats_sim{bi} = zeros(ncg, 1);
    highbar_stats_sim{bi} = zeros(ncg, 1);
    lowbar_stats_sim{bi} = zeros(ncg, 1);
end

for cgn = 1:ncg            % color group number
    thisp_ind = (cgn - 1) * npercg + (1:npercg);
    [all_stats, mean_stats, std_stats, ct_stats, err_stats, highbar_stats, lowbar_stats] = ...
        ltsburst_statistics(thisp_ind, cellID, real_ltst, real_spp, real_btime, real_spb);
%    [all_stats, mean_stats, std_stats, ct_stats, err_stats, highbar_stats, lowbar_stats] = ...
%        ltsburst_statistics(thisp_ind, cellID, outparams.ltspeaktime, outparams.spikesperpeak, outparams.bursttime, outparams.spikesperburst);
    for bi = 1:length(statstitle)
        all_stats_real{bi}{cgn} = all_stats{bi};
        mean_stats_real{bi}(cgn) = mean_stats(bi);
        std_stats_real{bi}(cgn) = std_stats(bi);
        ct_stats_real{bi}(cgn) = ct_stats(bi);
        err_stats_real{bi}(cgn) = err_stats(bi);
        highbar_stats_real{bi}(cgn) = highbar_stats(bi);
        lowbar_stats_real{bi}(cgn) = lowbar_stats(bi);
    end
    [all_stats, mean_stats, std_stats, ct_stats, err_stats, highbar_stats, lowbar_stats] = ...
        ltsburst_statistics(thisp_ind, cellID, sim_ltst, sim_spp, sim_btime, sim_spb);
    for bi = 1:length(statstitle)
        all_stats_sim{bi}{cgn} = all_stats{bi};
        mean_stats_sim{bi}(cgn) = mean_stats(bi);
        std_stats_sim{bi}(cgn) = std_stats(bi);
        ct_stats_sim{bi}(cgn) = ct_stats(bi);
        err_stats_sim{bi}(cgn) = err_stats(bi);
        highbar_stats_sim{bi}(cgn) = highbar_stats(bi);
        lowbar_stats_sim{bi}(cgn) = lowbar_stats(bi);
    end
end

% Plot statistics if plotStatisticsFlag == 1
if plotStatisticsFlag

    % Plot bar graph comparing LTS probabilities
    fprintf('Plotting bar graph comparing LTS probabilities ...\n');
    if outparams.showStatisticsFlag
        hFig.ltstp_bar = figure(201);
    else
        hFig.ltstp_bar = figure('Visible', 'off');
    end
    set(hFig.ltstp_bar, 'Name', 'Low threshold spike probability');
    clf(hFig.ltstp_bar);
    bar(1:size(colorMap, 1), ltsp);
    legend('Real data', 'Simulated data');
    xlabel('Pharm condition #');
    ylabel('LTS probability');
    title('Low threshold spike probability');
    figName = fullfile(outFolder, [prefix, '_ltstp_bar.png']);
    save_all_figtypes(hFig.ltstp_bar, figName);

    % Plot scatter plot of LTS onset times, don't save yet (axes not fixed)
    fprintf('Plotting scatter plot of LTS onset times ...\n');
    if outparams.showStatisticsFlag
        hFig.ltstcorr = figure(202);
    else
        hFig.ltstcorr = figure('Visible', 'off');
    end
    set(hFig.ltstcorr, 'Name', 'LTS onset times (ms)');
    clf(hFig.ltstcorr);
    for cgn = 1:size(colorMap, 1)            % color group number
        if both_has_lts_ct(cgn) > 0
            if ncg == 4
                subplot(2, 2, cgn); hold on;
                title(['Pharm condition ', num2str(cgn)])
            elseif ncg == 12
                subplot(4, 3, cgn); hold on;
                title(['Pharm condition ', num2str(floor(cgn/3) + 1)])
            end
            xlabel('Real data')
            ylabel('Simulated data')
            plot(ltst_real{cgn}, ltst_sim{cgn}, 'LineStyle', 'none', ...
                'Marker', 'o', 'MarkerEdgeColor', colorMap(cgn, :), 'MarkerFaceColor', colorMap(cgn, :));
        end
    end
    title('Correlation of LTS onset times (ms)')

    % Plot scatter plot of LTS max slopes, don't save yet (axes not fixed)
    fprintf('Plotting scatter plot of LTS max slopes ...\n');
    if outparams.showStatisticsFlag
        hFig.ltsdvdtvcorr = figure(203);
    else
        hFig.ltsdvdtvcorr = figure('Visible', 'off');
    end
    set(hFig.ltsdvdtvcorr, 'Name', 'LTS max slopes (mV/ms)');
    clf(hFig.ltsdvdtvcorr);
    for cgn = 1:size(colorMap, 1)            % color group number
        if both_has_lts_ct(cgn) > 0
            if ncg == 4
                subplot(2, 2, cgn); hold on;
                title(['Pharm condition ', num2str(cgn)])
            elseif ncg == 12
                subplot(4, 3, cgn); hold on;
                title(['Pharm condition ', num2str(floor(cgn/3) + 1)])
            end
            xlabel('Real data')
            ylabel('Simulated data')
            plot(ltsdvdtv_real{cgn}, ltsdvdtv_sim{cgn}, 'LineStyle', 'none', ...
                'Marker', 'o', 'MarkerEdgeColor', colorMap(cgn, :), 'MarkerFaceColor', colorMap(cgn, :));
        end
    end
    title('Correlation of LTS max slopes (mV/ms)')
    

    % Fix axes to make the scales consistent among subplots
    if sum(both_has_lts_ct) ~= 0
        ltst_max_dev_all_max = max(ltst_max_dev_all);
        ltst_std_all_max = max(ltst_std_all);
        if ltst_max_dev_all_max > 3 * ltst_std_all_max
            width = ltst_std_all_max;
        else
            width = ltst_max_dev_all_max;
        end

        % Don't let width be 0
        if width == 0
            width = 0.1;
        end

        set(0, 'CurrentFigure', hFig.ltstcorr);
        for cgn = 1:ncg            % color group number
            if both_has_lts_ct(cgn) > 0
                if ncg == 4
                    subplot(2, 2, cgn);
                elseif ncg == 12
                    subplot(4, 3, cgn);
                end                
                xmin = ltst_mean_all(cgn) * (1 - 1.1 * width);
                xmax = ltst_mean_all(cgn) * (1 + 1.1 * width);
                ymin = xmin;
                ymax = xmax;
                axis([xmin, xmax, ymin, ymax]);
            end
        end
        figName = fullfile(outFolder, [prefix, '_ltstcorr.png']);
        save_all_figtypes(hFig.ltstcorr, figName);

        ltsdvdtv_max_dev_all_max = max(ltsdvdtv_max_dev_all);
        ltsdvdtv_std_all_max = max(ltsdvdtv_std_all);
        if ltsdvdtv_max_dev_all_max > 3 * ltsdvdtv_std_all_max
            width = ltsdvdtv_std_all_max;
        else
            width = ltsdvdtv_max_dev_all_max;
        end

        % Don't let width be 0
        if width == 0
            width = 0.1;
        end

        set(0, 'CurrentFigure' , hFig.ltsdvdtvcorr);
        for cgn = 1:ncg            % color group number
            if both_has_lts_ct(cgn) > 0
                if ncg == 4
                    subplot(2, 2, cgn);
                elseif ncg == 12
                    subplot(4, 3, cgn);
                end                
                xmin = ltsdvdtv_mean_all(cgn) * (1 - 1.1 * width);
                xmax = ltsdvdtv_mean_all(cgn) * (1 + 1.1 * width);
                ymin = xmin;
                ymax = xmax;
                axis([xmin, xmax, ymin, ymax]);
            end
        end
        figName = fullfile(outFolder, [prefix, '_ltsdvdtvcorr.png']);
        save_all_figtypes(hFig.ltsdvdtvcorr, figName);
    end

    %% Create 2D bar graphs for LTS/burst statistics
    ltsburst_stats = cell(1, length(statstitle));    % for parfor
    parfor bi = 1:length(statstitle)
        if outparams.showStatisticsFlag
            ltsburst_stats{bi} = figure(210 + bi);
        else
            ltsburst_stats{bi} = figure('Visible', 'off');
        end
        fprintf('2D bar graph for %s ...\n', statsfilename{bi});
        set(ltsburst_stats{bi}, 'Name', [statsfilename{bi}]);
        clf(ltsburst_stats{bi});

        % Plot means with 95% confidence intervals
        bar_w_CI(ltsburst_stats{bi}, [mean_stats_real{bi} mean_stats_sim{bi}], ...
                    [lowbar_stats_real{bi} lowbar_stats_sim{bi}], ...
                    [highbar_stats_real{bi} highbar_stats_sim{bi}], ...
                    'XValues', (1:ncg)');

        set(gca, 'XTickLabel', pplabel2);
        if bi == 1
%            ylim([0 2500]);
        elseif bi == 2
%            ylim([0 800]);
        elseif bi == 3
%            ylim([0 1]);
        elseif bi == 4
%            ylim([0 6]);
        end
        legend('Real data', 'Simulated data');
        xlabel('Pharm Condition');
%        ylabel(statstitle{bi});
        title(statstitle{bi});
        figName = fullfile(outFolder, [prefix, '_', statsfilename{bi}]);
%        save_all_figtypes(ltsburst_stats{bi}, figName, {'png', 'fig'});
        save_all_figtypes(ltsburst_stats{bi}, figName, {'png'});
    end

    fprintf('\n');
    hFig.ltsburst_stats = ltsburst_stats;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hFig = plot_overlapped_traces (realData, simData, outparams, hFig, nSweeps, colorMap, ncg, npercg, xlimitsMax)
%% Update figure of all traces overlapped
%% TODO: Incorporate into plot_traces.m

%% Decide on the number of subplots
if outparams.cprflag
    nsubplots = 3;
else
    nsubplots = 5;
end

% See vMat of singleneuron4compgabab.hoc for simData columns
fprintf('Plotting figure of all traces overlapping for %s ...\n', outparams.prefix);
if outparams.showSweepsFlag
    if outparams.cprflag
        hFig.overlapped_traces = figure(101);
    else
        hFig.overlapped_traces = figure(111);
    end
else
    hFig.overlapped_traces = figure('Visible', 'off');
end
set(hFig.overlapped_traces, 'Name', 'Overlapped traces');
clf(hFig.overlapped_traces);
for iSwp = 1:nSweeps  
    % Find color group number
    cgn = ceil(iSwp/npercg);

    % Plot voltage traces from experiment
    subplot(nsubplots, 1, 1); hold on;
    plot(realData{iSwp}(:, 1), realData{iSwp}(:, 2), 'Color', colorMap(cgn, :), 'LineWidth', 1);

    % Plot voltage traces from simulations
    subplot(nsubplots, 1, 2); hold on;
    plot(simData{iSwp}(:, 1), simData{iSwp}(:, 2), 'Color', colorMap(cgn, :), 'LineWidth', 1);    % v of soma
    % plot(simData{iSwp}(:,1), simData{iSwp}(:,3), 'Color', 'c')                  % v of proximal dend
    % plot(simData{iSwp}(:,1), simData{iSwp}(:,4), 'Color', [0.4 0.3 0.8]);       % v of distal dend

    % Plot current traces from simulations
    % 20160722 - fixed from 9th column to 10th column
    subplot(nsubplots, 1, 3); hold on;
    if outparams.cprflag
        % Plot current pulse traces from simulations
        plot(simData{iSwp}(:, 1), simData{iSwp}(:, 11), 'Color', colorMap(cgn, :), ...
            'LineStyle', '-', 'LineWidth', 1);                 % cpi
    else
        % Plot GABAB IPSC current traces from simulations
        plot(simData{iSwp}(:, 1), simData{iSwp}(:, 9), 'Color', colorMap(cgn, :), ...
            'LineStyle', '-', 'LineWidth', 1);                 % GABABi
    end

    if ~outparams.cprflag
        % Plot activation & inactivation gates from simulations
        subplot(5,1,4); hold on;
        plot(simData{iSwp}(:, 1), simData{iSwp}(:, 5), 'Color', colorMap(cgn, :), ...
            'LineStyle', '-', 'LineWidth', 1);                % m_itGHK
        plot(simData{iSwp}(:, 1), simData{iSwp}(:, 6), 'Color', colorMap(cgn, :), ...
            'LineStyle', '--', 'LineWidth', 1);                % h_itGHK
        if ~outparams.plotMarkFlag
            plot(simData{iSwp}(:, 1), simData{iSwp}(:, 8), 'Color', colorMap(cgn, :), ...
                'LineStyle', ':', 'LineWidth', 1);            % m_Ih
        end

        % Plot calcium current from simulations
        subplot(5,1,5); hold on;
        plot(simData{iSwp}(:, 1), simData{iSwp}(:, 7), 'Color', colorMap(cgn, :), 'LineWidth', 1);        % ica
    end
end
subplot(nsubplots, 1, 1);
xlim(xlimitsMax);
if ~outparams.cprflag
    ylim([-110 -30]);
end
ylabel('Voltage (mV)', 'FontSize', 10);
legend('v\_soma (recorded)    !', 'Location', 'eastoutside');
if ~outparams.plotMarkFlag
    if outparams.cprflag
        title('All real and sim current pulse responses');
    else
        title('All real and sim traces');
    end
end
subplot(nsubplots, 1, 2);
xlim(xlimitsMax);
if ~outparams.cprflag
    ylim([-110 -30]);
end
ylabel('Voltage (mV)', 'FontSize', 10);
legend('v\_soma (simulated)   !', 'Location', 'eastoutside');

subplot(nsubplots, 1, 3);
ylabel('Current (nA)', 'FontSize', 10);
xlim(xlimitsMax);
if outparams.cprflag
    legend('cpi (simulated)           !', 'Location', 'eastoutside');
else
    legend('GABABi (simulated)   !', 'Location', 'eastoutside');
end

if ~outparams.cprflag
    subplot(5, 1, 4)
    xlim(xlimitsMax);
    ylabel('Probability', 'FontSize', 10);
    if outparams.plotMarkFlag
        legend('m\_itGHK (simulated) !', ...
            'h\_itGHK (simulated)', 'Location', 'eastoutside');
    else
        legend('m\_itGHK (simulated) !', ...
            'h\_itGHK (simulated)', 'm\_ih (simulated)', 'Location', 'eastoutside');
    end

    subplot(5, 1, 5)
    xlim(xlimitsMax);
    xlabel('Time (ms)');
    ylabel('Current (mA/cm^2)', 'FontSize', 10);
    legend('ica (simulated)          !', 'Location', 'eastoutside');
end

% Save figure
figName = fullfile(outparams.outFolder, [outparams.prefix, '_overlapped_traces']);
% save_all_figtypes(hFig.overlapped_traces, figName, {'png', 'fig'});
save_all_figtypes(hFig.overlapped_traces, figName, {'png'});
% Rename figure handle so it won't be overwritten
if outparams.cprflag
    hFig.cpr_overlapped_traces = hFig.overlapped_traces;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function hFig = plot_conductance_traces (realData, simData, outparams, hFig, nSweeps, colorMap, ncg, npercg, xlimitsMax)
%% Update figure comparing current traces between real data and simulations
%% TODO: Incorporate into plot_traces.m

% See vMat of singleneuron4compgabab.hoc for simData columns
fprintf('Plotting figure comparing current traces between real data and simulations for %s ...\n', outparams.prefix);
if outparams.showSweepsFlag
    if outparams.cprflag
        hFig.GABABi_comparison = figure(102);
    else
        hFig.GABABi_comparison = figure(112);
    end
else
    hFig.GABABi_comparison = figure('Visible', 'off');
end
set(hFig.GABABi_comparison, 'Name', 'GABAB currents');

% Plot GABAB IPSC current traces from experiments
clf(hFig.GABABi_comparison);
subplot(2,1,1); hold on;
for iSwp = 1:nSweeps  
    cgn = ceil(iSwp/npercg);        % color group number
    plot(realData{iSwp}(:, 1), realData{iSwp}(:, 3), 'Color', colorMap(cgn, :), 'LineStyle', '-');
end
if outparams.cprflag
    title('Current pulses, recorded')        %% TODO: To fix: It's not showing the current pulse!
else
    title('GABAB IPSC currents, recorded')
end
xlim(xlimitsMax);
if nSweeps == 4                % Legend only works if there are exactly 4 sweeps
    legend('Control', 'GAT1 Block', 'GAT3 Block', 'Dual Block')
end
xlabel('Time (ms)')
ylabel('Current (nA)')

% Plot current traces from simulations
subplot(2,1,2); hold on;
for iSwp = 1:nSweeps  
    cgn = ceil(iSwp/npercg);        % color group number
    if outparams.cprflag
        % Plot current pulse traces from simulations
        plot(simData{iSwp}(:, 1), simData{iSwp}(:, 11), 'Color', colorMap(cgn, :), 'LineStyle', '-');
    else
        % Plot GABAB IPSC current traces from simulations
        plot(simData{iSwp}(:, 1), simData{iSwp}(:, 9), 'Color', colorMap(cgn, :), 'LineStyle', '-');
    end
end
if outparams.cprflag
    title('Current pulses, simulated')
else
    title('GABAB IPSC currents, simulated')
end
xlim(xlimitsMax);
if nSweeps == 4                % Legend only works if there are exactly 4 sweeps
    legend('Control', 'GAT1 Block', 'GAT3 Block', 'Dual Block')
end
xlabel('Time (ms)')
ylabel('Current (nA)')
if outparams.cprflag
    figName = fullfile(outparams.outFolder, [outparams.prefix, '_cpi_comparison.png']);
else
    figName = fullfile(outparams.outFolder, [outparams.prefix, '_GABABi_comparison.png']);
end
save_all_figtypes(hFig.GABABi_comparison, figName);
% Copy figure handle so it won't be overwritten
if outparams.cprflag
    hFig.cpr_GABABi_comparison = hFig.GABABi_comparison;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hFig = plot_current_traces (realData, simData, outparams, hFig, nSweeps, colorMap, ncg, npercg, xlimitsMax)
%% Update figure comparing current traces between real data and simulations
%% TODO: Incorporate into plot_traces.m

% See vMat of singleneuron4compgabab.hoc for simData columns
fprintf('Plotting figure comparing current traces between real data and simulations for %s ...\n', outparams.prefix);
if outparams.showSweepsFlag
    if outparams.cprflag
        hFig.GABABg_comparison = figure(103);
    else
        hFig.GABABg_comparison = figure(113);
    end
else
    hFig.GABABg_comparison = figure('Visible', 'off');
end
set(hFig.GABABg_comparison, 'Name', 'GABAB conductances');
clf(hFig.GABABg_comparison);
subplot(2,1,1); hold on;
for iSwp = 1:nSweeps
    cgn = ceil(iSwp/npercg);        % color group number
    plot(realData{iSwp}(:, 1), realData{iSwp}(:, 4), 'Color', colorMap(cgn, :), 'LineStyle', '-');    
end
if outparams.cprflag
    title('Conductance during current pulse, recorded')
else
    title('GABAB IPSC conductances, recorded')
end
xlim(xlimitsMax);
if nSweeps == 4                % Legend only works if there are exactly 4 sweeps
    legend('Control', 'GAT1 Block', 'GAT3 Block', 'Dual Block')
end
xlabel('Time (ms)')
ylabel('Conductance (uS)')
subplot(2,1,2); hold on;
for iSwp = 1:nSweeps  
    cgn = ceil(iSwp/npercg);        % color group number
    plot(simData{iSwp}(:, 1), simData{iSwp}(:, 10), 'Color', colorMap(cgn, :), 'LineStyle', '-'); % 20160722 GABABg added
end
if outparams.cprflag
    title('Conductance during current pulse, simulated')
else
    title('GABAB IPSC conductances, simulated')
end
xlim(xlimitsMax);
if nSweeps == 4                % Legend only works if there are exactly 4 sweeps
    legend('Control', 'GAT1 Block', 'GAT3 Block', 'Dual Block')
end
xlabel('Time (ms)')
ylabel('Conductance (uS)')
if outparams.cprflag
    figName = fullfile(outparams.outFolder, [outparams.prefix, '_cpg_comparison.png']);
else
    figName = fullfile(outparams.outFolder, [outparams.prefix, '_GABABg_comparison.png']);
end
save_all_figtypes(hFig.GABABg_comparison, figName);
% Copy figure handle so it won't be overwritten
if outparams.cprflag
    hFig.cpr_GABABg_comparison = hFig.GABABg_comparison;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hFig = plot_voltage_residuals (tVecs, residuals, outparams, err, hFig, nSweeps, colorMap, ncg, npercg, xLimits)
%% Update figure of all residuals zoomed over fit region
%% TODO: Incorporate into plot_traces.m

if outparams.cprflag
    baseNoise = outparams.baseNoiseCpr;
    fitWindow = outparams.fitWindowCpr;
else
    baseNoise = outparams.baseNoiseIpscr;
    fitWindow = outparams.fitWindowIpscr;
end
swperr = err.swperr;

fprintf('Plotting figure of all residuals over fit region for %s ...\n', outparams.prefix);
if outparams.showSweepsFlag
    if outparams.cprflag
        hFig.alltraces_residuals = figure(105);
    else
        hFig.alltraces_residuals = figure(115);
    end
else
    hFig.alltraces_residuals = figure('Visible', 'off');
end
set(hFig.alltraces_residuals, 'Name', 'All traces zoomed');
clf(hFig.alltraces_residuals);
for iSwp = 1:nSweeps
    axesAll(iSwp) = subplot(ncg, npercg, iSwp); hold on;

    % Get vectors
    timeVec = tVecs{iSwp};
    residual = residuals{iSwp};
    
    % Plot traces
    plot(timeVec, residual, 'r');
    hold on;
    line([timeVec(1), timeVec(end)], [0, 0], 'Color', 'k', 'LineStyle', '--');
    % Adjust xLimits to show fitted region only
    xlim([xLimits(iSwp, 1), xLimits(iSwp, 2)]);
    if outparams.cprflag
        % ylim([-80 -60]);
    else
        % ylim([-110 -30]);
    end

    % If the trace is in use, show a green 'ON'
    if outparams.plotSwpWeightsFlag
        if sweepWeights(iSwp) ~= 0
            text('Units', 'normalized', 'Position', [0.1, 0.9], 'String', '\color{green} \bf ON');
        else
            text('Units', 'normalized', 'Position', [0.1, 0.9], 'String', '\color{gray} \bf OFF');
        end
    end

    % Show sweep info and error only if nSweeps <= 8
    if nSweeps <= 8
        title(['Noise = ', num2str(baseNoise(iSwp), 3), '; ', ...
                'RMSE = ', num2str(swperr(iSwp), 3)]);
    end

    % Plot fitWindow only if nSweeps <= 8
    if nSweeps <= 8
        yLimits = get(gca, 'YLim');
        line([fitWindow(1), fitWindow(1)], yLimits, 'Color', 'g', 'LineStyle', '--');
        line([fitWindow(2), fitWindow(2)], yLimits, 'Color', 'g', 'LineStyle', '--');
    end

    % Remove axes tick labels if nSweeps > 20
    if nSweeps > 20
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
    end
end

% Link x and y axes for all subplots
linkaxes(axesAll, 'xy');

% For suptitle
set(0, 'CurrentFigure', hFig.alltraces_residuals);

% If nSweeps > 20, expand all subplots by 1.2
% subplotsqueeze.m is in /home/Matlab/Downloaded_Functions
if nSweeps > 20
    subplotsqueeze(hFig.alltraces_residuals, 1.2);
end
suptitle(sprintf('All residuals for Experiment %s', strrep(outparams.prefix, '_', '\_')));        
                    % This is from the bioinformatics toolbox
% suplabel(sprintf('All traces for Experiment %s', strrep(outparams.prefix, '_', '\_')), 't');    
                    % This is in /home/Matlab/Downloaded_Functions
figName = fullfile(outparams.outFolder, [outparams.prefix, '_alltraces_residuals.png']);
save_all_figtypes(hFig.alltraces_residuals, figName);
% Copy figure handle so it won't be overwritten
if outparams.cprflag
    hFig.cpr_alltraces_residuals = hFig.alltraces_residuals;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%{
%% OLD CODE

% usedswps = outparams.sortedswpnum;
% if outparams.pfflag == 1
%     ccivamp = [0 0 0];
%     holdcurrent = outparams.holdcurrent;
%     holdpotential = outparams.holdpotential;
% else
%     %ccivamp = [0 0 0]; %-0.2 + 0.025 * (usedswps-1);  % in nA
%     ccivamp = outparams.ccivamp;
%     holdcurrent = outparams.holdcurrent;
%     holdpotential = outparams.holdpotential;
% end

% ccivamp = 0;

% assume -60 -65 -70
%% These will not be used
holdcurrent(find(holdpotential==-60)) = outparams.holdcurrent_minus60/1000;
holdcurrent(find(holdpotential==-65)) = outparams.holdcurrent_minus65/1000;
holdcurrent(find(holdpotential==-70)) = outparams.holdcurrent_minus70/1000;

%fprintf(fid,['sim("',outFolder,'/swp%g.out", '],usedswps(iSwp));  % params to be saved in .out file    

%% Build params.hoc file to be read by NEURON
%fid = fopen('params.hoc', 'w');
%    fprintf(fid,['sim("',outFolder,'/swp%d.out", '], iSwp);  % params to be saved in .out file    
%        fprintf(fid, '%6.4e, ',outparams.neuronparams(p));
%setfield(nowparams,outparams.neuronparamnames{p},outparams.neuronparams(p));       
%save (fullfile(outFolder,sprintf('/swp%g.mat',iSwp)), 'params');
%fprintf(fid,'%g, %g, %g )\n',ccivamp(iSwp),holdcurrent(iSwp),holdpotential(iSwp));

if vdcflag == 1
fprintf(fid,'%g, %g, %g, ', 0, holdcurrent(iSwp), holdpotential(iSwp));
fprintf(fid,'%g, %g, %g, %g, %g, %g)\n', ...
gababsyn.amp(iSwp)/1000,gababsyn.Trise(iSwp),gababsyn.TfallFast(iSwp), ...
gababsyn.TfallSlow(iSwp),gababsyn.w(iSwp),gababsyn.vHold(iSwp));
else

%    end
%fprintf(fid,'build()\n')
%    fprintf(fid,'%g, %g, %g, %g, %g, %g, %g, %g)\n', ...
%        0, holdcurrent(iSwp), holdpotential(iSwp), ...
%        gababsyn.amp(iSwp)/1000, gababsyn.Trise(iSwp), ...    % convert conductances from nS to uS
%        gababsyn.TfallFast(iSwp), gababsyn.TfallSlow(iSwp), gababsyn.w(iSwp), ...
%        outparams.tstop, outparams.currpulse/1000);    % convert current pulse from pA to nA
% movefile('params.hoc', fullfile(outFolder, 'params.backup'));

    %simData{iSwp} = load(fullfile(outFolder,['swp',num2str(usedswps(iSwp)),'.out']));

%    out = 1;
%    out = 0;

% 20131219 -> for passive param fit, pfflag = 1

%stati = cell(1, nSweeps);            % stores exit stati

%ltsWindow = [700 900]; %[6110 6400];
%ltsWindow = [3050 6500]; % - LTS window based on Christine's analysis of real data: 3050-6500 ms
%ltsWindow = [3000 9960]; % - This is the LTS window but is now defined in find_lts.m
%lb = round(ltsWindow(1)/di);
%ub = round(ltsWindow(2)/di);
%tlts = lb:di:ub; 
%ltsinterp1 = find((tvec0 > ltsWindow(1)),1,'first')-1;
%ltsinterp2 = find((tvec0 < ltsWindow(2)),1,'last')+1;
%ltsperiodind = (ltsinterp1:ltsinterp2);

%% Parameters used for finding LTS
mfw3 = 30;        % width in ms for the median filter for spikes (voltage traces)
mafw2 = 30;        % width in ms for the moving average filter for finding narrowest voltage peaks
minprom = 1;        % minimum LTS prominence in mV
maxNoiseprom = 0.1;    % maximum noise prominence in mV
lts_thr = -0.0023;    % 2nd derivative in V^2/s^2 below which defines an LTS peak
sp_thr = -30;        % amplitude threshold in mV for detecting a spike (the highest LTS peak is -34.01 mV)

% di = 2; %ms %%%% MAKE SURE THIS ISN'T TOO SHORT; DISTORTING LTS SHAPE?
% di = 1;

% vreal = zeros(length(indFitWin),numel(realData));
% vsim = zeros(length(indFitWin),numel(realData));
% vreal_lts = zeros(length(tlts),numel(realData));
% vsim_lts = zeros(length(tlts),numel(realData));

% swbnds = outparams.swpedges*1000; %[2000 8000]; % ms
% ltsbnds = [700 800]; %[6000 6600]; %ms

    clear val ind
    
    [val ind] = unique(simData{iSwp}(:,1)); % unique time points 
    vsim = interp1q(val,simData{iSwp}(ind,2),tvec0);
    vreal = realData{iSwp}(:,2);
    %%%% interpolate only the points we will be using to compute error
    %%%% (i.e. within error computing time window)  - 20131219
    
    clear val ind
    
    [val ind] = unique(simData{iSwp}(:,1)); % unique time points 
    vsim = interp1q(val,simData{iSwp}(ind,2),tvec0);
    vreal = realData{iSwp}(:,2);
    %%%% interpolate only the points we will be using to compute error
    %%%% (i.e. within error computing time window)  - 20131219
    

    
    %     indFitWin = indFitWin(1:2:end); % cut down on sample points of real data
    %         %%% MAKE SURE THIS IS OK!!  --> ok, verified 2011-01-30
    %     valind1 = find(val > swpreg(iSwp,1),1,'first')-1;
    %     valind2 = find(val < swpreg(iSwp,2),1,'last')+1;
    %     valperiodind = valind1:valind2;
    

        %mse_sst(iSwp) = (vreal-mean(vreal))'*(vreal-mean(vreal));

    
    treallong = tvec0(1:2:end);
    vreal = realData{iSwp}(1:2:end,2);
    swb = round(swpreg(1)/sims):round(swpreg(2)/sims);
    sse_sweep(iSwp) = (vreal(swb)-vsim(swb))'*(vreal(swb)-vsim(swb));
    sst_sweep(iSwp) = ((vreal(swb)-mean(vreal(swb))))'*((vreal(swb)-mean(vreal(swb))));
    sse_window(iSwp) = (vreal(indFitWin)-vsim(indFitWin))'*(vreal(indFitWin)-vsim(indFitWin));
    sst_window(iSwp) = (vreal(indFitWin)-mean(vreal(indFitWin)))'*(vreal(indFitWin)-mean(vreal(indFitWin)));
    

    % - compute lts error, EVEN IF NOT IN USE
    % compute lts error, IF IN USE
    %    if outparams.ltsErrorFlag == 1
    %    end

%%% compute rsquared (coefficient of determination)
%err.rsquared_swp = 1-(sse_sweep./sst_sweep);
%err.rsquared_window = 1-(sse_window./sst_window);

%swperr = rmse;

        legend('Real data', 'Simulated data');
%err.rsquared_window = rms; %1-(mse./mse_sst);  %%% 20131219
%ltserrv = sum((real_ltsv - sim_ltsv).^2);

sse_lts = a;
sst_lts = a;
b = cell(1, numel(realData));
real_ltsbnds = b;
sim_ltsbnds = b;

    err.rsquared_lts = 1 - (sse_lts./sst_lts);    % is this actually used??            % TO EXAMINE

sims = tvec0(2) - tvec0(1);            % sampling interval in ms

%    tbefore = tFit(1) - sims;                % time right before time vector of interest starts

    % Find low-threshold spikes
    [real_np2der(iSwp), real_spp(iSwp), real_ltst(iSwp), real_ltsv(iSwp), ...
        real_btime(iSwp), real_spb(iSwp), real_ltsbnds{iSwp}, real_ltsdvdtt(iSwp), real_ltsdvdtv(iSwp)] = ...
        find_lts(tFit, vreal, ireal, ...
            timeIPSC - tbefore, ipscpwin - tbefore, ltsWindow - tbefore, ...
            mfw3, mafw2, minprom, maxNoiseprom, lts_thr, sp_thr);        % this function is now find_lts_old.m
    [sim_np2der(iSwp), sim_spp(iSwp), sim_ltst(iSwp), sim_ltsv(iSwp), ...
        sim_btime(iSwp), sim_spb(iSwp), sim_ltsbnds{iSwp}, sim_ltsdvdtt(iSwp), sim_ltsdvdtv(iSwp)] = ...
        find_lts(tFit, vsim, isim, ...
            timeIPSC - tbefore, ipscpwin - tbefore, ltsWindow - tbefore, ...
            mfw3, mafw2, minprom, maxNoiseprom, lts_thr, sp_thr);

    % Use LTS bounds to find SSE & SST
    if ~isempty(real_ltsbnds{iSwp}) && ~isempty(sim_ltsbnds{iSwp})
        % Compute sse & sst in the LTS region
        real_ltsind = round(real_ltsbnds{iSwp}(1)/sims):round(real_ltsbnds{iSwp}(2)/sims);
        sim_ltsind = round(sim_ltsbnds{iSwp}(1)/sims):round(sim_ltsbnds{iSwp}(2)/sims);
        ltsind = intersect(real_ltsind, sim_ltsind);
        %ltsind = round(ltsbnds{iSwp}(1)/sims):round(ltsbnds{iSwp}(2)/sims);
        if isempty(ltsind) && ~isempty(real_ltsind)
            ltsind = real_ltsind;
        elseif isempty(ltsind) && ~isempty(sim_ltsind)
            ltsind = sim_ltsind;
        else
            ltsind = 1:length(vsim);
        end
        sse_lts(iSwp) = (vreal(ltsind) - vsim(ltsind))' * (vreal(ltsind) - vsim(ltsind));
        sst_lts(iSwp) = (vreal(ltsind) - mean(vreal(ltsind)))' * (vreal(ltsind) - mean(vreal(ltsind)));
    end

    
        vreal_lts(:,iSwp) = interp1q(tvec0(ltsperiodind)/di,realData{iSwp}(ltsperiodind,2),tlts');
        vsim_lts(:,iSwp) = interp1q(val/di,simData{iSwp}(ind,2),tlts');

        baseline = outparams.holdpotential(iSwp);              % - baseline membrane potential        
        [real_ltsv(iSwp) ind] = max(vreal_lts(:,iSwp));
        if real_ltsv(iSwp) > baseline + 5 % - only record lts peak if greater than baseline
                                                 % This is consistent with Christine's analysis of real data
            real_ltst(iSwp) = tlts(ind) - ltsWindow(1);
            real_diffvreal_lts = diff(vreal_lts(:,iSwp))/di;
            [real_ltsdvdtv(iSwp) ind] = max(real_diffvreal_lts);
            real_ltsdvdtt(iSwp) = tlts(ind) - ltsWindow(1);   % - not used
        else
            real_ltsv(iSwp) = baseline;   % - if there is no peak, assume the max is back to baseline
            real_ltst(iSwp) = -1;
            real_ltsdvdtv(iSwp) = -1;
        end
        
        [sim_ltsv(iSwp) ind] = max(vsim_lts(:,iSwp));
        if sim_ltsv(iSwp) > baseline + 5  % - only record lts peak if greater than baseline + 5 mV
                                                 % This is consistent with Christine's analysis of real data
            sim_ltst(iSwp) = tlts(ind) - ltsWindow(1);
            sim_diffvsim_lts = diff(vsim_lts(:,iSwp))/di;
            [sim_ltsdvdtv(iSwp) ind] = max(sim_diffvsim_lts);
            sim_ltsdvdtt(iSwp) = tlts(ind) - ltsWindow(1);
        else
            sim_ltsv(iSwp) = baseline;    % - if there is no peak, assume the max is back to baseline
            sim_ltst(iSwp) = -1;
            sim_ltsdvdtv(iSwp) = -1;
        end
    

    % check
    %     figure; plot(tlts,vreal_lts(:,iSwp)); hold on;
    %     plot(tlts(1:(end-1)),real_diffvreal_lts,'g');
    %     plot(tlts,vsim_lts(:,iSwp),'r');
    %     plot(real_ltst(iSwp),real_ltsv(iSwp),'iSwp*');
    %     plot(sim_ltst(iSwp),sim_ltsv(iSwp),'r*');
    %     plot(real_ltsdvdtt(iSwp),real_ltsdvdtv(iSwp),'ko');
    %     plot(sim_ltsdvdtt(iSwp),sim_ltsdvdtv(iSwp),'ro');

err_stats2_real{bi} = err_stats_real{bi} / 2;    % the function errorbar plots 2 * err_stat % That's the length of the error bar
err_stats2_sim{bi} = err_stats_sim{bi} / 2;    % the function errorbar plots 2 * err_stat % That's the length of the error bar

rmse = a;
    rmse(iSwp) = sqrt(mse(iSwp));                    % root mean squared error

    for p = 1:outparams.num_buildparams
        if p < outparams.num_buildparams
    for p = outparams.num_buildparams+1:numel(outparams.neuronparams)
            simCommands{iSwp} = [simCommands{iSwp}, sprintf('%g, ', outparams.neuronparams(p))];
            simCommands{iSwp} = [simCommands{iSwp}, sprintf('%g)\n', outparams.neuronparams(p))];

    simCommands{iSwp} = sprintf('build("%s", ', simMode);
    buildparams = outparams.neuronparams(outparams.neuronparamisbuild);
    for iSwp = 1:length(buildparams)
        if iSwp < length(buildparams)
            % Note: strcat would remove whitespacesm but [ , ] wouldn't
            simCommands{iSwp} = [simCommands{iSwp}, sprintf('%g, ', buildparams(iSwp))];
        else
            simCommands{iSwp} = [simCommands{iSwp}, sprintf('%g)\n', buildparams(iSwp))];
        end
    end

colorMap = {'k', 'b', 'r', 'g'};        % color groups: corresponding to 4 pharm conditions

mse = a;
mse(iSwp) = sum((vreal - vsim).^2) / length(indFitWin);        % mean squared error
err.swperr = mse;                % Use the mean squared errors as the sweep errors

err.totswperr = sweepWeights' .* err.swperr;        % weighted total sweep error
err.unweightedtotswperr = sum(err.swperr);    % unweighted total sweep error            %%% TO EXAMINE

a = zeros(1, nSweeps);

    err.ltsAmpErrors = real_ltsv - sim_ltsv;                % in mV
    err.ltsDelayErrors = (real_ltst - sim_ltst)/1000;            % in seconds
    err.ltsSlopeErrors = real_ltsdvdtv - sim_ltsdvdtv;            % in V/s
    rmse(iSwp) = sqrt(sum(((vreal - vsim)./vreal).^2) / length(indFitWin));    % root-mean-squared error for each sweep

    % Compute dimensionless LTS errors for each sweep
    err.ltsAmpErrors = (real_ltsv - sim_ltsv)./real_ltsv;
    err.ltsDelayErrors = (real_ltst - sim_ltst)./real_ltst;
    err.ltsSlopeErrors = (real_ltsdvdtv - sim_ltsdvdtv)./real_ltsdvdtv;

    % Convert NaN errors to outparams.isLtsError
    err.ltsAmpErrors(isnan(err.ltsAmpErrors)) = outparams.isLtsError;
    err.ltsDelayErrors(isnan(err.ltsDelayErrors)) = outparams.isLtsError;
    err.ltsSlopeErrors(isnan(err.ltsSlopeErrors)) = outparams.isLtsError;

    % Compute weighted-root-mean-squared average LTS errors
    err.avgLtsAmpError = sqrt(sum((sweepWeights .* err.ltsAmpErrors.^2) ./ totalSweepWeights));
    err.avgLtsDelayError = sqrt(sum((sweepWeights .* err.ltsDelayErrors.^2) ./ totalSweepWeights));
    err.avgLtsSlopeError = sqrt(sum((sweepWeights .* err.ltsSlopeErrors.^2) ./ totalSweepWeights));

if outparams.plotsweepsflag
    hFig = update_sweeps_figures(realData, simData, outparams, err, hFig);
end

if strcmp(outparams.modeselected,'modebutton_jitter') == 1

    for p = 1:outparams.numparams
        if ~outparams.neuronparamisbuild(p)
            simCommands{iSwp} = [simCommands{iSwp}, sprintf('%6.4e, ', outparams.neuronparams(p))];
        end
    end
        % Pass build parameters to the build() function
        if outparams.neuronparamisbuild(p)
            % Note: strcat would remove whitespacesm but [ , ] wouldn't
            simCommands{iSwp} = [simCommands{iSwp}, sprintf(', %g', outparams.neuronparams(p))];
        end

err.swperr = rmse;                        % use root-mean-squared errors [mV]
        legend('Real data', 'Simulated data');
if cprflag
    % weighted-root-mean-squared-averaged cpr sweep error [mV] or [1]
    [err.avgSwpError, outparams] = compute_avgerror ('cpr_avgswperr0', err.swperr, outparams, sweepWeights, totalSweepWeights);
else
    % weighted-root-mean-squared-averaged sweep error [mV] or [1]
    [err.avgSwpError, outparams] = compute_avgerror ('avgswperr0', err.swperr, outparams, sweepWeights, totalSweepWeights);
end

%% Compute sweep errors; make dimensionless by dividing by the holding potential
err.swperr = abs(rmse/mean(holdpotential));     % use root-mean-squared error 
                                                %   divided by the average holding potential

        % Compute dimensionless LTS error
        ltsswperr(iSwp) = (real_ltsdata(iSwp) - sim_ltsdata(iSwp))/real_ltsdata(iSwp);

    if isnan(real_ltsdata(iSwp)) && isnan(sim_ltsdata(iSwp))          % LTS does not exist in both cases
        % Do nothing (zero error)

    [err.ltsDelayErrors, ~] = ...
        compute_ltsdata_error (real_ltst, sim_ltst, ioffset, isLtsError, sweepWeights);

%        figName = fullfile(outFolder, [prefix, '_', statsfilename{bi}, '.png']);
%        figName = fullfile(outFolder, [prefix, '_', statsfilename{bi}, '.fig']);
%        saveas(ltsburst_stats{bi}, figName);
saveas(hFig.overlapped_traces, figName);

        % Plot means
        bar((1:ncg)', [mean_stats_real{bi} mean_stats_sim{bi}]); hold on;
        % Plot 95% confidence intervals
        errorbar((1:ncg)' - 0.13, mean_stats_real{bi}, err_stats_real{bi}, 'k');
        errorbar((1:ncg)' + 0.13, mean_stats_sim{bi}, err_stats_sim{bi}, 'k');


        legend('Real data', 'Simulated data');

set(hFig.ltstcorr, 'Visible', 'off');
set(ltsburst_stats{bi}, 'Visible', 'off');
set(hFig.alltraces_zoom, 'Visible', 'off');
set(hFig.GABABg_comparison, 'Visible', 'off');
set(hFig.GABABi_comparison, 'Visible', 'off');
set(hFig.overlapped_traces, 'Visible', 'off');

figure(hFig.alltraces_zoom);        % for suptitle
        figure(hFig.ltstcorr);
        figure(hFig.ltsdvdtvcorr);

% err.swperr = rmse./maxNoise';
% slopeUncertainty = (2*maxNoise ./ peakprom + 2*ioffset ./ peakwidth) ...
%                         .* real_ltsdvdtv';
% [err.ltsAmpErrors, outparams.lts_sweepWeights] = ...
%     compute_ltsdata_error (real_ltsv, sim_ltsv, maxNoise, isLtsError, sweepWeights);

% rmse(iSwp) = sqrt(sum((vreal - vsim).^2) / length(indFitWin));

% [err.avgSwpError, outparams] = ...
%     compute_avgerror ('cpr_avgswperr0', err.swperr, outparams, sweepWeights, totalSweepWeights);

squaredError = (vreal - vsim).^2;
meanSquaredError = nanmean(squaredError);
rmse(iSwp) = sqrt(meanSquaredError);
squaredError = (vreal - vsim).^2;
meanSquaredError = nanmean(squaredError);
rmse(iSwp) = sqrt(meanSquaredError);

% err.swperr = rmse./baseNoiseCpr';

vsimBase = cell(nSweeps, 1);
vrealBase = cell(nSweeps, 1);

baseInd = find(tvec0 >= baseWindow(1) & tvec0 <= baseWindow(2));
vsimBase{iSwp} = vsim{iSwp}(baseWindow(1):baseWindow(2))
vrealBase{iSwp} = vreal{iSwp}(baseWindow(1):baseWindow(2))

err.swperr = rmse ./ baseNoise;

if cprflag
    baseNoise = a;
end

if cprflag
    % Get the baseline indices
    baseInd = find(tFit{iSwp} >= baseWindow(1) & tFit{iSwp} <= baseWindow(2));

    % Get the traces for the baseline region
    vsimBase{iSwp} = vsim{iSwp}(baseInd);
    vrealBase{iSwp} = vreal{iSwp}(baseInd);

    % Compute the baseline rms error
    baseNoise(iSwp) = compute_rms_error(vrealBase{iSwp}, vsimBase{iSwp});
end

%% Compute sweep errors; make dimensionless by dividing by maximum noise
if cprflag
    % Use root-mean-squared error divided by the maximum noise
    err.swperr = rmse ./ baseNoiseCpr;

    % weighted-root-mean-squared-averaged cpr sweep error (dimensionless)
    [err.avgSwpError, outparams] = ...
        compute_avgerror ('cpr_avgswperr0', err.swperr, outparams, ...
                            ones(nSweeps, 1), nSweeps);
else
    % Use root-mean-squared error divided by the maximum noise
    err.swperr = rmse ./ baseNoise';

    % weighted-root-mean-squared-averaged sweep error (dimensionless)
    [err.avgSwpError, outparams] = ...
        compute_avgerror ('avgswperr0', err.swperr, outparams, sweepWeights, totalSweepWeights);
end


if ~cprflag
    maxNoise = 4 * baseNoise;
else
    maxNoise = NaN;                     % Not used but needed for parfor
end

sweepWeights = ones(nSweeps, 1)

title(['sweep ', num2str(iSwp), ': RMSE = ', num2str(err.swperr(iSwp), '%g')]);

if ~isfield(outparams, fieldnameIniterr)               % first time calculating this error

tstop = 10000;                  % Full simulation is 10 seconds == 10000 ms

if outparams.cprflag
    sweepWeights = outparams.sweepWeightsCpr;
else
    sweepWeights = outparams.sweepWeights;
end
title(['RMSE = ', num2str(swperr(iSwp), 3), '; ', ...
        'Weight = ', num2str(sweepWeights(iSwp), 3)]);

% Show limits of region fit (obsolete)
% line([xLimits(iSwp, 1) xLimits(iSwp, 1)], [-100 0], 'LineStyle', '--', 'Color', 'g');
% line([xLimits(iSwp, 2) xLimits(iSwp, 2)], [-100 0], 'LineStyle', '--', 'Color', 'g');

% Compute the residuals
residuals{iSwp} = vsim{iSwp} - vreal{iSwp};
outparams.residuals = residuals;

residuals = cellfun(@(x, y) x - y, vreal, vsim, 'UniformOutput', false);

tspace = linspace(xLimits(iSwp, 1), xLimits(iSwp, 2), length(outparams.residuals{iSwp}));
pr = plot(tspace, outparams.residuals{iSwp}, 'r');

% Note: The following doesn't work because each simDataNeuron
%           might be have different numbers of rows
%       In general, cellfun is slower than a for loop anyway
realTime = cellfun(@(x) x(:, 1), realData, 'UniformOutput', false);
simData = cellfun(@match_time_points, simDataNeuron, realTime, ...
                    'UniformOutput', false);

% Extract vectors from realData
realTime = cellfun(@(x) x(:, 1), realData, 'UniformOutput', false);
realVoltage = cellfun(@(x) x(:, 2), realData, 'UniformOutput', false);
realCurrent = cellfun(@(x) x(:, 3), realData, 'UniformOutput', false);

% Extract original vectors from simData
simTimeOrig = cellfun(@(x) x(:, 1), simData, 'UniformOutput', false);
simVoltageOrig = cellfun(@(x) x(:, 2), simData, 'UniformOutput', false);
simCurrentOrig = cellfun(@(x) x(:, 9), simData, 'UniformOutput', false);

% Interpolate model data and compute sum of squares error
parfor iSwp = 1:nSweeps
    % Interpolate model data (since CVODE (variable time step method) is applied in NEURON)
    indFitWin = find(tvec0 >= fitWindowLowerBound(iSwp) & tvec0 <= fitWindowUpperBound(iSwp));
    tFit{iSwp} = tvec0(indFitWin);                                 % time vector of interest
    [val, ind] = unique(simData{iSwp}(:, 1));                  % unique time points
    vsim{iSwp} = interp1q(val, simData{iSwp}(ind, 2), tFit{iSwp});  % simulated voltage vector of interest
    vreal{iSwp} = realData{iSwp}(indFitWin, 2);                        % recorded voltage vector of interest
    isim{iSwp} = interp1q(val, simData{iSwp}(ind, 9), tFit{iSwp});  % simulated current vector of interest
    ireal{iSwp} = realData{iSwp}(indFitWin, 3);                        % recorded current vector of interest

    % Compute root mean squared error over entire sweep region of interest
    %    Note: root-mean-squared error for each sweep (mV)
    rmse(iSwp) = compute_rms_error(vreal{iSwp}, vsim{iSwp});
end

tspace = linspace(xLimits(iSwp, 1), xLimits(iSwp, 2), length(residuals{iSwp}));

    slopeUncertainty = (2*baseNoise ./ peakprom + 2*ioffset ./ peakwidth) ...
                            .* real_ltsdvdtv';

% Save parameters as a .p file:
%    A 4 column matrix: parameter names, values, lower bounds, upper bounds
if saveParamsFlag
    pfilename = [prefix, '_neuronparams_backup.p'];
    neuronparams_backup = [outparams.neuronparamnames', num2cell(outparams.neuronparams'), ...
        num2cell(outparams.neuronparams_min'), num2cell(outparams.neuronparams_max')];
    save(fullfile(outFolder, pfilename), 'neuronparams_backup');
end

% disp(moduleLoadCommands);

% Unpack info from outparams
outFolder = outparams.outFolder;
prefix = outparams.prefix;
cprflag = outparams.cprflag;
gababsyn = outparams.gabab;
currpulseamp = outparams.currpulse;
SaveSimCmdsFlag = outparams.SaveSimCmdsFlag;
SaveSimOutFlag = outparams.SaveSimOutFlag;
CustomHoldCurrentFlag = outparams.CustomHoldCurrentFlag;
if cprflag
    simMode = 'passive';       % simulation mode: fitting of passive parameters
    holdpotential = outparams.holdpotentialCpr;
    holdcurrent = outparams.holdcurrentCpr;
    holdcurrentNoise = outparams.holdcurrentNoiseCpr;
else
    simMode = 'active';        % simulation mode: fitting of active parameters
    holdpotential = outparams.holdpotential;
    holdcurrent = outparams.holdcurrent;
    holdcurrentNoise = outparams.holdcurrentNoise;
end

%% RUN NEURON
results = cell(1, nSweeps);            % stores simulation standard outputs
tstart_NEURON = tic();
%##########
%##############
parfor iSwp = 1:nSweeps
    [status, results{iSwp}] = ...
        unix(sprintf(['%s\n', '%s - << here\n', ...
                      '%s\n', 'print "No_Errors!"\n', 'here'], ...
                     moduleLoadCommands, runNeuronCommand, simCommands{iSwp}));

    % Save simulation standard output in a text file if saveSimOutFlag is true 
    %   or if there is an error
    if saveSimOutFlag || isempty(strfind(results{iSwp}, 'No_Errors!'))
        fid = fopen(fullfile(outFolder, ...
            [prefix, '_simulation_', num2str(iSwp), '_output.txt']), 'w');
        fprintf(fid, ['Return status was: %d\n\n', ...
                      'Simulation output was:\n\n%s\n'], ...
                      status, results{iSwp});
        fclose(fid);
    end
end
%##############
%##########

%% Analyze simulation standard outputs
timeTakenNeuron = toc(tstart_NEURON);
if outparams.debugflag
    fprintf('It took %3.3g seconds to run all %d simulations with NEURON!!\n', ...
            timeTakenNeuron, nSweeps);
    fprintf('\n');
end
ranIntoErrors = isemptycell(strfind(results, 'No_Errors!'));
if sum(ranIntoErrors) == 0
    if outparams.debugflag
        fprintf('No Errors in NEURON!\n');
        fprintf('\n');

        % Display simulation outputs for debugging purposes 
        %     (these are also saved as text files if outparams.saveSimOutFlag is true)
        for iSwp = 1:nSweeps
            disp(results{iSwp});
            fprintf('\n');
        end
    end

else
    idxProblematic = find(ranIntoErrors > 0, 1);
    fprintf(['Simulation for Sweep #', num2str(idxProblematic), ...
            ' ran into errors with output:\n']);
%    error(results{idxProblematic});    %%% SOMETIMES DOESN'T PRINT STUFF FOR NO REASON!
    fprintf('%s\n', results{idxProblematic});
end

% Extract NEURON parameters to workspace
for p = 1:outparams.numparams
    % Extract from neuronparams unless the parameter is not fitted 
    %   for this column mode
    if (outparams.colmode == 1 && ~outparams.neuronparamsUseAcrossCells(p)) ...
        || (outparams.colmode == 2 && ~outparams.neuronparamsUseAcrossTrials(p))
        eval(sprintf('%s = %g;', ...
                outparams.neuronparamnames{p}, outparams.neuronparams(p)));
    end

    % Extract from neuronparamsAcrossCells if fitting across trials
    if outparams.colmode == 1 && outparams.neuronparamsUseAcrossCells(p)
        eval(sprintf('%s = %g;', outparams.neuronparamnames{p}, ...
                outparams.neuronparamsAcrossCells(p)));
    end
end

% Build simulation commands to be read by NEURON through the here-document
simCommands = cell(1, nSweeps);            % stores simulation commands
if saveSimCmdsFlag
    fid = fopen(fullfile(outFolder, ...
            [prefix, '_simulation_commands_backup.txt']), 'w');
end
for iSwp = 1:nSweeps
    % Supply NEURON parameters that were fitted across trials for this cell
    %   to workspace if fitting across cells
    if outparams.colmode == 2
        % Find the index of fitted cells for this sweep
        thisfiti = find(cell2mat(outparams.cellIdsToFit) ...
                        == outparams.cellID(iSwp));

        % Extract parameters from neuronparamsBest{thisfiti}
        for p = 1:outparams.numparams
            if outparams.neuronparamsUseAcrossTrials(p)
                eval(sprintf('%s = %g;', outparams.neuronparamnames{p}, ...
                        outparams.neuronparamsBest{thisfiti}(p)));
            end
        end
    end

    % Set output file name
    outFileName = sprintf('%s/%s_swp%d.out', outFolder, prefix, iSwp);

    % Start with the build() command in singleneuron4compgabab.hoc
    simCommands{iSwp} = sprintf('build("%s", %g, %g, %g)\n', ...
                            simMode, diamSoma, LDend, diamDend);

    % Command to adjust global passive parameters
    simCommands{iSwp} = [simCommands{iSwp}, sprintf('adjust_globalpas(%g, %g, %g)\n', ...
                                            cm, Ra, corrD)];

    % Command to adjust passive leak channel parameters
    simCommands{iSwp} = [simCommands{iSwp}, sprintf('adjust_leak(%g, %g, %g)\n', ...
                                            gpas, epas, corrD)];

    if strcmp(simMode, 'active')      % only perform these if in active fit mode
        % Command to adjust T-type calcium channel parameters
        simCommands{iSwp} = [simCommands{iSwp}, sprintf(['adjust_IT(%g, %g, ', ...
                                            '%g, %g, %g, ', ...
                                            '%g, %g, %g)\n'], ...
                                            pcabarITSoma, pcabarITDend1, ...
                                            pcabarITDend2, shiftmIT, shifthIT, ...
                                            slopemIT, slopehIT, corrD)];
        % Command to adjust H channel parameters
        simCommands{iSwp} = [simCommands{iSwp}, sprintf(['adjust_Ih(%g, %g, ', ...
                                            '%g, %g, %g, %g)\n'], ...
                                            ghbarIhSoma, ghbarIhDend1, ...
                                            ghbarIhDend2, ehIh, shiftmIh, corrD)];
        % Command to adjust inward-rectifying potassium channel parameters
        simCommands{iSwp} = [simCommands{iSwp}, sprintf(['adjust_IKir(%g, %g, ', ...
                                            '%g, %g)\n'], ...
                                            gkbarIKirSoma, gkbarIKirDend1, ...
                                            gkbarIKirDend2, corrD)];
        % Command to adjust A-Type potassium channel parameters
        simCommands{iSwp} = [simCommands{iSwp}, sprintf(['adjust_IA(%g, %g, ', ...
                                            '%g, %g)\n'], ...
                                            gkbarIASoma, gkbarIADend1, ...
                                            gkbarIADend2, corrD)];
        % Command to adjust persistent sodium channel parameters
        simCommands{iSwp} = [simCommands{iSwp}, sprintf(['adjust_INaP(%g, %g, ', ...
                                            '%g, %g)\n'], ...
                                            gnabarINaPSoma, gnabarINaPDend1, ...
                                            gnabarINaPDend2, corrD)];
    end

    % The sim() command in singleneuron4compgabab.hoc
    %   NOTE: GABA-B conductance amplitudes are converted from nS to uS
    simCommands{iSwp} = [simCommands{iSwp}, sprintf(['sim("%s", "%s", %d, %g, %g, ', ...
                                    '%g, %g, %g, %g, %g, %d, %g, %g)\n'], ...
                                    simMode, outFileName, tstop, ...
                                    holdpotential(iSwp), currpulseamp(iSwp), ...
                                    gababsyn.amp(iSwp)/1000, gababsyn.Trise(iSwp), ...
                                    gababsyn.TfallFast(iSwp), ...
                                    gababsyn.TfallSlow(iSwp), gababsyn.w(iSwp), ...
                                    customHoldCurrentFlag, holdcurrent(iSwp), ...
                                    holdcurrentNoise(iSwp))];

    % Print simulation command to a text file
    if saveSimCmdsFlag
        fprintf(fid, '%s\n', simCommands{iSwp});
    end
end
if saveSimCmdsFlag
    fclose(fid);
end

% Unpack info from outparams
gababsyn = outparams.gabab;
currpulseamp = outparams.currpulse;

% Set end time of simulations
cprflag = outparams.cprflag;
if cprflag
    simMode = 'passive';       % simulation mode: fitting of passive parameters
    holdpotential = outparams.holdpotentialCpr;
    holdcurrent = outparams.holdcurrentCpr;
    holdcurrentNoise = outparams.holdcurrentNoiseCpr;
    % Simulate until the end of the current pulse response
    tstop = outparams.cprWindow(2);
else
    simMode = 'active';        % simulation mode: fitting of active parameters
    holdpotential = outparams.holdpotential;
    holdcurrent = outparams.holdcurrent;
    holdcurrentNoise = outparams.holdcurrentNoise;
    % Simulate until the end of the IPSC response
    tstop = outparams.ipscrwin(2);    
end

% If in jitter mode and if parameter is checked, add jitter (default is 10%)
%% TODO: Need to return jittered parameters to main gui
if outparams.runmode == 3
    for iSwp = 1:nSweeps
        for p = 1:numel(outparams.neuronparams)
            if outparams.neuronparams_use(p) ~= 0        % if parameter is checked
                if outparams.neuronparamislog(p)        % These parameters might need a larger range, 
                                        % but make it the same for now
    %                outparams.neuronparams(p) = outparams.neuronparams(p) ^ ...
    %                    (1 + (outparams.neuronparams_jit(p)/100) * (-1 + 2*rand));
                    outparams.neuronparams(p) = outparams.neuronparams(p) * ...
                        (1 + (outparams.neuronparams_jit(p)/100) * (-1 + 2*rand));            
                else
                    outparams.neuronparams(p) = outparams.neuronparams(p) * ...
                        (1 + (outparams.neuronparams_jit(p)/100) * (-1 + 2*rand));            
                end
            end
        end
    end
end

% Save parameters as a .csv file:
%    A 4 column table: parameter names, values, lower bounds, upper bounds
if saveParamsFlag
    % Create path to parameters file
    sheetPath = fullfile(outFolder, [prefix, '_params.csv']);

    % Save the parameters file
    save_params(sheetPath, ...
        outparams.neuronparamnames, outparams.neuronparams, ...
        outparams.neuronparams_min, outparams.neuronparams_max);
end


% Load .out files created by NEURON
outFileName = cell(nSweeps, 1);       % .out file name for each sweep
simDataNeuron = cell(nSweeps, 1);     % simData saved by NEURON for each sweep
parfor iSwp = 1:nSweeps
    % Create the .out file name
    outFileName{iSwp} = fullfile(outFolder, ...
                        [prefix, '_swp', num2str(iSwp), '.out']);

    % Load the simData saved by NEURON into an array
    %   Note: this is an array with 12 columns
    %   TODO: log column specs here
    simDataNeuron{iSwp} = load(outFileName{iSwp});
end

residuals = cell(1, nSweeps);         % residuals (sim - real)
parfor iSwp = 1:nSweeps
    % Get the voltage vector from the simulated data
    vVecsSim = simData{iSwp}(:, 2);

    % Get the voltage vector from the real data
    vVecsRec = realData{iSwp}(:, 2);

    % Compute the residual
    residuals{iSwp} = vVecsSim - vVecsRec;
end

%       ~/Adams_Functions/m3ha_extract_and_match_time_points.m
[simData, tVecs] = ...
    m3ha_extract_and_match_time_points(realData, simDataNeuron);
% Get all time vectors from the real data
tVecs = cellfun(@(x) x(:, 1), realData, 'UniformOutput', false);
% Get the voltage vector from the real data
vVecsRec = realData{iSwp}(:, 2);x
% Get the voltage vector from the simulated data
vVecsSim = simData{iSwp}(:, 2);

[colorMap, ncg, npercg] = define_color_groups(realData, outparams);

function [colorMap, ncg, npercg] = define_color_groups (nSweeps, simMode)
%% Define the color groups

ncg = size(colorMap, 1);          % number of color groups
npercg = ceil(nSweeps/ncg);     % number of elements in a color group

% Define the color groups
[colorMap, ncg, npercg] = define_color_groups(nSweeps, outparams);

% Analyze, plot traces and display error value
[err, outparams] = errorcalc(realData, simData, outparams, ...
                                nSweeps, colorMap, ncg, npercg);

%                       Note: if 1 column, each row is a pharm condition
%                             if 2 columns, each row is a pharm-gincr condition

if strcmpi(simMode, 'passive')
    % Color groups correspond to 3 vHold conditions
    colorMap = [rgb('Blue'); rgb('Red'); rgb('Purple')];
elseif size(rowConditions, 2) == 1
    % Color groups correspond to 4 pharm conditions
    colorMap = [rgb('Black'); rgb('Blue'); ...
                rgb('Red'); rgb('Purple')];
elseif size(rowConditions, 2) == 2
    % Color groups corresponding to pharm-g incr pairs
    colorMap = colormap(jet(size(rowConditions, 1)));
end

function rowConditions = decide_from_simMode (simMode, ...
                                        rowConditions, rowConditionsCpr)
%% Decide on things based on simMode

if strcmpi(simMode, 'passive')
    if isempty(rowConditionsCpr)
        % Label the rows 1, 2, ..., nSweeps
        rowConditions = transpose(1:nSweeps);
    else
        % Use provided row conditions
        rowConditions = rowConditionsCpr;
    end
elseif strcmpi(simMode, 'active')
    if isempty(rowConditions)
        % Label the rows 1, 2, ..., nSweeps
        rowConditions = transpose(1:nSweeps);
    end
end

%% Extract info from real data
tvec0 = realData{1}(:, 1);              % time vector in ms

[err, outparams] = errorcalc(realData, simData, outparams, ...
                                nSweeps, colorMap, ncg, npercg);

xLimits = repmat(outparams.cprWindow, nSweeps, 1);
xLimits = repmat([2800, 4500], nSweeps, 1);
xLimits = repmat(outparams.ltsWindow, nSweeps, 1);
xlimitsMax = max(xLimits, [], 1);    % this is the limits used in overlapped & conductance & current figures

%    xLimits = cprWindow;
%        xLimits = ltsWindow;

% Decide on the column numbers to extract
colsToExtractReal = [TIME_COL_REC, VOLT_COL_REC, CURR_COL_REC];

% Decide on the column numbers to extract
colsToExtractSim = [TIME_COL_SIM, VOLT_COL_SIM, CURR_COL_SIM];


% Extract info from outparams
cpStartWindowOrig = outparams.cpStartWindowOrig;
vHold = outparams.vHold;
cprWindow = outparams.cprWindow;

% Update outparams
outparams.baseNoiseCpr = baseNoiseCpr;
outparams.sweepWeightsCpr = sweepWeightsCpr;
outparams.nSweeps = nSweeps;

%       ~/Adams_Functions/compute_rms_error.m

%% Recompute baseline noise and sweep weights
% Find the new baseline noise and sweep weights
% Count the number of sweeps for the current pulse response
nSweeps = numel(realData);

% Get the sampling interval in ms
simsCpr = realData{1}(2, 1) - realData{1}(1, 1);

% Find the expected current pulse start time
cpstartExpectedOrig = mean(cpStartWindowOrig);

% Find the indices for the baseline
nIndToPadBeforeCpr = round(cprWindow(1)/simsCpr);
cprbasewinLength = round(cpstartExpectedOrig/simsCpr);
indBaselineCpr = nIndToPadBeforeCpr + (1:cprbasewinLength);

% Find the holding potential, holding current and baseline noise
baseNoise = zeros(nSweeps, 1);
for iSwp = 1:nSweeps
    % Get the current voltage trace
    vRealThis = realData{iSwp}(:, 2);

    % Find the baseline noise
    baselineThis = vRealThis(indBaselineCpr);
    baseNoiseThis = compute_rms_error(baselineThis);

    % Save in arrays
    baseNoise(iSwp) = baseNoiseThis;
end

if outparams.cprflag
    baseNoise = outparams.baseNoiseCpr;
    fitWindow = outparams.fitWindowCpr;
else
    baseNoise = outparams.baseNoiseIpscr;
    fitWindow = outparams.fitWindowIpscr;
end
swperr = err.swperr;

if outparams.cprflag
    hFig.alltraces_zoom = figure(104);
else
    hFig.alltraces_zoom = figure(114);
end

% Plot traces
p1 = plot(realData{iSwp}(:, 1), realData{iSwp}(:, 2), 'k');
p2 = plot(simData{iSwp}(:, 1), simData{iSwp}(:, 2), 'r'); 
% p3 = plot(simData{iSwp}(:, 1), simData{iSwp}(:, 3), 'color', 'c');              % prox dend    %%% TO EXAMINE
% p4 = plot(simData{iSwp}(:, 1), simData{iSwp}(:, 4), 'color', [0.4 0.3 0.8])     % dist dend    %%% TO EXAMINE

% Adjust xLimits to show fitted region only
xlim(xLimits(iSwp, :));
if outparams.cprflag
    % ylim([-80 -60]);
else
    % ylim([-110 -30]);
end

if sweepWeights(iSwp) ~= 0
    text('Units', 'normalized', 'Position', [0.1 0.9], 'String', '\color{green} \bf ON');
else
    text('Units', 'normalized', 'Position', [0.1 0.9], 'String', '\color{gray} \bf OFF');
end

% Remove axes tick labels if nSweeps > 20
if nSweeps > 20
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
end

% Generate a title over all subplots
% This is from the bioinformatics toolbox
suptitle();
                    
% suplabel(sprintf('All traces for Experiment %s', strrep(outparams.prefix, '_', '\_')), 't');    
                    % This is in /home/Matlab/Downloaded_Functions

% Copy figure handle so it won't be overwritten
if outparams.cprflag
    hFig.cpr_alltraces_zoom = hFig.alltraces_zoom;
end

% Get the y-axis limits
yLimits = get(gca, 'YLim');

% Plot lines spanning the y-axis
line([fitWindow(1), fitWindow(1)], yLimits, ...
        'Color', 'g', 'LineStyle', '--');
line([fitWindow(2), fitWindow(2)], yLimits, ...
        'Color', 'g', 'LineStyle', '--');

% Re-compute baseline noise 
baseNoise = compute_baseline_noise(vVecsRec, tVecs, baseWindow);

% Re-compute sweep weights based on baseline noise
sweepWeights = 1 ./ baseNoise;

baseNoiseCprDefault = 1;        % (mV)
baseNoiseIpscrDefault = 1;      % (mV)
sweepWeightsCprDefault = 1;     % equal weights by default
sweepWeightsIpscrDefault = 1;   % equal weights by default

@(x) validateattributes(x, {'numeric'}, {'2d'}));

[tVecs, vVecsSim, iVecsSim] = ...
    extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, CURR_COL_SIM]);

% Restrict vectors to xLimits to save time on plotting
[tVecs, vVecsSim, vVecsRec, residuals, ...
    iVecsSim, vVecsDend1, vVecsDend2] = ...
    argfun(@(x) extract_subvectors(x, 'Endpoints', endPointsForPlots), ...
            tVecs, vVecsSim, vVecsRec, residuals, ...
            iVecsSim, vVecsDend1, vVecsDend2);

% Combine vectors into matrices
[tVecs, vVecsRec, vVecsSim, residuals, ...
    iVecsSim, vVecsDend1, vVecsDend2] = ...
    argfun(@(x) force_matrix(x, 'AlignMethod', 'leftAdjustPad'), ...
            tVecs, vVecsRec, vVecsSim, residuals, ...
            iVecsSim, vVecsDend1, vVecsDend2);

% CURR_COL_SIM = iCP_COL_SIM;
% CURR_COL_SIM = IDCLAMP_COL_SIM;

[realData, simData] = ...
    m3ha_average_by_vhold(realData, simData, vHold, ...
                        cpStartWindowOrig, cprWindow);

%}
