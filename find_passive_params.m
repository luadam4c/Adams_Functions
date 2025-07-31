function [passiveParams, fitResults, fitObject, ...
            goodnessOfFit, algorithmInfo, decision, ...
            tVecFitted, vVecFitted, allResults] = ...
                find_passive_params (tvec0, ivec0s, vvec0s, varargin)
%% Extract passive parameters from both the rising and falling phase of a current pulse response
% Usage: [passiveParams, fitResults, fitObject, ...
%           goodnessOfFit, algorithmInfo, decision, ...
%           tVecFitted, vVecFitted, allResults] = ...
%               find_passive_params (tvec0, ivec0s, vvec0s, varargin)
% Outputs:
%       passiveParams   - passive parameters returned by 
%                           fit_and_estimate_passive_params.m
%                           for the phase given by decision
%                         as well as the following estimates 
%                           from I-V relationships:
%                           epasEstimate: resting membrane potential (mV)
%                           RinEstimate: input resistance (MOhm)
%       fitResults      - fitting results returned by 
%                           fit_and_estimate_passive_params.m
%                           for the phase given by decision
%       fitObject       - cfit object returned by 
%                           fit_and_estimate_passive_params.m
%                           for the phase given by decision
%       goodnessOfFit   - goodness of fit structure returned by
%                           fit_and_estimate_passive_params.m
%                           for the phase given by decision
%       algorithmInfo   - algorithm info structure returned by
%                           fit_and_estimate_passive_params.m
%                           for the phase given by decision
%       decision        - a string of what phase/averaging method combination
%                           was used for fitting and parameter estimation
%       tVecFitted      - time vector(s) that was used for the curve fit
%       vVecFitted      - voltage vector(s) that was used for the curve fit
%       allResults  - a structure of other results along the way, with fields:
%                       pulseAmplitude - vector of current pulse amplitudes (pA)
%                                           for each sweep
%                       pulseWidth     - vector of pulse widths (ms) 
%                                           for each sweep
%                       voltageChange  - vector of overall changes in 
%                                           membrane potential recorded (mV) 
%                                           for each sweep
%                       rmseRising     - vector of root-mean-squared errors (mV) 
%                                           in the rising phase for each sweep
%                       rmseFalling    - vector of root-mean-squared errors (mV) 
%                                           in the falling phase for each sweep
%                     and all outputs from different methods
%                       
% Arguments:
%       tvec0       - original time vector (ms)
%                   must be a numeric vector
%       ivec0s      - original current vector(s) (pA)
%                   must be a numeric array with the same # of rows as tvec0
%       vvec0s      - original voltage vector(s) (mV)
%                   must be a numeric array with the same # of rows as tvec0
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'HoldCurrent': holding current (pA) applied for each trace
%                   must be a numeric vector
%                   default == []
%                   - 'PulseWindow': window (ms) in which the current 
%                                       pulse would lie, e.g. [95, 115]
%                   must be within range of [timeBase timeFinal]
%                   default == [timeBase, tvec0(end)]
%                   - 'PulseResponseWindow': window (ms) in which the 
%                                           current pulse response would lie, 
%                                           e.g. [95, 500]
%                   must be within range of [timeBase timeFinal]
%                   default == [timeBase, tvec0(end)]
%                   - 'PlotFlag': whether to plot figures
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OutFolder': directory to place outputs
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FileBase': base of output file names (without extension)
%                   must be a string scalar or a character vector
%                   default == 'unnamed'
%                   - 'Ivec1s': median-filtered current vector(s) in pA
%                   must be a numeric array with the same # of rows as tvec0
%                   default == medfilt1(ivec0s, round(medFiltWindowMs/siMs))
%                   - 'Suffix': suffix for file and directory names
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'TitleMod': modifier for titles
%                   must be a string scalar or a character vector
%                   default == ''
%
% Requires:
%       cd/argfun.m
%       cd/check_subdir.m
%       cd/compute_average_trace.m
%       cd/estimate_resting_potential.m
%       cd/extract_subvectors.m
%       cd/find_window_endpoints.m
%       cd/fit_and_estimate_passive_params.m
%       cd/isnumericvector.m
%       cd/locate_functionsdir.m
%       cd/medianfilter.m
%       cd/parse_pulse.m
%       cd/parse_pulse_response.m
%       cd/plot_ball_stick.m
%       cd/plot_cfit_pulse_response.m
%       cd/plot_pulse.m
%       cd/plot_pulse_response.m
%       cd/print_structure.m
%       cd/suptitle.m (same as built-in but with fixed typo)
%       ~/Downloaded_Functions/subplotsqueeze.m
%
% Used by:
%       cd/m3ha_parse_dclamp_data.m
%       cd/m3ha_plot_figure03.m
%       cd/m3ha_estimate_passive_params.m

% File History:
% 2016-10-27 - Adapted from m3ha_find_ipsc_peak.m and dclampDataExtractor.m
% 2016-10-31 - Changed equation form of aFittype to 
%                   'a*(exp(-x/b)-1)+c*(exp(-x/d)-1)' 
%                   from 'a*exp(-x/b)+c*exp(-x/d)+e'
% 2016-10-31 - Removed '- round(0.5/siMs)' from the definition of indBaseline
% 2016-10-31 - Made sure tau1 >= tau2
% 2016-11-01 - Added dataMode so that the titles and filenames can change
% 2016-11-01 - Exclude faulty traces, including those with spontaneous spikes, 
%                   from the fitting
% 2016-11-01 - Added pulse width
% 2016-11-01 - Added long pulse response
% 2016-11-01 - Added L, rho, taum
% 2016-11-01 - Plot ivec1s only (ivec0s is not directly used in the analyses)
% 2016-11-02 - Moved check directories to check_subdir.m
% 2016-11-02 - Added the falling phase of the current pulse response
% 2016-11-12 - Now returns estimates from pooled data if those from averaged 
%               data don't exist
% 2016-11-12 - Now outputs all estimates
% 2016-12-04 - Added series resistance Rs and changed the way somatic and 
%               dendritic resistances are computed
% 2016-12-04 - Increase the upper bound of tau to 500 ms
% 2016-12-04 - Added rmseRising && rmseFalling, the root-mean-squared errors of 
%               the rising and falling phase, respectively, for each sweep
% 2016-12-04 - Added the functions fit_setup_2exp && measure_error
% 2016-12-04 - Fixed measure_error so that sweeps yielding nonsensical 
%               responses have rmse = Inf
% 2016-12-05 - Added tau1Range = [20, 500] and tau2Range = [0, 20]
% 2016-12-05 - Removed typtau and tau_max, changed initial conditions 
%               to tau1Range(1) and tau2Range(2)
% 2016-12-05 - Added tau1RangeRising and tau2RangeRising
% 2016-12-05 - Added tau1Init, tau2Init, tau1InitRising, tau2InitRising
% 2017-12-21 - SpecsForFitmode() -> specs_for_fitmode()
% 2018-01-24 - Added isdeployed
% 2018-10-03 - Changed tabs to spaces
% 2018-10-10 - Changed voltageChange to be steady - base
% 2018-10-12 - Moved code to its own functions
% 2018-10-14 - Added the combined phase
% 2018-10-15 - Removed dependence to dataMode and added Suffix and TitleMod
%               as optional arguments instead
% 2018-11-13 - Updated usage of parse_pulse_response.m
% 2019-12-20 - Added tVecFitted, vVecFitted as outputs

%% Flags

% TODO: Make optional argument
printFieldsFlag = 0; %1;

%% Hard-coded parameters
% Parameters to be consistent with m3ha_import_raw_traces.m
meanVoltageWindow = 0.5;    % width in ms for calculating mean voltage 
                            %   for input resistance calculations

% Parameters used for data analysis
medFiltWindowMs = 10;       % width in ms for the median filter for corrupted data 
                            %   (current traces)
spikeThresholdInit = -45;   % initial amplitude threshold (mV) for detecting a spike
maxScalingFactor = 10;      % maximum scaling factor from initial voltage 
                            %   amplitude value
%{
tau1InitRising = 23;        % typical tau1 in rising phase ~ 23 ms
tau1RangeRising = [7, 200]; % range of tau1 in rising phase
tau2InitRising = 0.8;       % typical tau2 in rising phase ~ 0.8 ms
tau2RangeRising = [0, 7];   % range of tau2 in rising phase
tau1InitFalling = 50;           % typical tau1 in falling phase ~ 50 ms
tau1RangeFalling = [20, 500];   % range of tau1 in falling phase
tau2InitFalling = 6;            % typical tau2 in falling phase ~ 6 ms
tau2RangeFalling = [0, 20];     % range of tau2 in falling phase
%}

tau1InitRising = 10;            % typical tau1 in rising phase 
tau1RangeRising = [0, Inf];     % range of tau1 in rising phase
tau2InitRising = 1;             % typical tau2 in rising phase
tau2RangeRising = [0, Inf];     % range of tau2 in rising phase
tau1InitFalling = 10;           % typical tau1 in falling phase 
tau1RangeFalling = [0, Inf];    % range of tau1 in falling phase
tau2InitFalling = 1;            % typical tau2 in falling phase
tau2RangeFalling = [0, Inf];    % range of tau2 in falling phase
tau1InitCombined = 10;          % typical tau1 in combined phase 
tau1RangeCombined = [0, Inf];   % range of tau1 in combined phase
tau2InitCombined = 1;           % typical tau2 in combined phase
tau2RangeCombined = [0, Inf];   % range of tau2 in combined phase

%% Fixed parameters for all cells
Cm = 0.88;              % specific membrane capacitance (uF/cm^2)
Ra = 173;               % axial resistivity (Ohm-cm)
Rs = 10;                % series (pipette) resistance (MOhm)

%% Directories for placing figures
directories = {'cprRising', 'cprFalling', 'cprCombined', 'passive'};

%% Default values for optional arguments
verboseDefault = false;             % don't print to standard output by default
holdCurrentDefault = [];            % no holding current supplied by default
pulseWindowDefault = [];            % set later
pulseResponseWindowDefault = [];    % set later
plotFlagDefault = false;
outFolderDefault = pwd;
fileBaseDefault = 'unnamed';
ivec1sDefault = [];                 % set later
suffixDefault = '';
titleModDefault = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% If not compiled, add directories to search path for required functions
if exist('subplotsqueeze.m', 'file') ~= 2 && ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for subplotsqueeze.m
    addpath_custom(fullfile(functionsDirectory, 'Downloaded_Functions')); 
end

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'tvec0', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'ivec0s', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addRequired(iP, 'vvec0s', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'HoldCurrent', holdCurrentDefault, ...
    @(x) assert(isnumericvector(x), 'HoldCurrent must be a numeric vector!'));
addParameter(iP, 'PulseWindow', pulseWindowDefault, ...
    @(x) isempty(x) || isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'PulseResponseWindow', pulseResponseWindowDefault, ...
    @(x) isempty(x) || isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'PlotFlag', plotFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileBase', fileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Ivec1s', ivec1sDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'Suffix', suffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TitleMod', titleModDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, tvec0, ivec0s, vvec0s, varargin{:});
verbose = iP.Results.Verbose;
holdCurrent = iP.Results.HoldCurrent;
pulseWindow = iP.Results.PulseWindow;
pulseResponseWindow = iP.Results.PulseResponseWindow;
plotFlag = iP.Results.PlotFlag;
outFolder = iP.Results.OutFolder;
fileBase = iP.Results.FileBase;
ivec1s = iP.Results.Ivec1s;
suffix = iP.Results.Suffix;
titleMod = iP.Results.TitleMod;

% Check relationships between arguments
if ~isequal(length(tvec0), size(ivec0s, 1))
    error('Time and current vectors do not have the same length!');
elseif ~isequal(length(tvec0), size(vvec0s, 1))
    error('Time and voltage vectors do not have the same length!');
elseif ~isempty(pulseWindow) && length(pulseWindow) < 2
    error('PulseWindow must have a start and an end!');
elseif ~isempty(pulseResponseWindow) && length(pulseResponseWindow) < 2
    error('PulseResponseWindow must have a start and an end!');
elseif ~isempty(ivec1s) && ~isequal(size(ivec1s, 1), size(tvec0, 1))
    error('ivec1s must have the same column length as tvec0!');
end

%% Preparation
% Initialize outputs
passiveParams = [];
fitResults = [];
fitObject = [];
goodnessOfFit = [];
algorithmInfo = [];
decision = 'none';
allResults = [];

% Force tvec0 to be a column vector
tvec0 = tvec0(:);

% Get the sampling interval in ms
siMs = tvec0(2) - tvec0(1);

% Get the number of sweeps
nSweeps = size(ivec0s, 2);

% Get the time just before the start
timeBase = tvec0(1) - siMs;

% Get the time just beyond the end
% Note: in case siMs is not divisible by pulseResponseWindow(2)
timeFinal = tvec0(end) + siMs;

% Check value ranges for the windows
if ~isempty(pulseWindow) && pulseWindow(1) < timeBase
    error('pulseWindow(1) == %g < timeBase == %g is out of range!', ...
            pulseWindow(1), timeBase);
elseif ~isempty(pulseWindow) && pulseWindow(2) > timeFinal
    error('pulseWindow(2) == %g > timeFinal == %g is out of range!', ...
            pulseWindow(2), timeFinal);
elseif ~isempty(pulseResponseWindow) && pulseResponseWindow(1) < timeBase
    error('pulseResponseWindow(1) == %g < timeBase == %g is out of range!', ...
            pulseResponseWindow(1), timeBase);
elseif ~isempty(pulseResponseWindow) && pulseResponseWindow(2) > timeFinal
    error('pulseResponseWindow(2) == %g > timeFinal == %g is out of range!', ...
            pulseResponseWindow(2), timeFinal);
end

% Set defaults for the windows
if isempty(pulseWindow)
    pulseWindow = [timeBase, tvec0(end)];
end
if isempty(pulseResponseWindow)
    pulseResponseWindow = [timeBase, tvec0(end)];
end

% If not provided, median filter current traces to get rid of corrupted data
if isempty(ivec1s)
    ivec1s = medianfilter(ivec0s, medianFilterWindowMs, siMs);
end

% Append suffix to directory names
directories = strcat(directories, suffix);

% Print to standard output
if verbose
    fprintf('ANALYZING passive parameters for %s ...\n', fileBase);
    fprintf('Sampling interval == %g ms\n', siMs);
    fprintf('\n');
end

%% Restrict the original vectors to the current pulse response window for speed
% Find the indices for the current pulse response
endPoints = find_window_endpoints(pulseResponseWindow, tvec0);

% Extract the regions to fit
[tvecCpr, vvecsCpr, ivecsCpr] = ...
    argfun(@(x) extract_subvectors(x, 'Endpoints', endPoints), ...
            tvec0, vvec0s, ivec1s);

%% Parse the pulses and pulse responses
% Parse the pulse
[pulseParams, ~] = parse_pulse(ivecsCpr);

% Extract the current pulse widths and amplitudes
pulseWidth = pulseParams.pulseWidthSamples * siMs;
pulseAmplitude = pulseParams.pulseAmplitude;

% Parse the pulse response, using the indices from the pulse
[responseParams, responseData] = ...
    parse_pulse_response(vvecsCpr, siMs, 'PulseVector', ivecsCpr, ...
                            'SameAsPulse', true, 'FitResponse', false, ...
                            'MeanValueWindowMs', meanVoltageWindow);

% Extract the holding potentials (mV)
holdPotential = responseParams.baseValue;

% Extract the maximum voltage values (mV)
maxVoltage = responseParams.maxValue;

% Extract the voltage changes (mV)
voltageChange = responseParams.steadyAmplitude;

% Decide whether each trace will be used
toUse = pulseWidth >= 0 & ...
        sign(pulseAmplitude) == sign(voltageChange) & ...
        maxVoltage <= spikeThresholdInit;

% Count the number of traces to be used
nToUse = sum(toUse);

% If no trace used, return
if nToUse == 0
    fprintf('No trace matches the 3 criteria below:\n');
    fprintf('\t1. pulseWidth >= 0\n');
    fprintf('\t2. sign(pulseAmplitude) == sign(voltageChange)\n');
    fprintf('\t3. maxVoltage > %g\n', spikeThresholdInit);
    return
end

%% Compute average values
% Compute average pulse width (ms)
meanPulseWidth = mean(pulseWidth(toUse));

% Compute average pulse amplitude (pA)
meanPulseAmplitude = mean(pulseAmplitude(toUse));

% Compute average recorded voltage change (mV)
meanVoltageChange = mean(voltageChange(toUse));

if verbose
    fprintf('Average recorded voltage change (mV) is %g\n', meanVoltageChange);
    fprintf('Average current pulse amplitude (pA) is %g\n', meanPulseAmplitude);
    fprintf('Average pulse width (ms) is %g\n', meanPulseWidth);
    fprintf('\n');
end

%% Restrict to vectors to use 
% Extract shifted rising/falling phase vectors
tvecsRising = responseData.tvecsRising;
vvecsRising = responseData.vvecsRising;
tvecsFalling = responseData.tvecsFalling;
vvecsFalling = responseData.vvecsFalling;
tvecsCombined = responseData.tvecsCombined;
vvecsCombined = responseData.vvecsCombined;

% Restrict to vectors to use
tvecsRisingToUse = tvecsRising(toUse);
vvecsRisingToUse = vvecsRising(toUse);
tvecsFallingToUse = tvecsFalling(toUse);
vvecsFallingToUse = vvecsFalling(toUse);
tvecsCombinedToUse = tvecsCombined(toUse);
vvecsCombinedToUse = vvecsCombined(toUse);

%% Estimate the resting membrane potential
if ~isempty(holdCurrent)
    % Restrict computed values to vectors to use
    holdCurrentToUse = holdCurrent(toUse);
    holdPotentialToUse = holdPotential(toUse);

    % Estimate the resting membrane potential (mV) 
    %   and the input resistance (MOhm)
    [epasEstimate, RinEstimate] = ...
        estimate_resting_potential(holdPotentialToUse, holdCurrentToUse);
end

%% Put data from all current pulse responses together in two different ways
% Method 1: Put all data points together
tAllRising = vertcat(tvecsRisingToUse{:});
vAllRising = vertcat(vvecsRisingToUse{:});
tAllFalling = vertcat(tvecsFallingToUse{:});
vAllFalling = vertcat(vvecsFallingToUse{:});
tAllCombined = vertcat(tvecsCombinedToUse{:});
vAllCombined = vertcat(vvecsCombinedToUse{:});

% Method 2: Average all current pulse responses
% Compute the average voltage traces, truncating and aligning to the left
vAvgRising = compute_average_trace(vvecsRisingToUse, ...
                                    'AlignMethod', 'leftadjust');
vAvgFalling = compute_average_trace(vvecsFallingToUse, ...
                                    'AlignMethod', 'leftadjust');
vAvgCombined = compute_average_trace(vvecsCombinedToUse, ...
                                    'AlignMethod', 'leftadjust');

% Count the number of samples in each average trace
nSamplesRising = length(vAvgRising);
nSamplesFalling = length(vAvgFalling);
nSamplesCombined = length(vAvgCombined);

% Truncate the first time vector accordingly
tAvgRising = tvecsRisingToUse{1}(1:nSamplesRising);
tAvgFalling = tvecsFallingToUse{1}(1:nSamplesFalling);
tAvgCombined = tvecsCombinedToUse{1}(1:nSamplesCombined);

%% Fit a double exponential to each voltage trace
% Fit to the rising phase of pooled data
[paramsAllRising, fitResultsAllRising, cfitAllRising, ...
    gofAllRising, algInfoAllRising] = ...
    fit_and_estimate_passive_params(tAllRising, vAllRising, ...
                            meanPulseWidth, meanPulseAmplitude, ...
                            'PhaseName', 'rising', ...
                            'AmplitudeEstimate', meanVoltageChange, ...
                            'MaxScalingFactor', maxScalingFactor, ...
                            'Tau1Init', tau1InitRising, ...
                            'Tau2Init', tau2InitRising, ...
                            'Tau1Range', tau1RangeRising, ...
                            'Tau2Range', tau2RangeRising, ...
                            'MembraneCapacitance', Cm, ...
                            'AxialResistivity', Ra, ...
                            'SeriesResistance', Rs);

% Print the structure if requested
if printFieldsFlag
    print_structure(paramsAllRising);
end

% Fit to the falling phase of pooled data
[paramsAllFalling, fitResultsAllFalling, cfitAllFalling, ...
    gofAllFalling, algInfoAllFalling] = ...
    fit_and_estimate_passive_params(tAllFalling, vAllFalling, ...
                            meanPulseWidth, meanPulseAmplitude, ...
                            'PhaseName', 'falling', ...
                            'AmplitudeEstimate', meanVoltageChange, ...
                            'MaxScalingFactor', maxScalingFactor, ...
                            'Tau1Init', tau1InitFalling, ...
                            'Tau2Init', tau2InitFalling, ...
                            'Tau1Range', tau1RangeFalling, ...
                            'Tau2Range', tau2RangeFalling, ...
                            'MembraneCapacitance', Cm, ...
                            'AxialResistivity', Ra, ...
                            'SeriesResistance', Rs);

% Print the structure if requested
if printFieldsFlag
    print_structure(paramsAllFalling);
end

% Fit to the combined rising and falling phase of pooled data
[paramsAllCombined, fitResultsAllCombined, cfitAllCombined, ...
    gofAllCombined, algInfoAllCombined] = ...
    fit_and_estimate_passive_params(tAllCombined, vAllCombined, ...
                            meanPulseWidth, meanPulseAmplitude, ...
                            'PhaseName', 'combined', ...
                            'AmplitudeEstimate', meanVoltageChange, ...
                            'MaxScalingFactor', maxScalingFactor, ...
                            'Tau1Init', tau1InitCombined, ...
                            'Tau2Init', tau2InitCombined, ...
                            'Tau1Range', tau1RangeCombined, ...
                            'Tau2Range', tau2RangeCombined, ...
                            'MembraneCapacitance', Cm, ...
                            'AxialResistivity', Ra, ...
                            'SeriesResistance', Rs);

% Print the structure if requested
if printFieldsFlag
    print_structure(paramsAllCombined);
end

% Fit to the rising phase of average trace
[paramsAvgRising, fitResultsAvgRising, cfitAvgRising, ...
    gofAvgRising, algInfoAvgRising] = ...
    fit_and_estimate_passive_params(tAvgRising, vAvgRising, ...
                            meanPulseWidth, meanPulseAmplitude, ...
                            'PhaseName', 'rising', ...
                            'AmplitudeEstimate', meanVoltageChange, ...
                            'MaxScalingFactor', maxScalingFactor, ...
                            'Tau1Init', tau1InitRising, ...
                            'Tau2Init', tau2InitRising, ...
                            'Tau1Range', tau1RangeRising, ...
                            'Tau2Range', tau2RangeRising, ...
                            'MembraneCapacitance', Cm, ...
                            'AxialResistivity', Ra, ...
                            'SeriesResistance', Rs);

% Print the structure if requested
if printFieldsFlag
    print_structure(paramsAvgRising);
end

% Fit to the falling phase of average trace
[paramsAvgFalling, fitResultsAvgFalling, cfitAvgFalling, ...
    gofAvgFalling, algInfoAvgFalling] = ...
    fit_and_estimate_passive_params(tAvgFalling, vAvgFalling, ...
                            meanPulseWidth, meanPulseAmplitude, ...
                            'PhaseName', 'falling', ...
                            'AmplitudeEstimate', meanVoltageChange, ...
                            'MaxScalingFactor', maxScalingFactor, ...
                            'Tau1Init', tau1InitFalling, ...
                            'Tau2Init', tau2InitFalling, ...
                            'Tau1Range', tau1RangeFalling, ...
                            'Tau2Range', tau2RangeFalling, ...
                            'MembraneCapacitance', Cm, ...
                            'AxialResistivity', Ra, ...
                            'SeriesResistance', Rs);

% Print the structure if requested
if printFieldsFlag
    print_structure(paramsAvgFalling);
end

% Fit to the combined rising and falling phase of average trace
[paramsAvgCombined, fitResultsAvgCombined, cfitAvgCombined, ...
    gofAvgCombined, algInfoAvgCombined] = ...
    fit_and_estimate_passive_params(tAvgCombined, vAvgCombined, ...
                            meanPulseWidth, meanPulseAmplitude, ...
                            'PhaseName', 'combined', ...
                            'AmplitudeEstimate', meanVoltageChange, ...
                            'MaxScalingFactor', maxScalingFactor, ...
                            'Tau1Init', tau1InitCombined, ...
                            'Tau2Init', tau2InitCombined, ...
                            'Tau1Range', tau1RangeCombined, ...
                            'Tau2Range', tau2RangeCombined, ...
                            'MembraneCapacitance', Cm, ...
                            'AxialResistivity', Ra, ...
                            'SeriesResistance', Rs);

% Print the structure if requested
if printFieldsFlag
    print_structure(paramsAvgCombined);
end

%% Decide on final parameters in the following order: 
%       1. combined phase of average trace
%       2. combined phase of pooled data
%       3. falling phase of average trace
%       4. falling phase of pooled data
%       5. rising phase of average trace
%       6. rising phase of pooled data
if isstruct(paramsAvgCombined)
    tVecFitted = tAvgCombined;
    vVecFitted = vAvgCombined;
    passiveParams = paramsAvgCombined;
    fitResults = fitResultsAvgCombined;
    fitObject = cfitAvgCombined;
    goodnessOfFit = gofAvgCombined;
    algorithmInfo = algInfoAvgCombined;
    decision = 'combined phase of average trace';
elseif isstruct(paramsAllCombined)
    tVecFitted = tAllCombined;
    vVecFitted = vAllCombined;
    passiveParams = paramsAllCombined;
    fitResults = fitResultsAllCombined;
    fitObject = cfitAllCombined;
    goodnessOfFit = gofAllCombined;
    algorithmInfo = algInfoAllCombined;
    decision = 'combined phase of pooled data';
elseif isstruct(paramsAvgFalling)
    tVecFitted = tAvgFalling;
    vVecFitted = vAvgFalling;
    passiveParams = paramsAvgFalling;
    fitResults = fitResultsAvgFalling;
    fitObject = cfitAvgFalling;
    goodnessOfFit = gofAvgFalling;
    algorithmInfo = algInfoAvgFalling;
    decision = 'falling phase of average trace';
elseif isstruct(paramsAllFalling)
    tVecFitted = tAllFalling;
    vVecFitted = vAllFalling;
    passiveParams = paramsAllFalling;
    fitResults = fitResultsAllFalling;
    fitObject = cfitAllFalling;
    goodnessOfFit = gofAllFalling;
    algorithmInfo = algInfoAllFalling;
    decision = 'falling phase of pooled data';
elseif isstruct(paramsAvgRising)
    tVecFitted = tAvgRising;
    vVecFitted = vAvgRising;
    passiveParams = paramsAvgRising;
    fitResults = fitResultsAvgRising;
    fitObject = cfitAvgRising;
    goodnessOfFit = gofAvgRising;
    algorithmInfo = algInfoAvgRising;
    decision = 'rising phase of average trace';
elseif isstruct(paramsAllRising)
    tVecFitted = tAllRising;
    vVecFitted = vAllRising;
    passiveParams = paramsAllRising;
    fitResults = fitResultsAllRising;
    fitObject = cfitAllRising;
    goodnessOfFit = gofAllRising;
    algorithmInfo = algInfoAllRising;
    decision = 'rising phase of pooled data';
end    

% Add estimates of epas and Rin if holding currents provided
if ~isempty(holdCurrent)
    passiveParams.epasEstimate = epasEstimate;
    passiveParams.RinEstimate = RinEstimate;
end

%% Find a goodness-of-fit measure (root-mean-squared error) for each sweep
%{
rmseRising = measure_error(tvecsRising, vvecsRising, 'rising', ...
                            meanVoltageChange, maxScalingFactor, ...
                            tau1InitRising, tau2InitRising, ...
                            tau1RangeRising, tau2RangeRising, meanPulseWidth);
rmseFalling = measure_error(tvecsFalling, vvecsFalling, 'falling', ...
                            meanVoltageChange, maxScalingFactor, ...
                            tau1InitFalling, tau2InitFalling, ...
                            tau1RangeFalling, tau2RangeFalling, meanPulseWidth);
%}

%% Plot results
if plotFlag
    % Check if needed directories exist in outFolder
    check_subdir(outFolder, directories);

    % TODO: make three functions
    % plot_cprRising;
    % plot_cprFalling;
    % plot_passive;

    % Plot rising phase of current pulse response
    h = figure(2000);
    set(h, 'Visible', 'off');
    set(h, 'Name', 'Rising phase of current pulse response');
    clf(h);
    subplot(2,2,1);
    plot_pulse(tvecCpr, ivecsCpr, 'PulseParams', pulseParams, ...
                                'ToUse', toUse, 'XLimits', pulseWindow);
    subplot(2,2,2);
    plot_pulse_response(tvecCpr, vvecsCpr, 'ResponseParams', responseParams, ...
                                'ToUse', toUse, 'XLimits', pulseWindow);
    subplot(2,2,3);
    plot_cfit_pulse_response(tAllRising, vAllRising, ...
                            'FitObject', cfitAllRising, ...
                            'FitResults', fitResultsAllRising, ...
                            'GoodnessOfFit', gofAllRising, ...
                            'PassiveParams', paramsAllRising, ...
                            'LegendLocation', 'suppress');
    subplot(2,2,4);
    plot_cfit_pulse_response(tAvgRising, vAvgRising, ...
                            'FitObject', cfitAvgRising, ...
                            'FitResults', fitResultsAvgRising, ...
                            'GoodnessOfFit', gofAvgRising, ...
                            'PassiveParams', paramsAvgRising, ...
                            'LegendLocation', 'suppress');
    figname = fullfile(outFolder, directories{1}, ...
                        [fileBase, '_cprRising', suffix, '.png']);
    % subplotsqueeze(h, 1.1);
    % suptitle(strjoin({'Rising phase of current pulse response for', ...
    %                     strrep(fileBase, '_', '\_'), titleMod}));
    saveas(h, figname);
    % close(h);

    % Plot falling phase of current pulse response
    h = figure(2001);
    set(h, 'Visible', 'off');
    set(h, 'Name', 'Falling phase of current pulse response');
    clf(h);
    subplot(2,2,1);
    plot_pulse(tvecCpr, ivecsCpr, 'PulseParams', pulseParams, ...
                                'ToUse', toUse, 'XLimits', pulseResponseWindow);
    subplot(2,2,2);
    plot_pulse_response(tvecCpr, vvecsCpr, 'ResponseParams', responseParams, ...
                                'ToUse', toUse, 'XLimits', pulseResponseWindow);
    subplot(2,2,3);
    plot_cfit_pulse_response(tAllFalling, vAllFalling, ...
                            'FitObject', cfitAllFalling, ...
                            'FitResults', fitResultsAllFalling, ...
                            'GoodnessOfFit', gofAllFalling, ...
                            'PassiveParams', paramsAllFalling, ...
                            'LegendLocation', 'suppress');
    subplot(2,2,4);
    plot_cfit_pulse_response(tAvgFalling, vAvgFalling, ...
                            'FitObject', cfitAvgFalling, ...
                            'FitResults', fitResultsAvgFalling, ...
                            'GoodnessOfFit', gofAvgFalling, ...
                            'PassiveParams', paramsAvgFalling, ...
                            'LegendLocation', 'suppress');
    figname = fullfile(outFolder, directories{2}, ...
                        [fileBase, '_cprFalling', suffix, '.png']);
    % subplotsqueeze(h, 1.1);
    % suptitle(strjoin({'Falling phase of current pulse response for', ...
    %                     strrep(fileBase, '_', '\_'), titleMod}));
    saveas(h, figname);
    % close(h);

    % Plot combined phases of current pulse response
    h = figure(2002);
    set(h, 'Visible', 'off');
    set(h, 'Name', 'Combined phases of current pulse response');
    clf(h);
    subplot(2,2,1);
    plot_pulse(tvecCpr, ivecsCpr, 'PulseParams', pulseParams, ...
                                'ToUse', toUse, 'XLimits', pulseResponseWindow);
    subplot(2,2,2);
    plot_pulse_response(tvecCpr, vvecsCpr, 'ResponseParams', responseParams, ...
                                'ToUse', toUse, 'XLimits', pulseResponseWindow);
    subplot(2,2,3);
    plot_cfit_pulse_response(tAllCombined, vAllCombined, ...
                            'FitObject', cfitAllCombined, ...
                            'FitResults', fitResultsAllCombined, ...
                            'GoodnessOfFit', gofAllCombined, ...
                            'PassiveParams', paramsAllCombined, ...
                            'LegendLocation', 'suppress');
    subplot(2,2,4);
    plot_cfit_pulse_response(tAvgCombined, vAvgCombined, ...
                            'FitObject', cfitAvgCombined, ...
                            'FitResults', fitResultsAvgCombined, ...
                            'GoodnessOfFit', gofAvgCombined, ...
                            'PassiveParams', paramsAvgCombined, ...
                            'LegendLocation', 'suppress');
    figname = fullfile(outFolder, directories{3}, ...
                        [fileBase, '_cprCombined', suffix, '.png']);
    subplotsqueeze(h, 1.1);
    suptitle(strjoin({'Combined phases of current pulse response for', ...
                        strrep(fileBase, '_', '\_'), titleMod}));
    saveas(h, figname);
    % close(h);

    % Plot input resistance analysis with passive parameter fitting
    scrsz = get(groot, 'ScreenSize');        % screen size
    h = figure(2003);
    set(h, 'Position', [1, scrsz(4)/4, scrsz(3)/2, scrsz(4)*3/4]);
    set(h, 'PaperPositionMode', 'auto')
    set(h, 'Visible', 'off');
    set(h, 'Name', 'Passive parameter fitting');
    clf(h);
    subplot(3,2,1);
    plot_cfit_pulse_response(tAllCombined, vAllCombined, ...
                            'FitObject', cfitAllCombined, ...
                            'FitResults', fitResultsAllCombined, ...
                            'GoodnessOfFit', gofAllCombined, ...
                            'PassiveParams', paramsAllCombined, ...
                            'LegendLocation', 'suppress');
    text(0.1, 1.1, sprintf('Total number of sweeps: %d', nSweeps), ...
        'Units', 'normalized');
    subplot(3,2,2);
    plot_cfit_pulse_response(tAvgCombined, vAvgCombined, ...
                            'FitObject', cfitAvgCombined, ...
                            'FitResults', fitResultsAvgCombined, ...
                            'GoodnessOfFit', gofAvgCombined, ...
                            'PassiveParams', paramsAvgCombined, ...
                            'LegendLocation', 'suppress');
    subplot(3,2,3);
    plot_cfit_pulse_response(tAllFalling, vAllFalling, ...
                            'FitObject', cfitAllFalling, ...
                            'FitResults', fitResultsAllFalling, ...
                            'GoodnessOfFit', gofAllFalling, ...
                            'PassiveParams', paramsAllFalling, ...
                            'LegendLocation', 'suppress');
    subplot(3,2,4);
    plot_cfit_pulse_response(tAvgFalling, vAvgFalling, ...
                            'FitObject', cfitAvgFalling, ...
                            'FitResults', fitResultsAvgFalling, ...
                            'GoodnessOfFit', gofAvgFalling, ...
                            'PassiveParams', paramsAvgFalling, ...
                            'LegendLocation', 'suppress');
    subplot(3,2,5);
    plot_geometry_and_passive_params(paramsAllCombined, paramsAllFalling);
    % plot_geometry_and_passive_params(paramsAllCombined, []);
    subplot(3,2,6);
    plot_geometry_and_passive_params(paramsAvgCombined, paramsAvgFalling);
    % plot_geometry_and_passive_params(paramsAvgCombined, []);
    % figname = fullfile(outFolder, directories{4}, ...
    %                       [fileBase, '_passive', suffix, '.png']);
    figNameNoPng = fullfile(outFolder, directories{4}, ...
                            [strrep(fileBase, '.', 'p'), '_passive', suffix]);
    % subplotsqueeze(h, 1.1);
    suptitle(sprintf('Passive parameter fitting for %s %s\n', ...
                        strrep(fileBase, '_', '\_'), titleMod));
    print(h, figNameNoPng, '-dpng', '-r0');
    % saveas(h, figname);
    % close(h);
end

%% Output results
% Store fixed parameters in algorithmInfo
algorithmInfo.medFiltWindowMs = medFiltWindowMs;
algorithmInfo.spikeThresholdInit = spikeThresholdInit;
algorithmInfo.maxScalingFactor = maxScalingFactor;
algorithmInfo.tau1InitRising = tau1InitRising;
algorithmInfo.tau1RangeRising = tau1RangeRising;
algorithmInfo.tau2InitRising = tau2InitRising;
algorithmInfo.tau2RangeRising = tau2RangeRising;
algorithmInfo.tau1InitFalling = tau1InitFalling;
algorithmInfo.tau1RangeFalling = tau1RangeFalling;
algorithmInfo.tau2InitFalling = tau2InitFalling;
algorithmInfo.tau2RangeFalling = tau2RangeFalling;
algorithmInfo.tau1InitCombined = tau1InitCombined;
algorithmInfo.tau1RangeCombined = tau1RangeCombined;
algorithmInfo.tau2InitCombined = tau2InitCombined;
algorithmInfo.tau2RangeCombined = tau2RangeCombined;
algorithmInfo.Cm = Cm;
algorithmInfo.Ra = Ra;
algorithmInfo.Rs = Rs;
algorithmInfo.pulseResponseWindow = pulseResponseWindow;

% Store other measured outputs in allResults
allResults.pulseAmplitude = pulseAmplitude;
allResults.pulseWidth = pulseWidth;
allResults.holdPotential = holdPotential;
allResults.maxVoltage = maxVoltage;
allResults.voltageChange = voltageChange;
allResults.toUse = toUse;

allResults.tAvgRising = tAvgRising;
allResults.vAvgRising = vAvgRising;
allResults.paramsAvgRising = paramsAvgRising;
allResults.fitResultsAvgRising = fitResultsAvgRising;
allResults.cfitAvgRising = cfitAvgRising;
allResults.gofAvgRising = gofAvgRising;
allResults.algInfoAvgRising = algInfoAvgRising;

allResults.tAllRising = tAllRising;
allResults.vAllRising = vAllRising;
allResults.paramsAllRising = paramsAllRising;
allResults.fitResultsAllRising = fitResultsAllRising;
allResults.cfitAllRising = cfitAllRising;
allResults.gofAllRising = gofAllRising;
allResults.algInfoAllRising = algInfoAllRising;

allResults.tAvgFalling = tAvgFalling;
allResults.vAvgFalling = vAvgFalling;
allResults.paramsAvgFalling = paramsAvgFalling;
allResults.fitResultsAvgFalling = fitResultsAvgFalling;
allResults.cfitAvgFalling = cfitAvgFalling;
allResults.gofAvgFalling = gofAvgFalling;
allResults.algInfoAvgFalling = algInfoAvgFalling;

allResults.tAllFalling = tAllFalling;
allResults.vAllFalling = vAllFalling;
allResults.paramsAllFalling = paramsAllFalling;
allResults.fitResultsAllFalling = fitResultsAllFalling;
allResults.cfitAllFalling = cfitAllFalling;
allResults.gofAllFalling = gofAllFalling;
allResults.algInfoAllFalling = algInfoAllFalling;

allResults.tAllCombined = tAllCombined;
allResults.vAllCombined = vAllCombined;
allResults.paramsAllCombined = paramsAllCombined;
allResults.fitResultsAllCombined = fitResultsAllCombined;
allResults.cfitAllCombined = cfitAllCombined;
allResults.gofAllCombined = gofAllCombined;
allResults.algInfoAllCombined = algInfoAllCombined;

allResults.tAvgCombined = tAvgCombined;
allResults.vAvgCombined = vAvgCombined;
allResults.paramsAvgCombined = paramsAvgCombined;
allResults.fitResultsAvgCombined = fitResultsAvgCombined;
allResults.cfitAvgCombined = cfitAvgCombined;
allResults.gofAvgCombined = gofAvgCombined;
allResults.algInfoAvgCombined = algInfoAvgCombined;

%{
allResults.rmseRising = rmseRising;
allResults.rmseFalling = rmseFalling;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rmseAll = measure_error(xvecs, yvecs, phaseName, ...
                                    amplitudeEstimate, maxScalingFactor, ...
                                    tau1Init, tau2Init, ...
                                    tau1Range, tau2Range, pulseWidth)
%% Returns root-mean-squared error for each trace fitted individually
% TODO TODO TODO

% Setup fitting type, initial conditions and boundaries
[equationForm, aFittype, coeffInit, coeffLower, coeffUpper] = ...
    fit_setup_2exp(phaseName, amplitudeEstimate, maxScalingFactor, ...
                   tau1Init, tau2Init, tau1Range, tau2Range);

% For each trace, fit to curve and compute sum-of-squares error
nvecs = numel(xvecs);
fitObject_swp = cell(1, nvecs);        % curve fits for each individual sweep
gof_swp = cell(1, nvecs);        % goodness-of-fit statistics for each individual sweep
rmseAll = zeros(1, nvecs);        % root-mean-squared error (mV) for each individual sweep
parfor k = 1:nvecs
    % Fit data to equation
    if ~isempty(xvecs{k}) && ~isempty(yvecs{k})
        [cfit_swp{k}, gof_swp{k}] = fit(xvecs{k}, yvecs{k}, aFittype, 'StartPoint', coeffInit, ...
                        'Lower', coeffLower, 'Upper', coeffUpper); 
        rmseAll(k) = gof_swp{k}.rmse;
    else
        rmseAll(k) = Inf;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_geometry_and_passive_params (params1, params2)
%% Plots geometry of model cell from both the combined phase parameters and the falling phase parameters

% Plot model cell #1
plot_ball_stick('GeomParams', params1, 'EdgeColor', 'r', 'FaceColor', 'none');

% Plot model cell #2 if exists
if ~isempty(params2)
    plot_ball_stick('GeomParams', params2, 'EdgeColor', 'm', ...
                        'FaceColor', 'none');
end

% Add axes labels
xlabel('Position (um)')
ylabel('Position (um)')

% Adjust axes
xlim([-50, 50]);
ylim([-50, 50]);
% xlim([-100, 100]);
% ylim([-100, 100]);

% Get axes limits
ax = gca;
xLimits = get(ax, 'Xlim');
yLimits = get(ax, 'Ylim');

% Get axes ranges
xRange = xLimits(2) - xLimits(1);
yRange = yLimits(2) - yLimits(1);

% Get starting x and y positions
xpos = xLimits(1) + (1/30) * xRange;
ypos = yLimits(1) + (29/30) * yRange;

% Show parameter values
text(xpos, ypos, ['L = ', num2str(params1.L)], ...
    'FontSize', 8, 'Color', 'r');
text(xpos, ypos - (1/15) * yRange, ...
    ['rho_1 = ', num2str(params1.rho)], ...
    'FontSize', 8, 'Color', 'r');
text(xpos, ypos - (2/15) * yRange, ...
    ['Rsoma_1 = ', num2str(params1.Rsoma, '%.0f'), ' MOhm'], ...
    'FontSize', 8, 'Color', 'r');
text(xpos, ypos - (3/15) * yRange, ...
    ['Rdend_1 = ', num2str(params1.Rdend, '%.0f'), ' MOhm'], ...
    'FontSize', 8, 'Color', 'r');
text(xpos, ypos - (4/15) * yRange, ...
    ['taum_1 = ', num2str(params1.taum), ' ms'], ...
    'FontSize', 8, 'Color', 'r');
text(xpos, ypos - (5/15) * yRange, ...
    ['Rm_1 = ', num2str(params1.Rm), ' Ohm-cm^2'], ...
    'FontSize', 8, 'Color', 'r');

text(xpos, ypos - (9/15) * yRange, ...
    ['L_2 = ', num2str(params2.L)], ...
    'FontSize', 8, 'Color', 'm');
text(xpos, ypos - (10/15) * yRange, ...
    ['rho_2 = ', num2str(params2.rho)], ...
    'FontSize', 8, 'Color', 'm');
text(xpos, ypos - (11/15) * yRange, ...
    ['Rsoma_2 = ', num2str(params2.Rsoma, '%.0f'), ' MOhm'], ...
    'FontSize', 8, 'Color', 'm');
text(xpos, ypos - (12/15) * yRange, ...
    ['Rdend_2 = ', num2str(params2.Rdend, '%.0f'), ' MOhm'], ...
    'FontSize', 8, 'Color', 'm');
text(xpos, ypos - (13/15) * yRange, ...
    ['taum_2 = ', num2str(params2.taum), ' ms'], ...
    'FontSize', 8, 'Color', 'm');
text(xpos, ypos - (14/15) * yRange, ...
    ['Rm_2 = ', num2str(params2.Rm), ' Ohm-cm^2'], ...
    'FontSize', 8, 'Color', 'm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% OLD CODE

ft_R = fittype('a*(exp(-x/b)-1)+c*(exp(-x/d)-1)');    
ft_F = fittype('-a*exp(-x/b)-c*exp(-x/d)');        

dvss_S = paramsRising.responseAmplitude;
Rinput_S = paramsRising.Rinput;
if strcmp(phaseName, 'falling')
    dvss_L = paramsFalling.responseAmplitude;
    Rinput_L = paramsFalling.Rinput;
end

    dvss_min = min(dvss_L, dvss_S);        % the more negative of the two steady-state voltage changes
    ylim([(dvss_min - 0.5) 0.5]);        % Add 0.5 mV above and below

line(xLimits, [dvss_S dvss_S], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.5);
if strcmp(phaseName, 'falling')
    line(xLimits, [dvss_L dvss_L], 'Color', 'm', 'LineStyle', '--', 'LineWidth', 0.5);
end

if nargin < 10 || (nargin >= 10 && isempty(dataMode))
    dataMode = 0;
end

fitObject = cell(1, nSweeps);                % double exponential fit

        indBaseline(iSwp, :) = cpstart - round(0.5/siMs) - fliplr(indMeanVoltageWindow);    % base indices

aFittype = fittype('a*exp(-x/b)+c*exp(-x/d)+e');    % double exponential
% coeff = coeffnames(aFittype);            % {'a'; 'b'; 'c'; 'd'; 'e'}
    cfit_R1 = fit(tAllRising, vAllRising, aFittype, ...
        'StartPoint', [meanVoltageChange, tau1Init, tau2Init, meanVoltageChange, tau1Init, tau2Init, -2 * meanVoltageChange], ...
        'Lower', [0, 0, 0, 0, -20 * meanVoltageChange], ...
        'Upper', [10 * meanVoltageChange, 100, 10 * meanVoltageChange, 100, 0]); 

h = figure('Visible', 'off');

        plot(tvec0(indCpr), ivec0s(indCpr, iSwp), 'b'); hold on; 

ypos = responseAmplitude + 1;

    if voltageChange(iSwp) > 0 && pulseAmplitude(iSwp) < 0    % exclude faulty traces, but this is really unnecessary
    end

% L = double(solve(abs(C1/(2*C0*tau2/tau1-C1)) == cot(alpha1*x)*(cot(alpha1*x)-1/(alpha1*x)), x));    

eqnS = [num2str(CS1, 2), '*(exp(-x/', num2str(tau1, 2), ...
    ')-1)+', num2str(CS2, 2), '*(exp(-x/', num2str(tau2, 2), ')-1)'];
eqnL = [num2str(CL1, 2), '*(exp(-x/', num2str(tau1, 2), ...
    ')-1)+', num2str(CL2, 2), '*(exp(-x/', num2str(tau2, 2), ')-1)'];

typtau = 2;        % typical tau ~ 2 ms
tau_max = 500;        % tau cannot exceed 500 ms

        coeffInit(iCoeff) = typtau;

% Set up initial conditions and boundary conditions for fit()
for iCoeff = 1:nCoeffs
    coeffLower(iCoeff) = 0;
    if coeffNames{iCoeff} == 'a' || coeffNames{iCoeff} == 'c'
        coeffInit(iCoeff) = amplitudeEstimate;
        coeffUpper(iCoeff) = maxScalingFactor * amplitudeEstimate;
    elseif coeffNames{iCoeff} == 'b' || coeffNames{iCoeff} == 'd'
        coeffInit(iCoeff) = mean(tau_min, tau_max);
        coeffLower(iCoeff) = tau_min;
        coeffUpper(iCoeff) = tau_max;
    end
end

    % This won't work
        coeffInit(iCoeff) = mean(tau1Range(1), tau1Range(2));
        coeffInit(iCoeff) = mean(tau2Range(1), tau2Range(2));

        coeffInit(iCoeff) = tau1Range(1);
        coeffInit(iCoeff) = tau2Range(2);

%% Check arguments
if nargin < 2
    error('Not enough input arguments, type ''help find_passive_params'' for usage');
elseif isempty(tvec0) || isempty(ivec0s) || isempty(vvec0s)
    error('First three inputs cannot be empty!');
elseif ~isnumeric(tvec0) || ~isnumeric(ivec0s) || ~isnumeric(vvec0s)
    error('First three inputs must be numbers or numeric arrays!');
elseif ~isequal(size(tvec0, 1), size(ivec0s, 1))
    error('Time and current vectors do not have the same length!');
elseif ~isequal(size(tvec0, 1), size(vvec0s, 1))
    error('Time and voltage vectors do not have the same length!');
elseif size(tvec0, 2) > 1
    error('tvec0 must be a column vector!');
elseif nargin >= 4 && length(pulseWindow) < 2
    error('pulseWindow must have a start and an end!');
elseif nargin >= 5 && length(pulseResponseWindow) < 2
    error('pulseResponseWindow must have a start and an end!');
end

elseif nargin >= 7 && (plotFlag ~= 1 && plotFlag ~= 0 && plotFlag ~= false && plotFlag ~= true)
    error('plotFlag is out of range!');
elseif nargin >= 8 && ~ischar(outFolder)
    error('outFolder must be a character array!');
elseif nargin >= 9 && ~ischar(fileBase)
    error('fileBase must be a character array!');
elseif nargin >= 10 && ~isnumeric(ivec1s)
    error('ivec1s must be a numeric array!');

if nargin < 7
    plotFlag = 0;
end
if nargin < 8
    outFolder = pwd;
end
if nargin < 9
    fileBase = 'unnamed';
end

if nargin >= 11
else
    suffix = '';
    titleMod = '';
end

directories = {'/cprRising/', '/cprFalling/', '/passive/'};
for k = 1:numel(directories)
    temp = directories{k}(2:end);        % remove the leading '/'
    directories{k} = [ '/', strrep(temp, '/', [suffix, '/']) ];
end

%% Only do the following if dataMode is used

function [passiveParams, pulseAmplitude, pulseWidth, voltageChange, rmseRising, rmseFalling, params_L_F2, params_S_R2, params_L_F1, paramsAllRising] = find_passive_params (tvec0, ivec0s, vvec0s, pulseWindow, pulseResponseWindow, pulseMidpoint, plotFlag, outFolder, fileBase, ivec1s, dataMode)

% Median filter current vectors
ivec1s = zeros(nSamples, nSweeps);
parfor iSwp = 1:nSweeps
    ivec1s(:, iSwp) = medfilt1(ivec0s(:, iSwp), medianFilterWindowSamples);
end

indCpStartWin = idxStart:idxMid;            % indices for finding current pulse start
indCpEndWin = idxMid:idxEnd;                % indices for finding current pulse end
idxMid = find(tvec0 <= pulseMidpoint, 1, 'last');

idxStart = find(tvec0 >= pulseResponseWindow(1), 1);
idxEnd = find(tvec0 <= pulseResponseWindow(2), 1, 'last');

idxCpStart = zeros(nSweeps, 1);
idxCpEnd = zeros(nSweeps, 1);
parfor iSwp = 1:nSweeps
    [idxCpStart(iSwp), idxCpEnd(iSwp)] = ...
        find_pulse_endpoints(ivecCpr(:, iSwp));
end

%% TODO Move
% Calculate the mean voltage calculation window in samples
meanVoltageWindow = 0.5;    % width in ms for calculating mean voltage 
                            %   for input resistance calculations
meanVoltageWindowSamples = round(meanVoltageWindow/siMs);

% Find the current pulse amplitudes
indCpStartWin = idxStart:idxMid;            % indices for finding current pulse start
indCpEndWin = idxMid:idxEnd;                % indices for finding current pulse end
idxMid = find(tvec0 <= pulseMidpoint, 1, 'last');
indBaseline = zeros(nSweeps, length(indMeanVoltageWindow));
indSteady = zeros(nSweeps, length(indMeanVoltageWindow));
ivec1sPart1 = ivec1s(indCpStartWin, :);            % Use median-filtered current trace
ivec1sPart2 = ivec1s(indCpEndWin, :);
ivec1sPart1Begin = indCpStartWin(1);
ivec1sPart2Begin = indCpEndWin(1);

% Compute the indices for taking the mean of voltages
indMeanVoltageWindow = 1:meanVoltageWindowSamples;

ivec1sPart3 = ivec1s(idxStart - 1 + indMeanVoltageWindow, :);     % indices for measuring cp baseline
ivec1sPart4 = ivec1s(idxMid - 1 + indMeanVoltageWindow, :);       % indices for measuring cp peak

%% Find current pulses for each sweep
firstDipPoint = cell(1, nSweeps);   % index of current pulse start
pulseAmplitude = zeros(1, nSweeps);    % current pulse amplitude (pA)
pulseWidth = zeros(1, nSweeps);        % pulse width (ms)
voltageChange = zeros(1, nSweeps);  % overall change in membrane potential recorded (mV)
tvecsRising = cell(1, nSweeps);          % adjusted rising phase time vectors (start at 0)
vvecsRising = cell(1, nSweeps);          % adjusted rising phase voltage vectors (start at around 0)
tvecsFalling = cell(1, nSweeps);          % adjusted falling phase time vectors (start at 0)
vvecsFalling = cell(1, nSweeps);          % adjusted falling phase voltage vectors (asymptote to around 0)
used = zeros(1, nSweeps);           % whether the sweep is used for fitting
% for iSwp = 1:nSweeps              % FOR each sweep
parfor iSwp = 1:nSweeps             % FOR each sweep
    pulseAmplitude(iSwp) = mean(ivec1sPart4(:, iSwp)) - mean(ivec1sPart3(:, iSwp));    % current pulse amplitude (pA)
    firstDipPoint{iSwp} = find(ivec1sPart1(:, iSwp) > pulseAmplitude(iSwp)/4, 1, 'last');
    if isempty(firstDipPoint{iSwp})                % exclude faulty traces
        tvecsRising{iSwp} = [];
        vvecsRising{iSwp} = [];
        tvecsFalling{iSwp} = [];
        vvecsFalling{iSwp} = [];
        used(iSwp) = false;
    else
        % Find the voltage trace corresponding to the current pulse response
        cpstart = (ivec1sPart1Begin - 1) + firstDipPoint{iSwp};    % index of current pulse start
        before_rise_pt = find(ivec1sPart2(:, iSwp) < pulseAmplitude(iSwp) * 3/4, 1, 'last');
        cpend = (ivec1sPart2Begin - 1) + before_rise_pt;    % index of current pulse end
        pulseWidth(iSwp) = (cpend - cpstart) * siMs;             % pulse width (ms)
        indBaseline(iSwp, :) = cpstart - fliplr(indMeanVoltageWindow);        % base indices
        indSteady(iSwp, :) = cpend - fliplr(indMeanVoltageWindow);        % last indices of the current pulse
        basev = mean(vvec0s(indBaseline(iSwp, :), iSwp));
        lastv = mean(vvec0s(indSteady(iSwp, :), iSwp));
        voltageChange(iSwp) = basev - lastv;         % change in membrane potential on voltage trace (mV)
        vmax = max(vvec0s(:, iSwp));        % maximum voltage in current pulse response window (mV)
        if voltageChange(iSwp) <= 0 || pulseAmplitude(iSwp) >= 0 || pulseWidth(iSwp) <= 0 ...
            || vmax > spikeThresholdInit        % exclude faulty traces, including those with spontaneous spikes
            tvecsRising{iSwp} = [];
            vvecsRising{iSwp} = [];
            tvecsFalling{iSwp} = [];
            vvecsFalling{iSwp} = [];
            used(iSwp) = false;
        else
            tvecsRising{iSwp} = tvec0(cpstart:cpend) - tvec0(cpstart);
            vvecsRising{iSwp} = vvec0s(cpstart:cpend, iSwp) - basev;
            tvecsFalling{iSwp} = tvec0(cpend:idxEnd) - tvec0(cpend);
            vvecsFalling{iSwp} = vvec0s(cpend:idxEnd, iSwp) - basev;
            used(iSwp) = true;
        end
    end
end

if ~isempty(find(voltageChange > 0, 1))
    meanVoltageChange = mean(voltageChange(voltageChange > 0));        % average recorded voltage change (mV), excluding faulty traces
else
    meanVoltageChange = NaN;
end
if ~isempty(find(pulseAmplitude < 0, 1))
    meanPulseAmplitude = mean(pulseAmplitude(pulseAmplitude < 0));            % average current pulse amplitude (pA), excluding faulty traces
else
    meanPulseAmplitude = NaN;
end
if ~isempty(find(pulseWidth > 0, 1))
    meanPulseWidth = mean(pulseWidth(pulseWidth > 0));            % average pulse width (ms), excluding faulty traces
else
    meanPulseWidth = NaN;
end

coeffLower(iCoeff) = 0;
if coeffNames{iCoeff} == 'a' || coeffNames{iCoeff} == 'c'
    coeffInit(iCoeff) = amplitudeEstimate;
    coeffUpper(iCoeff) = maxScalingFactor * amplitudeEstimate;
elseif coeffNames{iCoeff} == 'b'
    coeffInit(iCoeff) = tau1Init;
    coeffLower(iCoeff) = tau1Range(1);
    coeffUpper(iCoeff) = tau1Range(2);
elseif coeffNames{iCoeff} == 'd'
    coeffInit(iCoeff) = tau2Init;
    coeffLower(iCoeff) = tau2Range(1);
    coeffUpper(iCoeff) = tau2Range(2);
end

eqFormRising = '-a*(1-exp(-x/b))+-c*(1-exp(-x/d))';
                            % double exponential equation form for rising phase
eqFormFalling = '-a*exp(-x/b)-c*exp(-x/d)';
                            % double exponential equation form for falling phase

function [fitObject, goodnessOfFit, algorithmInfo] = ...
                fit_2exp(xvec, yvec, phaseName, ...
                            amplitudeEstimate, maxScalingFactor, ...
                            tau1Init, tau2Init, tau1Range, tau2Range)
%% Fits a rising or falling double exponential curve to data

% Setup fitting type, initial conditions and boundaries
[equationForm, aFittype, coeffInit, coeffLower, coeffUpper] = ...
    fit_setup_2exp(phaseName, ...
                    'AmplitudeEstimate', amplitudeEstimate, ...
                    'MaxScalingFactor', maxScalingFactor, ...
                    'Tau1Init', tau1Init, ...
                    'Tau2Init', tau2Init, ...
                    'Tau1Range', tau1Range, ...
                    'Tau2Range', tau2Range);

% Fit data
[fitObject, goodnessOfFit, algorithmInfo] = ...
    fit(xvec, yvec, aFittype, 'StartPoint', coeffInit, ...
        'Lower', coeffLower, 'Upper', coeffUpper); 

% Store the equation form
algorithmInfo.equationForm = equationForm;

tAllRising = [];
vAllRising = [];
tAllFalling = [];
vAllFalling = [];
for iSwp = 1:nSweeps
    tAllRising = cat(1, tAllRising, tvecsRising{iSwp});
    vAllRising = cat(1, vAllRising, vvecsRising{iSwp});
    tAllFalling = cat(1, tAllFalling, tvecsFalling{iSwp});
    vAllFalling = cat(1, vAllFalling, vvecsFalling{iSwp});
end

tvecsRisingToUse = tvecsRising(~cellfun(@isempty, tvecsRising));        % Remove all elements that are empty
vvecsRisingToUse = vvecsRising(~cellfun(@isempty, vvecsRising));        % ditto
tvecsFallingToUse = tvecsFalling(~cellfun(@isempty, tvecsFalling));        % ditto
vvecsFallingToUse = vvecsFalling(~cellfun(@isempty, vvecsFalling));        % ditto

for iSwp = 1:numel(vvecsRisingToUse)
    vvecsRisingToUse{iSwp} = vvecsRisingToUse{iSwp}(1:minLengthRising);
end

vAvgRising = mean(cell2mat(vvecsRisingToUseAligned), 2);

% Compute the minimum length of all rising phase traces
minLengthRising = min(cellfun(@length, tvecsRisingToUse));

% Compute an averaged rising phase trace
if minLengthRising > 0
    % Adjust the first time vector to this minimum length
    tAvgRising = tvecsRisingToUse{1}(1:minLengthRising);

    % Adjust voltage vectors to this minimum length
    vvecsRisingToUseAligned = ...
        cellfun(@(x) x(1:minLengthRising), vvecsRisingToUse, ...
                'UniformOutput', false);

    % Compute the mean of all adjusted voltage vectors
    vAvgRising = mean(horzcat(vvecsRisingToUseAligned{:}), 2);
else
    tAvgRising = [];
    vAvgRising = [];
end

% Compute the minimum length of all falling phase traces
minLengthFalling = min(cellfun(@length, tvecsFallingToUse));

% Compute an averaged falling phase trace
if minLengthFalling > 0
    tAvgFalling = tvecsFallingToUse{1}(1:minLengthFalling);            % ditto for falling phase
    for iSwp = 1:numel(vvecsRisingToUse)
        vvecsFallingToUse{iSwp} = vvecsFallingToUse{iSwp}(1:minLengthFalling);    % ditto for falling phase
    end
    vAvgFalling = mean(cell2mat(vvecsFallingToUse), 2);        % ditto for falling phase
else
    tAvgFalling = [];
    vAvgFalling = [];
end

%% Fit a double exponential to each voltage trace
if ~isempty(tAllRising)
    % Fit pooled rising phase data with double exponential
    [cfit_R1, eqnS_R1, eqnL_R1, CS0_R1, tau1_R1, CS1_R1, tau2_R1, ~, ~] = ...
        fit_pulse_response (tAllRising, vAllRising, 'rising', meanVoltageChange, maxScalingFactor, tau1InitRising, tau2InitRising, tau1RangeRising, tau2RangeRising, meanPulseWidth);

    % Estimate parameters from short-pulse coefficients
    [paramsAllRising] = coeff2params (CS0_R1, tau1_R1, CS1_R1, tau2_R1, meanPulseAmplitude, Cm, Ra, Rs);
    if printFieldsFlag
        print_structure(paramsAllRising);
    end
else
    paramsAllRising = NaN;
end
if ~isempty(tAllFalling)
    % Fit pooled falling phase data with double exponential
    [cfit_F1, eqnS_F1, eqnL_F1, ~, tau1_F1, ~, tau2_F1, CL0_F1, CL1_F1] = ...
        fit_pulse_response (tAllFalling, vAllFalling, 'falling', meanVoltageChange, maxScalingFactor, tau1InitFalling, tau2InitFalling, tau1RangeFalling, tau2RangeFalling, meanPulseWidth);

    % Estimate parameters from long-pulse coefficients
    [params_L_F1] = coeff2params (CL0_F1, tau1_F1, CL1_F1, tau2_F1, meanPulseAmplitude, Cm, Ra, Rs);
    if printFieldsFlag
        print_structure(params_L_F1);
    end
else
    params_L_F1 = NaN;
end
if ~isempty(tAvgRising)
    % Fit averaged rising phase data with double exponential
    [cfit_R2, eqnS_R2, eqnL_R2, CS0_R2, tau1_R2, CS1_R2, tau2_R2, ~, ~] = ...
        fit_pulse_response (tAvgRising, vAvgRising, 'rising', meanVoltageChange, maxScalingFactor, tau1InitRising, tau2InitRising, tau1RangeRising, tau2RangeRising, meanPulseWidth);

    % Estimate parameters from short-pulse coefficients
    [params_S_R2] = coeff2params (CS0_R2, tau1_R2, CS1_R2, tau2_R2, meanPulseAmplitude, Cm, Ra, Rs);
    if printFieldsFlag
        print_structure(params_S_R2);
    end
else
    params_S_R2 = NaN;
end
if ~isempty(tAllFalling)
    % Fit averaged falling phase data with double exponential
    [cfit_F2, eqnS_F2, eqnL_F2, ~, tau1_F2, ~, tau2_F2, CL0_F2, CL1_F2] = ...
        fit_pulse_response (tAllFalling, vAllFalling, 'falling', meanVoltageChange, maxScalingFactor, tau1InitFalling, tau2InitFalling, tau1RangeFalling, tau2RangeFalling, meanPulseWidth);

    % Estimate parameters from long-pulse coefficients
    [params_L_F2] = coeff2params (CL0_F2, tau1_F2, CL1_F2, tau2_F2, meanPulseAmplitude, Cm, Ra, Rs);
    if printFieldsFlag
        print_structure(params_L_F2);
    end
else
    params_L_F2 = NaN;
end

if exist('params_L_F2', 'var') == 1 && isstruct(params_L_F2)
    passiveParams = params_L_F2;
elseif exist('params_S_R2', 'var') == 1 && isstruct(params_S_R2)
    passiveParams = params_S_R2;
elseif exist('params_L_F1', 'var') == 1 && isstruct(params_L_F1)
    passiveParams = params_L_F1;
elseif exist('paramsAllRising', 'var') == 1 && isstruct(paramsAllRising)
    passiveParams = paramsAllRising;
else
    passiveParams = NaN;
end   

if nToUse > 0
else
    meanPulseWidth = NaN;
    meanPulseAmplitude = NaN;
    meanVoltageChange = NaN;
end

[pulseParams, pulseData] = parse_pulse(ivecsCpr);

% Extract the indices of the endpoints and midpoint of the current pulses
idxCpStart = pulseParams.idxBeforeStart;
idxCpEnd = pulseParams.idxBeforeEnd;
idxCpMid = pulseParams.idxMidpoint;

% Get the number of samples
nSamples = length(tvec0);

%                   - 'PulseMidpoint': approximate midpoint of the 
%                                       current pulse (ms)
%                   must be within range of pulseResponseWindow
%                   default == mean(pulseWindow)
pulseMidpointDefault = [];          % set later
addParameter(iP, 'PulseMidpoint', pulseMidpointDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
pulseMidpoint = iP.Results.PulseMidpoint;
elseif ~isempty(pulseMidpoint) && ...
        (pulseMidpoint < pulseResponseWindow(1) || ...
            pulseMidpoint > pulseResponseWindow(2))
    error('Pulse Midpoint is out of range!');
% TODO: Use this info to help detect current pulse?
if isempty(pulseMidpoint)
    pulseMidpoint = mean(pulseWindow);
end

toUse = pulseWidth >= 0 & pulseAmplitude <= 0 & ...
        voltageChange <= 0 & maxVoltage > spikeThresholdInit;

% Extract the separated response indices and the numbers of samples in them
nSamplesRising = responseParams.nSamplesRising;
nSamplesFalling = responseParams.nSamplesFalling;
indRising = responseData.indRising;
indFalling = responseData.indFalling;

% Extract the vectors already transformed to cell arrays
vvecsCprCell = responseData.vectors;

% Extract the mean baseline values
meanBaseValue = responseParams.meanBaseValue;

% Generate shifted rising/falling phase vectors so that time starts at zero
%   and steady state value is zero
% Note: This will make curve fitting easier
tvecsRising = arrayfun(@(x) transpose(1:x) * siMs, nSamplesRising, ...
                        'UniformOutput', false);
vvecsRising = cellfun(@(x, y, z) x(y) - z, ...
                        vvecsCprCell, indRising, {meanBaseValue}, ...
                        'UniformOutput', false);
tvecsFalling = arrayfun(@(x) transpose(1:x) * siMs, nSamplesFalling, ...
                        'UniformOutput', false);
vvecsFalling = cellfun(@(x, y, z) x(y) - z, ...
                        vvecsCprCell, indFalling, {meanBaseValue}, ...
                        'UniformOutput', false);

plot_pulse(nSweeps, toUse, indCpr, tvec0, ivec1s, firstDipPoint, indBaseline, indSteady, pulseWindow, meanPulseAmplitude);
plot_pulse_response(nSweeps, toUse, indCpr, tvec0, vvec0s, firstDipPoint, indBaseline, indSteady, pulseWindow)
plot_pulse(nSweeps, toUse, indCpr, tvec0, ivec1s, firstDipPoint, indBaseline, indSteady, pulseResponseWindow, meanPulseAmplitude);
plot_pulse_response(nSweeps, toUse, indCpr, tvec0, vvec0s, firstDipPoint, indBaseline, indSteady, pulseResponseWindow)

%function plot_pulse (nSweeps, toUse, indCpr, tvec0, ivec1s, firstDipPoint, indBaseline, indSteady, xLimits, meanPulseAmplitude)

function plot_pulse_response (nSweeps, toUse, indCpr, tvec0, vvec0s, firstDipPoint, indBaseline, indSteady, xLimits)
%% Plots voltage responses to current pulses

for iSwp = 1:nSweeps            % Plot each voltage trace
    if toUse(iSwp)
        plot(tvec0(indCpr), vvec0s(indCpr, iSwp), '-'); hold on; 
    else                % if not used, plot as dotted line
        plot(tvec0(indCpr), vvec0s(indCpr, iSwp), '--'); hold on; 
    end
    if ~isempty(firstDipPoint{iSwp})
        plot(tvec0(indBaseline(iSwp, 1)), vvec0s(indBaseline(iSwp, 1), iSwp), '>');
        plot(tvec0(indBaseline(iSwp, end)), vvec0s(indBaseline(iSwp, end), iSwp), '<');
        plot(tvec0(indSteady(iSwp, 1)), vvec0s(indSteady(iSwp, 1), iSwp), '>');
        plot(tvec0(indSteady(iSwp, end)), vvec0s(indSteady(iSwp, end), iSwp), '<');
    end
end
xlabel('Time (ms)')
ylabel('Voltage (mV)')
xlim(xLimits);

if ~isempty(cfitAllRising)
end
if ~isempty(tAvgRising)
end
if ~isempty(tAllFalling)
end
if ~isempty(tAvgFalling)
end
if ~isempty(tAllRising)
end
if ~isempty(tAvgRising)
end
if ~isempty(tAllFalling)
end
if ~isempty(tAvgFalling)
end
if ~isempty(tAllRising) && ~isempty(tAllFalling)
end
if ~isempty(tAvgRising) && ~isempty(tAvgFalling)
end

subplot(2,2,3);
    plot_cfit_pulse_response('falling', pulseResponseWindow, cfit_F1, tAllFalling, vAllFalling, params_L_F1, ...
            eqnS_F1, eqnL_F1, tau1_F1, tau2_F1, CL0_F1, CL1_F1, meanPulseWidth);
subplot(2,2,4);
    plot_cfit_pulse_response('falling', pulseResponseWindow, cfit_F2, tAvgFalling, vAvgFalling, params_L_F2, ...
            eqnS_F2, eqnL_F2, tau1_F2, tau2_F2, CL0_F2, CL1_F2, meanPulseWidth);

medianFilterWindowSamples = round(medFiltWindowMs/siMs);

subplot(3,2,1);
plot_cfit_pulse_response(tAllRising, vAllRising, ...
                        'FitObject', cfitAllRising, ...
                        'FitResults', fitResultsAllRising, ...
                        'PassiveParams', paramsAllRising, ...
                        'LegendLocation', 'suppress');
text(0.1, 1.1, sprintf('Total number of sweeps: %d', nSweeps), ...
    'Units', 'normalized');
subplot(3,2,2);
plot_cfit_pulse_response(tAvgRising, vAvgRising, ...
                        'FitObject', cfitAvgRising, ...
                        'FitResults', fitResultsAvgRising, ...
                        'PassiveParams', paramsAvgRising, ...
                        'LegendLocation', 'suppress');
subplot(3,2,3);
plot_cfit_pulse_response(tAllFalling, vAllFalling, ...
                        'FitObject', cfitAllFalling, ...
                        'FitResults', fitResultsAllFalling, ...
                        'PassiveParams', paramsAllFalling, ...
                        'LegendLocation', 'suppress');
subplot(3,2,4);
plot_cfit_pulse_response(tAvgFalling, vAvgFalling, ...
                        'FitObject', cfitAvgFalling, ...
                        'FitResults', fitResultsAvgFalling, ...
                        'PassiveParams', paramsAvgFalling, ...
                        'LegendLocation', 'suppress');
subplot(3,2,5);
plot_geometry_and_passive_params(paramsAllRising, paramsAllFalling);
subplot(3,2,6);
plot_geometry_and_passive_params(paramsAvgRising, paramsAvgFalling);

%                   - 'DataMode': data mode
%                   must be one of the following
%                       -1 - none
%                       0 - all data
%                       1 - all of g incr = 100%, 200%, 400%
%                       2 - all of g incr = 100%, 200%, 400% 
%                               but exclude cell-pharm-g_incr sets 
%                               containing problematic sweeps
%                   default == -1

%       cd/m3ha_specs_for_datamode.m (if dataMode is used)

addParameter(iP, 'DataMode', dataModeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
dataModeDefault = 0;
dataMode = iP.Results.DataMode;
elseif dataMode ~= 0 && dataMode ~= 1 && dataMode ~= 2
    error('dataMode out of range!');

% Set suffix and title modification according to dataMode
[suffix, titleMod] = m3ha_specs_for_datamode(dataMode);

[idxStart, idxEnd] = find_window_endpoints(pulseResponseWindow, tvec0);
indCpr = transpose(idxStart:idxEnd);

indCpr = transpose(endPoints(1):endPoints(2));

% Restrict the vectors to the current pulse response window
tvecCpr = tvec0(indCpr);
vvecsCpr = vvec0s(indCpr, :);
ivecsCpr = ivec1s(indCpr, :);

%} 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
