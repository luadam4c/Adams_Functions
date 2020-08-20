function [eventInfo, eventClass, dataLowpass, noiseLevel] = ...
                minEASE_detect_gapfree_events(data, siMs, directionPsc, varargin)
%% Spontaneous EPSC/IPSC detector and classifyer
% Usage: [eventInfo, eventClass, dataLowpass, noiseLevel] = ...
%               minEASE_detect_gapfree_events(data, siMs, directionPsc, varargin)
% Explanation: 
%   (1) Lowpass filter data
%   (2) Compute a RMS noise of the "Gaussian parts" of the lowpass-filtered data
%   (3) Find all directional events in the data after applying a 
%           moving-average filter, but compute event statistics on the 
%           lowpass-filtered data
%   (4) Classify events and identify PSCs
%       "full decay time" and "half decay time" would ignore removed events
%   (5) Compute an averaged Type I PSC trace of the original data
%
% Outputs:
%       eventInfo   - TODO: Description of eventInfo
%                   specified as a TODO
%       eventClass  - TODO: Description of eventInfo
%                   specified as a TODO
%       dataLowpass - TODO: Description of dataLowpass
%                   specified as a TODO
%
% Side Effects:
%       Plots TODO
%
% Arguments:    
%       data            - vector of current trace data
%                       must be a numeric vector
%       siMs            - sampling interval in ms/sample
%                       must be a numeric positive scalar
%       directionPsc    - direction of post-synaptic current to detect
%                       must be an unambiguous, case-insensitive match to one of: 
%                           'Excitatory' or 'E' - downward peaks (EPSCs)
%                           'Inhibitory' or 'I' - upward peaks (IPSCs)
%       varargin    - 'LowpassCutoff': 3 dB cutoff frequency (Hz) 
%                                       of lowpass filter
%                   must be a nonnegative scalar
%                   default == 3000
%                   - 'LowpassNpoles': order of Butterworth lowpass filter
%                                   (degree of transfer function denominator)
%                   must be a positive integer scalar
%                   default == 8
%                   - 'NoiseWindowSize': size in samples for a noise window
%                   must be a positive integer scalar
%                   default == minimum of 5 samples or 
%                               however many to make 1000 windows
%                   - 'ZSkewnessThres': skewness z-score threshold for deciding 
%                                       if a noise window is Gaussian
%                   must be a nonnegative scalar
%                   default == 0.2
%                   - 'ZExcessKurtosisThres': excess kurtosis z-score threshold
%                                   for deciding if a noise window is Gaussian
%                   must be a nonnegative scalar
%                   default == 0.2
%                   - 'Signal2Noise': Signal/noise ratio of PSCs to noise
%                   must be a positive scalar
%                   default == 2
%                   - 'MinEventThreshold': minimum amplitude threshold (pA)
%                                           for a directional event
%                   must be a positive scalar
%                   default == 8 pA
%                   - 'BaselineWindowMs': window for computing baseline (ms)
%                   must be a nonnegative scalar
%                   default == 5 ms
%                   - 'MaxBelowBasePerc': maximum below baseline percentage (%)
%                   must be a nonnegative scalar
%                   default == 100 %
%                   - 'SmoothWindowMs': moving average filter window (ms)
%                   must be a nonnegative scalar
%                   default == 0.5 ms
%                   - 'MinPscAmpAbs': minimum allowable PSC amplitude (pA)
%                   must be a nonnegative scalar
%                   default == 10 pA
%                   - 'MinPscAmpRel': relative minimum allowable PSC amplitude 
%                                       (% of maximum)
%                   must be a nonnegative scalar less than 100
%                   default == 10 %
%                   - 'MaxPscRiseMs': maximum allowable PSC rise time (ms)
%                   must be a nonnegative scalar
%                   default == 8 ms
%                   - 'MinPscDecayMs': minimum allowable PSC decay time (ms)
%                   must be a nonnegative scalar
%                   default == 0.3 ms
%                   - 'MaxPscDecayMs': maximum allowable PSC decay time (ms)
%                   must be a nonnegative scalar
%                   default == 50 ms
%                   - 'SealTestWindowMs': seal test window (ms)
%                   must be a nonnegative, nondecreasing vector of 2 elements
%                   default == [1000, 1050] ms
%                   - 'TraceLengthMs': total trace length for each PSC (ms)
%                   must be a numeric positive scalar
%                   default == 50
%                   - 'BeforePeakMs': trace length before break point 
%                                       for each PSC (ms)
%                   must be a numeric positive scalar
%                   default == 3
%                   - 'OutputDirectory': directory to place outputs
%                   must be a string scalar or a character vector
%                   default == '/media/shareX/share/minEASE/output/'
%                   - 'FigTypes': figure type(s) for saving; 
%                       e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by the saveas() function 
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - 'FileNo': file number (sweep number)
%                   must be an integer scalar
%                   default == 0
%                   - 'OutputLabel': Output label
%                   must be a string scalar or a character vector
%                   default == ['Swp', num2str(fileNo)]
%                   - 'PlotEventDetectionFlag': whether to plot event detection
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotAverageTraceFlag': whether to plot average trace
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'MessageMode' - how message boxes are shown
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'wait'  - stops program and waits for the user
%                                   to close the message box
%                       'show'  - does not stop program but still show the
%                                   message box
%                       'none'  - neither stop program nor show a message box
%                   default == 'wait'
%                   - 'SkipDetection': whether to skip event detection
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Verbose' - whether to print to standard output
%                                   regardless of message mode
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   
% Requires:
%       cd/minEASE_plot_gapfree_events.m
%       cd/minEASE_compute_plot_average_psc.m
%       cd/compute_rms_Gaussian.m
%       cd/find_directional_events.m
%       cd/identify_bursts.m
%       cd/isfigtype.m
%       cd/print_or_show_message.m
%       cd/save_all_figtypes.m
%
% Used by:
%       cd/minEASE.m

% File History:
% ---------- Created by Koji Takahashi & Mark P Beenhakker
% 2017-05-21 AL - Renamed gapfree_analysis_101115_trunc.m -> detect_events.m
% 2017-05-23 AL - Changed argument dt -> siMs in milliseconds
% 2017-05-24 AL - Now uses butter() instead of butterfilter()
% 2017-05-25 AL - Renamed detect_events.m -> minEASE_detect_gapfree_events.m
% 2017-05-25 AL - Added directionFactor so that both thresholds 
%                   can be plotted correctly
% 2017-05-30 AL - Simplified code for identifying PSCs
% 2017-06-01 AL - pscInfo now has same number of columns as eventInfo
% 2017-06-05 AL - Now uses filtfilt() instead of filter()
% 2017-06-05 AL - Removed burst detection and now detects Type II and 
%                   Type III events instead
% 2017-06-06 AL - Added SweepLabel as a parameter
% 2017-06-13 AL - Now returns noiseLevel as an output
% 2017-06-14 AL - Made plotEventDetectionFlag & plotAverageTraceFlag parameters
% 2017-07-24 AL - Added Average Type II & Type III PSC Traces
% 2017-07-24 AL - Added 'DealWithTooShort' as an optional argument
% 2017-07-24 AL - The parameter 'BeforeBreakMs' is changed to 'BeforePeakMs'
% 2017-07-25 AL - Moved code to minEASE_compute_plot_average_psc.m
% 2017-10-16 AL - Added baselineWindowMs, maxBelowBasePerc and minPscDecayMs
% 2018-01-28 AL - Now marks event threshold at peaks and at all peaks
% 2018-01-28 AL - Now plots all breakpoints and peaks
% 2018-01-28 AL - Added markerSize and markerEdgeWidth
% 2018-01-28 AL - Added isdeployed
% 2018-01-29 AL - Added sweepLabel to plot title
% 2018-01-29 AL - Now automatically detects region of interest (roiPlot)
% 2018-02-02 AL - Move message to minEASE.m
% 2018-02-08 AL - Added showMessage and now uses print_or_show_message
% 2018-02-08 AL - Now accepts either minPscAmpAbs or minPscAmpRel
% 2018-02-08 AL - Added Hard-coded constants section
% 2018-02-15 AL - Now plots event detection even if no events are detected
% 2018-02-15 AL - Fixed bug on centering to the peak of an event of interest
%                   Previously it was just centering to the center of the trace
% 2018-02-15 AL - Fixed bug on plotting eventAmpThreshold
% 2018-02-16 AL - Now prints message on how many windows 
%                   were used to compute Gaussian noise
% 2018-02-16 AL - Now plots the direction-filtered data in magenta 
%                   and the Gaussian part of the data in green
% 2018-02-16 AL - Changed skewnessCutoff and excessKurtosisCutoff to 
%                   zSkewnessThres and zExcessKurtosisThres
% 2018-02-27 AL - Changed showMessages to messageMode with possible values:
%                   'wait', 'show', 'none'
% 2018-03-02 MD - Defined verbose parameter for print_or_show_message
% 2018-08-03 AL - Renamed sweepLabel -> outputLabel
% 2018-09-17 AL - Now uses save_all_figtypes.m
% 

%% Hard-coded parameters
validMessageModes = {'wait', 'show', 'none'};
markerSize = 6;                 % marker size for plotting 
                                %   detected breakpoints and peaks
markerEdgeWidth = 2;            % marker edge width for plotting 
                                %   detected breakpoints and peaks

%% Hard-coded constants
%   Used by minEASE.m
TYPEONE_CLASSNUM    = 1;
TYPETWO_CLASSNUM    = 2;
TYPETHREE_CLASSNUM  = 3;
SLOWRISE_CLASSNUM   = 4;
WRONGDECAY_CLASSNUM = 5;
TOOSMALL_CLASSNUM   = 6;
INSEALTEST_CLASSNUM = 7;
REMOVED_CLASSNUM    = 8;

%% Constants to be consistent with find_directional_events.m
IDXBREAK_COLNUM      = 1;
IDXPEAK_COLNUM       = 2;
VALBREAK_COLNUM      = 3;
VALPEAK_COLNUM       = 4;
EVENTAMP_COLNUM      = 5;
TOTALRISE_COLNUM     = 6;
TENNINETYRISE_COLNUM = 7;
IEI_COLNUM           = 8;
ISI_COLNUM           = 9;
HALFDECAY_COLNUM     = 10;
FULLDECAY_COLNUM     = 11;

%% Conversion constants
MS_PER_S = 1e3;                 % milliseconds per second

%% Default flags
plotEventDetectionFlagDefault = true;   % whether to plot event detection
plotAverageTraceFlagDefault = true;     % whether to plot average trace
verboseDefault = false;             % default: Program does not print message
                                    %   even if message box is shown

%% Default parameters for Butterworth lowpass filter
lowpassCutoffDefault = 3000;    % 3 dB cutoff frequency (Hz) of lowpass filter
lowpassNpolesDefault = 8;       % order of Butterworth lowpass filter
                                % (degree of transfer function denominator)

%% Default parameters for computing Gaussian RMS noise
noiseWindowSizeDefault = 100;   % default noise window size (samples)
zSkewnessThresDefault = 0.2;    % default Gaussian noise skewness cutoff
zExcessKurtosisThresDefault = 0.2;  % default noise excess kurtosis cutoff

%% Default parameters for finding directional events
signal2NoiseDefault = 2;        % default signal/noise ratio of events to noise
                                % TODO: why 2?
minEventThresholdDefault = 8;   % default minimum event amplitude threshold (pA)
                                % TODO: why 8 pA?
baselineWindowMsDefault = 5;    % window for computing baseline (ms)
                                % TODO: why 5 ms?
maxBelowBasePercDefault = 100;  % default maximum below baseline percentage (%)
smoothWindowMsDefault = 0.5;    % default moving average filter window (ms)

%% Default parameters for classifying events and identifying PSCs
minPscAmpAbsDefault = 10;       % default minimum allowable PSC amplitude (pA)
                                % TODO: why 10 pA?
minPscAmpRelDefault = 10;       % default relative minimum allowable 
                                %   PSC amplitude (% of maximum)
                                % works well for Paula's calcium imaging data
maxPscRiseMsDefault = 8;        % default maximum allowable PSC rise time (ms)
                                % TODO: why 8 ms?
minPscDecayMsDefault = 0.3;     % default minimum allowable PSC decay time (ms)
                                % TODO: why 0.3 ms?
maxPscDecayMsDefault = 50;      % default maximum allowable PSC decay time (ms)
                                % TODO: why 50 ms?
sealTestWindowMsDefault = [1000, 1050];      
                                % default seal test window (ms)
                                % TODO: why 50 ms?

%% Default parameters for averaging PSC traces
traceLengthMsDefault = 50;      % default total trace length for each PSC (ms)
                                % TODO: why 50 ms?
beforePeakMsDefault = 3;        % default trace length before breakpoint (ms)
                                % TODO: Why 3 ms?

%% Default parameters for plotting event detection
roiPlotDefault = [];            % default region of interest to plot 
                                %   for event detection (ms)

%% Default parameters for saving outputs
outputDirectoryDefault = '/media/shareX/share/minEASE/output/';
figTypesDefault = 'png';        % default figure type(s) for saving
fileNoDefault = 0;              % default file number for saving output
outputLabelDefault = '';        % temporary default output label

%% Other default parameters
messageModeDefault = 'none';    % print to standard output by default
skipDetectionDefault = false;   % whether to skip event detection by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to an Input Parser
addRequired(iP, 'data', ...                     % vector of current data
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'siMs', ...                     % sampling interval in ms/sample
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addRequired(iP, 'directionPsc', ...             % PSC direction to detect
    @(x) any(validatestring(x, {'Excitatory', 'Inhibitory', 'E', 'I'})));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'LowpassCutoff', lowpassCutoffDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'LowpassNpoles', lowpassNpolesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'NoiseWindowSize', noiseWindowSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'ZSkewnessThres', zSkewnessThresDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'ZExcessKurtosisThres', zExcessKurtosisThresDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'Signal2Noise', signal2NoiseDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'MinEventThreshold', minEventThresholdDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'BaselineWindowMs', baselineWindowMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'MaxBelowBasePerc', maxBelowBasePercDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'SmoothWindowMs', smoothWindowMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'MinPscAmpAbs', [], ... 
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'MinPscAmpRel', [], ... 
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'MaxPscRiseMs', maxPscRiseMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'MinPscDecayMs', minPscDecayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'MaxPscDecayMs', maxPscDecayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'SealTestWindowMs', sealTestWindowMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, ...
            {'vector', 'nonnegative', 'numel', 2, 'nondecreasing'}));
addParameter(iP, 'TraceLengthMs', traceLengthMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'BeforePeakMs', beforePeakMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'OutputDirectory', outputDirectoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) min(isfigtype(x, 'ValidateMode', true)));
addParameter(iP, 'FileNo', fileNoDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
addParameter(iP, 'OutputLabel', outputLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(iP, 'PlotEventDetectionFlag', plotEventDetectionFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotAverageTraceFlag', plotAverageTraceFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MessageMode', messageModeDefault, ...
    @(x) any(validatestring(x, validMessageModes)));
addParameter(iP, 'SkipDetection', skipDetectionDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RoiPlot', roiPlotDefault,...
    @(x) isempty(x) || (isnumeric(x) && isvector(x) && all(x >= 0)...
                        && diff(x) >= 0 && numel(x) == 2));

% Read from the Input Parser
parse(iP, data, siMs, directionPsc, varargin{:});
directionPsc           = validatestring(directionPsc, ...
                         {'Excitatory', 'Inhibitory', 'E', 'I'});
lowpassCutoff          = iP.Results.LowpassCutoff;
lowpassNpoles          = iP.Results.LowpassNpoles;
noiseWindowSize        = iP.Results.NoiseWindowSize;
zSkewnessThres         = iP.Results.ZSkewnessThres;
zExcessKurtosisThres   = iP.Results.ZExcessKurtosisThres;
signal2Noise           = iP.Results.Signal2Noise;
minEventThreshold      = iP.Results.MinEventThreshold;
baselineWindowMs       = iP.Results.BaselineWindowMs;
maxBelowBasePerc       = iP.Results.MaxBelowBasePerc;
smoothWindowMs         = iP.Results.SmoothWindowMs;
minPscAmpAbs           = iP.Results.MinPscAmpAbs;
minPscAmpRel           = iP.Results.MinPscAmpRel;
maxPscRiseMs           = iP.Results.MaxPscRiseMs;
minPscDecayMs          = iP.Results.MinPscDecayMs;
maxPscDecayMs          = iP.Results.MaxPscDecayMs;
sealTestWindowMs       = iP.Results.SealTestWindowMs;
traceLengthMs          = iP.Results.TraceLengthMs;
beforePeakMs           = iP.Results.BeforePeakMs;
outputDirectory        = iP.Results.OutputDirectory;
[~, figTypes]          = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
fileNo                 = iP.Results.FileNo;
outputLabel             = iP.Results.OutputLabel;
plotEventDetectionFlag = iP.Results.PlotEventDetectionFlag;
plotAverageTraceFlag   = iP.Results.PlotAverageTraceFlag;
messageMode            = validatestring(iP.Results.MessageMode, ...
                                        validMessageModes);
skipDetection          = iP.Results.SkipDetection;
verbose                = iP.Results.Verbose;
roiPlot                = iP.Results.RoiPlot;       % TODO: Make this an optional argument

% Set dependent argument defaults
if isempty(outputLabel)
    outputLabel = ['Swp', num2str(fileNo)];
end

% Give the absolute amplitude threshold priority
if ~isempty(minPscAmpAbs)       % if absolute amplitude threshold given
    % Print warning message
    if ~isempty(minPscAmpRel)   % if relative amplitude threshold also given
        message = {'Absolute minimum PSC amplitude threshold exists.', ...
                    'Relative minimum PSC amplitude threshold ignored!'};
        mTitle = 'Input Ignored Warning';
        icon = 'warn';
        print_or_show_message(message, 'MessageMode', messageMode, ...
                                'MTitle', mTitle, 'Icon', icon, ...
                                'Verbose', verbose);
    end

    % Remove relative threshold if any
    minPscAmpRel = [];
elseif ~isempty(minPscAmpRel)   % if absolute amplitude threshold not given
                                %   but relative amplitude threshold given
    % Use relative amplitude threshold
else                            % if no amplitude threshold given
    % Use the default absolute amplitude threshold
    minPscAmpAbs = minPscAmpAbsDefault;
end

% Display warning message if some inputs are unmatched
if ~isempty(fieldnames(iP.Unmatched))
    fprintf('WARNING: The following name-value pairs could not be parsed: \n');
    disp(iP.Unmatched);
end

%% Extract info from inputs
nSamples = length(data);                        % number of sample points

% Prepare for plotting if plotEventDetectionFlag is true
if plotEventDetectionFlag
    % Set up time vector (ms)
    tVec = linspace(siMs, nSamples*siMs, nSamples);

    % Find maximum range of x axis
    xRange = [0, tVec(end)];

    % Choose region of interest to plot
    % TODO: Make this a parameter to accept user input
    %       If user input is not provided, try finding a region with 
    %       detected events. Otherwise, use default
    if ~isempty(roiPlot) && ...
        (roiPlot(1) < xRange(1) || roiPlot(end) > xRange(end))
        error('Region of interest cannot be greater than the maximum range!');
    end
end

%% First filter the data with a zero-phase lowpass Butterworth filter
%   with 3 dB cutoff frequency lowpassCutoff
%   and order lowpassNpoles each way (or 2*lowpassNpoles in total)

% Find the normalized cutoff frequency Wn = fc/(fs/2), 
%   where fs = sampling frequency (Hz) = 1/siMs 
%   and fs/2 is the Nyquist frequency
Wn = lowpassCutoff * 2 * (siMs / MS_PER_S);
                        % normalized cutoff frequency (half-cycles/sample)

% Find the transfer function coefficients of a lowpass Butterworth filter
%   with order npoles and normalized cutoff frequency Wn
[numeratorCoeff, denominatorCoeff] = butter(lowpassNpoles, Wn);

% Check the order of the filter
orderFilter = filtord(numeratorCoeff, denominatorCoeff);
if nSamples <= 3 * orderFilter
    error(['Not enough data points to apply a ', ...
            'Butterworth filter of order %d twice!\n'], ...
            orderFilter);
end

% Lowpass-filter data twice (forward & reverse directions)
dataLowpass = filtfilt(numeratorCoeff, denominatorCoeff, data);

%% Find root-mean-square noise of the "Gaussian parts" of dataLowpass
%   Algorithm borrowed from matlabscripts/John
[noiseLevel, indGauss, message, ~, nKeptWindows] = ...
        compute_rms_Gaussian(dataLowpass, ...
                    'WindowSize', noiseWindowSize, ...
                    'ZSkewnessThres', zSkewnessThres, ...
                    'ZExcessKurtosisThres', zExcessKurtosisThres);
mTitle = 'Root-mean square level of Gaussian part';
if nKeptWindows == 0
    icon = 'warn';
else
    % TODO: custom icon
    icon = 'none';
end
print_or_show_message(message, 'MessageMode', messageMode, ...
                        'MTitle', mTitle, 'Icon', icon, 'Verbose', verbose);

% Get the "Gaussian parts" of dataLowpass
if plotEventDetectionFlag
    tVecGauss = tVec(indGauss);
    dataGauss = dataLowpass(indGauss);
end

%% Exit program if skipDetection is true
if skipDetection
    eventInfo = 'skipped';
    eventClass = 'skipped';
    return
end

%% Find all directional events in dataLowpass
% Compute windows in terms of samples
baselineWindowSamples = floor(baselineWindowMs / siMs);
smoothWindowSamples = floor(smoothWindowMs / siMs);

% Decide on the direction of the event
switch directionPsc
case {'Excitatory', 'E'}
    directionEvent = 'Downward';
    directionFactor = -1;
case {'Inhibitory', 'I'}
    directionEvent = 'Upward';
    directionFactor = 1;
otherwise
    error('This PSC direction is not supported!\n');
end

% Use find_directional_events.m to find all directional events
[eventInfo, dataSmooth, dataDirFilt, eventAmpThreshold, ...
    idxBreaks, valBreaks, idxPeaks, valPeaks] = ...
    find_directional_events(dataLowpass, directionEvent, ...
                'NoiseLevel', noiseLevel, ...
                'Signal2Noise', signal2Noise, ...
                'MinAmpThreshold', minEventThreshold, ...
                'BaselineWindow', baselineWindowSamples, ...
                'MaxBelowBasePerc', maxBelowBasePerc, ...
                'SmoothWindow', smoothWindowSamples);

% Check if an event is detected
nEvents = size(eventInfo, 1);           % number of directional events detected
if nEvents == 0
    % Nothing to classify
    eventClass = [];
    idxEventPeak = [];
    idxBreaks = [];
    valBreaks = [];
    idxPeaks = [];
    valPeaks = [];

    % Plot figure for event detection
    if plotEventDetectionFlag
        plot_event_detection (roiPlot, siMs, xRange, outputDirectory, outputLabel, ...
                                markerSize, markerEdgeWidth, figTypes, ...
                                tVec, data, dataLowpass, dataSmooth, ...
                                tVecGauss, dataGauss, dataDirFilt, ...
                                eventAmpThreshold, directionFactor, ...
                                eventInfo, eventClass, idxEventPeak, ...
                                idxBreaks, valBreaks, idxPeaks, valPeaks);
    end

    % Exit program
    return;
end

%% Extract information from eventInfo
% Column assignments for eventInfo:
%   1 = index at event breakpoint
%   2 = index at event peak
%   3 = data value at event breakpoint
%   4 = data value at event peak
%   5 = amplitude of the event
%   6 = 0-100% rise time (samples)
%   7 = 10-90% rise time (samples)
%   8 = inter-event interval from this event 
%       peak to next event peak (samples)
%   9 = interstimulus interval from this event 
%       peak to next event breakpoint (samples)
%   10 = 50% decay time (samples)
%   11 = "full decay" time (samples):
%           time to return within noiseLevel 
%           of breakpoint value
idxEventBreak         = eventInfo(:, IDXBREAK_COLNUM);
idxEventPeak          = eventInfo(:, IDXPEAK_COLNUM);
valEventBreak         = eventInfo(:, VALBREAK_COLNUM);
valEventPeak          = eventInfo(:, VALPEAK_COLNUM);
eventAmplitude        = eventInfo(:, EVENTAMP_COLNUM);
totalRiseTime         = eventInfo(:, TOTALRISE_COLNUM);
tenNinetyRiseTime     = eventInfo(:, TENNINETYRISE_COLNUM);
interEventInterval    = eventInfo(:, IEI_COLNUM);
interStimulusInterval = eventInfo(:, ISI_COLNUM);
halfDecayTime         = eventInfo(:, HALFDECAY_COLNUM);
fullDecayTime         = eventInfo(:, FULLDECAY_COLNUM);

%% Classify events and identify all post-synaptic currents (EPSCs or IPSCs)
%   The series of operations will try to pick out just "type I events"
%   through kinetics and amplitude assay
eventClass = zeros(nEvents, 1);                 % initialize class vector

% Omit events with amplitude that are too small
if ~isempty(minPscAmpAbs)
    % Omit events with amplitude less than minPscAmpAbs
    isTooSmall = isnan(eventAmplitude) | eventAmplitude < minPscAmpAbs;
elseif ~isempty(minPscAmpRel)
    % Find the largest event amplitude
    largestAmplitude = max(eventAmplitude);

    % Omit events with amplitude less than minPscAmpRel % of largest amplitude
    isTooSmall = isnan(eventAmplitude) | ...
                    eventAmplitude < minPscAmpRel * largestAmplitude / 100;
end
eventInfoTooSmall = eventInfo(isTooSmall, :);   % all events that are too small
eventClass(isTooSmall & eventClass == 0) = TOOSMALL_CLASSNUM;
                                                % CLASS 6 events are too small
                                                %   and possibly slow rise and 
                                                %   slow decay

% Omit events with 10-90% rise times greater than maxPscRiseMs
isSlowRise = isnan(tenNinetyRiseTime) | tenNinetyRiseTime > maxPscRiseMs/siMs;
eventInfoSlowRise = eventInfo(isSlowRise, :);   % all events that are slow rise
eventClass(isSlowRise & eventClass == 0) = WRONGDECAY_CLASSNUM;
                                                % CLASS 5 events are slow rise
                                                %   and possibly slow decay

% Omit all events with half decay times smaller than minPscDecayMs 
%       or greater than maxPscDecayMs
isWrongDecay = isnan(halfDecayTime) | ...
                    halfDecayTime < minPscDecayMs/siMs | ...
                    halfDecayTime > maxPscDecayMs/siMs;
eventInfoWrongDecay = eventInfo(isWrongDecay, :);   
                                                % all events that are wrong decay
eventClass(isWrongDecay & eventClass == 0) = SLOWRISE_CLASSNUM;
                                                % CLASS 4 events are wrong decay

% Find incomplete events:
%       events that have not decayed fully
%           i.e.,
%       events with interstimulus interval (peak to next breakpoint)
%           less than the "full decay time"  
%   TODO: may need to change to "PSC full decay time"
%   Note: In this case, fullDecayTime cannot be calculated so is NaN
isIncomplete = isnan(fullDecayTime);
isIncomplete(end) = false;                      % the completeness of the last 
                                                %   event cannot be determined
eventInfoIncomplete = eventInfo(isIncomplete, :);   % all incomplete events

% Find premature events:
%       events that begin before the preceding event has decayed fully,
%           i.e.,
%       events with preceding interstimulus interval (prev peak to breakpoint)
%           less than the preceding "full decay time"  
%   TODO: may need to change to preceding "PSC full decay time"
%   Note: In this case, the preceding fullDecayTime is NaN
isPremature = circshift(isnan(fullDecayTime), 1);   % shift to right by 1
isPremature(1) = false;                         % the first event cannot 
                                                %   be premature
eventInfoPremature = eventInfo(isPremature, :); % all premature events

% Omit Type II PSCs:
%       incomplete events that are not also premature events
%       and not filtered out by amplitude, rise time nor decay time criteria
eventClass(isIncomplete & ~isPremature & eventClass == 0) = TYPETWO_CLASSNUM; 
                                                % CLASS 2 events = Type II PSCs

% Omit Type III PSCs:
%       premature events that follow Type II PSCs or other Type III PSCs
%       and not filtered out by amplitude, rise time nor decay time criteria
for iEvent = 2:nEvents
    if isPremature(iEvent) && eventClass(iEvent) == 0 && ...
        (eventClass(iEvent - 1) == TYPETWO_CLASSNUM || ...
            eventClass(iEvent - 1) == TYPETHREE_CLASSNUM)
        eventClass(iEvent) = TYPETHREE_CLASSNUM;
                                                % CLASS 3 events = Type III PSCs
    end
end

% Omit events within the seal test window
idxSealTestOn = round(sealTestWindowMs(1)/siMs);
idxSealTestOff = round(sealTestWindowMs(2)/siMs);
isInSealTest = idxEventPeak >= idxSealTestOn & idxEventPeak <= idxSealTestOff;
eventInfoInSealTest = eventInfo(isInSealTest, :);
                            % all events that are within the seal test window
eventClass(isInSealTest) = INSEALTEST_CLASSNUM; % CLASS 7 events are within 
                                                %   the seal test window

% Check if there are any events left and make them PSCs
isTypeOne = eventClass == 0;
typeOneInfo = eventInfo(isTypeOne, :);
eventClass(isTypeOne) = TYPEONE_CLASSNUM;       % CLASS 1 events are Type I PSCs

%% Plot figure for event detection
if plotEventDetectionFlag
    plot_event_detection (roiPlot, siMs, xRange, outputDirectory, outputLabel, ...
                            markerSize, markerEdgeWidth, figTypes, ...
                            tVec, data, dataLowpass, dataSmooth, ...
                            tVecGauss, dataGauss, dataDirFilt, ...
                            eventAmpThreshold, directionFactor, ...
                            eventInfo, eventClass, idxEventPeak, ...
                            idxBreaks, valBreaks, idxPeaks, valPeaks);
end

%% TODO: Visualize the classes

%% Detect bursts
% TODO: Put Type II & Type III events together to make bursts and count them

%% Compute mean data
% meanAmp = mean(pscInfo(:, 5));
% meanRise = mean(pscInfo(:, 7) * siMs);
% meanHW = mean(pscInfo(:, 9) * siMs);

%% Compute burst statistics  %% TODO: What for?
%{
justBursts = eventInfoExpanded(eventInfoExpanded(:, 10) ~= 0, :);
                                            % information for events in bursts
if ~isempty(justBursts)
    % Compute time that this cell spent bursting (ms),
    %   i.e. sum of IEI between spikes
    time_in_burst = sum(justBursts(:, 8)); % sum of the iei % TODO: Shouldn't included IEI of last spike
    time_in_burst = time_in_burst * siMs; % sum of the iei (in ms)

    % Percentage of the total time of the recording that the cell was bursting
    percent_in_burst = time_in_burst / (nSamples*siMs); % percentage of total (ms)
    amp_in_burst = sum(justBursts(:, 5));   % TODO: What does the total amplitude sum mean?
    number_of_burst_perMin = max(justBursts(:,10)) / (( nSamples*siMs / 1000) / 60);
else
    time_in_burst = 0;
    percent_in_burst = 0;
    amp_in_burst = 0;
    number_of_burst_perMin = 0;
end
%}

% [aaa bbb] = fileparts(abf_file);
% display(sprintf('%s Means: Rise: %f Amp: %f, HW: %f',...
%     bbb, meanRise, meanAmp, meanHW));
% display(sprintf('%s Bursts: Percent in burst: %f, Cumulative Amp in Burst: %f',...
%     bbb, percent_in_burst, amp_in_burst));

%% Plot figure for burst detection
%{ 
TODO
if plotBurstDetectionFlag
    fig2 = figure(2);
    clf(fig2);
    title('Burst Detection', 'FontSize', 14);
    hold on;
    
    % Plot original data in black
    plot(tVec, data, 'k');

    % Plot all event breakpoints in blue
    plot(eventInfo(:, 1) * siMs, eventInfo(:, 3), 'b.');

    % Plot event breakpoints in "crude burst regions" in green
    plot(eventInfoExpanded(eventInfoExpanded(:, 11) > 0, 1) * siMs, ...
            eventInfoExpanded(eventInfoExpanded(:, 11) > 0, 3), ...
            'g.');

    % Plot event breakpoints in bursts in yellow
    plot(eventInfoExpanded(eventInfoExpanded(:, 10) > 0, 1) * siMs, ...
            eventInfoExpanded(eventInfoExpanded(:, 10) > 0, 3), ...
            'y.');

    % Set axis limits to be the range of data
    axis tight

    % Save figure
    figName = fullfile(outputDirectory, ['Burst_detection_', outputLabel]);
    save_all_figtypes(fig2, figName, figTypes);

    % Zoom to region of interest
    xlim(ROI_PLOT);
    ylim('auto');

    % Save figure
    figName = fullfile(outputDirectory, ...
            ['Burst_detection_', outputLabel, '_zoom']);
    save_all_figtypes(fig2, figName, figTypes);

    % Close figure
    % close(fig2);
end
%}

%% Plot figure for averaged PSC traces
%{
if plotAverageTraceFlag
    minEASE_compute_plot_average_psc(eventInfo, eventClass, data, siMs, ...
                                    traceLengthMs, beforePeakMs, 'none', ...
                                    outputDirectory, outputLabel, ...
                                    'FigTypes', figTypes);
    minEASE_compute_plot_average_psc(eventInfo, eventClass, data, siMs, ...
                                    traceLengthMs, beforePeakMs, 'padright', ...
                                    outputDirectory, outputLabel, ...
                                    'FigTypes', figTypes);
    minEASE_compute_plot_average_psc(eventInfo, eventClass, data, siMs, ...
                                    traceLengthMs, beforePeakMs, 'padboth', ...
                                    outputDirectory, outputLabel, ...
                                    'FigTypes', figTypes);
    minEASE_compute_plot_average_psc(eventInfo, eventClass, data, siMs, ...
                                    traceLengthMs, beforePeakMs, 'omit', ...
                                    outputDirectory, outputLabel, ...
                                    'FigTypes', figTypes);
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_event_detection (roiPlot, siMs, xRange, outputDirectory, outputLabel, markerSize, markerEdgeWidth, figTypes, tVec, data, dataLowpass, dataSmooth, tVecGauss, dataGauss, dataDirFilt, eventAmpThreshold, directionFactor, eventInfo, eventClass, idxEventPeak, idxBreaks, valBreaks, idxPeaks, valPeaks)

% If region of interest not given, try to provide one
if isempty(roiPlot)        
    % Decide on the center of the region of interest
    center = [];                % initialize center
    for iClass = 1:7
        % If there is an event of this class, use the time of the middle 
        %   event and stop (Class 1 has priority over Class 2, and so on.)
        isThisClass = eventClass == iClass;
        if any(isThisClass)         % if there is any event in this class
            % Find the median number of total number of events in this class
            halfNThisClass = round(sum(isThisClass)/2);

            % Get all the event number in eventInfo for this class
            thisEventNumbers = find(isThisClass);

            % Find the event number for the median event
            medianEventNumber = thisEventNumbers(halfNThisClass);
            
            % Find the time corresponding to this event peak
            center = idxEventPeak(medianEventNumber) * siMs;
            break;
        end
    end

    % If no events are detected, use the center of the trace
    if isempty(center)
        maxROI = xRange(end) - xRange(1);
        center = xRange(1) + maxROI/2;
    end

    % Find a 100 ms region around the ROI center
    %   Make sure the region of interest doesn't exceed boundaries
    roiPlot = [max(center - 50, xRange(1)), ...
                min(center + 50, xRange(end))];
end

% Create figure
fig1 = figure(1);
clf(fig1);

% Create first subplot
fig1Sp1Ax = subplot(4, 3, [1, 3]);
title(['Data Processing for ', outputLabel], 'FontSize', 14, ...
        'interpreter', 'none');
hold on;

% Plot original data in black
plot(tVec, data, '-k', 'LineWidth', 0.5);

% Plot lowpass filtered data in blue
plot(tVec, dataLowpass, '-b', 'LineWidth', 0.5);

% Plot moving-average-filtered data in magenta
plot(tVec, dataSmooth, '-m', 'LineWidth', 0.5);

% Plot Gaussian part of the data in green
%   Can't plot with a line, because the Gaussian windows might not be continuous
plot(tVecGauss, dataGauss, '.g', 'MarkerSize', 3);

% Plot event amplitude threshold (pA) as a dotted line in cyan
if ~isempty(eventInfo)
    plot(tVec(idxPeaks), ...
         valBreaks + eventAmpThreshold * directionFactor, ...
         'Color', [0 1 1], 'LineStyle', ':', 'LineWidth', 0.5);
end
ylim('auto');

% Remove X tick label to save space
set(gca, 'XTickLabel', []);

% Create second subplot
fig1Sp2Ax = subplot(4, 3, [4, 6]);
hold on;

% Plot direction-filtered data in magenta
plot(tVec, dataDirFilt, '-m', 'LineWidth', 0.5);

% Plot event amplitude threshold (pA) as a dotted line in cyan
line(xRange, [eventAmpThreshold, eventAmpThreshold] * directionFactor, ...
     'Color', [0 1 1], 'LineStyle', ':', 'LineWidth', 0.5);

% Create third subplot
fig1Sp3Ax = subplot(4, 3, [7, 12]);
minEASE_plot_gapfree_events(tVec, data, dataLowpass, ...
                    eventInfo, eventClass, outputLabel, ...
                    'DataSmooth', dataSmooth, ...
                    'MarkerSize', markerSize, ...
                    'LineWidth', markerEdgeWidth, ...
                    'NoTitleFlag', true);

% Link the x axes for the 3 subplots
linkaxes([fig1Sp1Ax, fig1Sp2Ax, fig1Sp3Ax], 'x');

% Save figure
figName = fullfile(outputDirectory, ['Event_detection_', outputLabel]);
save_all_figtypes(fig1, figName, figTypes);

% Zoom to region of interest
xlim(roiPlot);
ylim('auto');

% Also plot all breakpoints (black crosses) and peaks (red circles)
subplot(4, 3, [1, 3]);          % go to first subplot
if ~isempty(eventInfo)
    plot(tVec(idxBreaks), valBreaks, 'kx', ...
         'MarkerSize', markerSize, 'LineWidth', markerEdgeWidth);
    plot(tVec(idxPeaks), valPeaks, 'ro', ...
         'MarkerSize', markerSize, 'LineWidth', markerEdgeWidth);
end

% Save figure
figName = fullfile(outputDirectory, ['Event_detection_', outputLabel, '_zoom']);
save_all_figtypes(fig1, figName, figTypes);

% Close figure
% close(fig1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%{
OLD CODE:

dataLowpass = filter(numeratorCoeff, denominatorCoeff, data);
dataLowpass = zof_mark(data, lowpassCutoff, lowpassNpoles, siMs/MS_PER_S);

%                   - 'MinBaselineDiff': minimum baseline difference (pA)
%                                        for identifying a "crude burst region"
%                   must be a positive scalar
%                   default == 40
%                   - 'CrudeRegionSize': crude burst region size (events)
%                   must be a positive integer scalar
%                   default == 50
%                   - 'MinSpikesPerBurst': minimum number of spikes in a burst
%                   must be a positive integer scalar
%                   default == 8
%                   - 'MaxIsiMs': maximum inter-spike interval (ms) 
%                   must be a positive scalar
%                   default == 50
%% Default parameters for identifying bursts
minBaselineDiffDefault = 40;    % default minimum baseline difference (pA)
                                %   for identifying a crude burst region
                                % TODO: Why 40 pA?
crudeRegionSizeDefault = 50;    % default burst region size (events)
                                % TODO: Why 50 spikes?
minSpikesPerBurstDefault = 8;   % default minimum number of spikes in a burst
                                % TODO: Why 8 spikes?
maxIsiMsDefault = 50;           % default maximum inter-spike interval (ms) 
                                %   within a burst
                                % TODO: Why 50 ms?
addParameter(iP, 'MinBaselineDiff', minBaselineDiffDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'CrudeRegionSize', crudeRegionSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'MinSpikesPerBurst', minSpikesPerBurstDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'MaxIsiMs', maxIsiMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
minBaselineDiff      = iP.Results.MinBaselineDiff;
crudeRegionSize      = iP.Results.CrudeRegionSize;
minSpikesPerBurst    = iP.Results.MinSpikesPerBurst;
maxIsiMs             = iP.Results.MaxIsiMs;

% Determine whether directional events are in bursts:
%   bursts should have lots of PSCs with high variability 
%   in their baseline currents
[eventInfoExpanded] = ...
    identify_bursts(eventInfo, 'MinBaselineDiff', minBaselineDiff, ...
                                'CrudeRegionSize', crudeRegionSize, ...
                                'MinSpikesPerBurst', minSpikesPerBurst, ...
                                'MaxIsi', floor(maxIsiMs / siMs));

    % Plot event breakpoints in "crude burst regions" in green
    plot(eventInfoExpanded(eventInfoExpanded(:, 11) > 0, 1) * siMs, ...
            eventInfoExpanded(eventInfoExpanded(:, 11) > 0, 3), ...
            'g.');

    % Display text for region numbers
    % TODO

    % Plot event breakpoints in bursts in yellow
    plot(eventInfoExpanded(eventInfoExpanded(:, 10) > 0, 1) * siMs, ...
            eventInfoExpanded(eventInfoExpanded(:, 10) > 0, 3), ...
            'y.');

% TODO: Make the following into a function
%       [pscInfo] = identify_PSCs (eventInfoExpanded);
%       or [pscInfo] = identify_PSCs (eventInfo, inBurst);

% Omit all events that are in bursts
%   Initialize pscInfo as the events that are NOT in bursts
%   Note: bursts will interfere with the kinetics and amplitude assay
eventInfoInBurst = eventInfo(eventInfoExpanded(:, 10) == 1, :);
pscInfo = eventInfo(eventInfoExpanded(:, 10) == 0, :);

% Omit all events with amplitude less than minPscAmp
isTooSmall = pscInfo(:, 5) < minPscAmp;
if ~isempty(isTooSmall)
    % Store small events in its own matrix eventInfoTooSmall
    eventInfoTooSmall = pscInfo(isTooSmall, :);

    % Delete small events from pscInfo
    pscInfo(isTooSmall, :) = [];
end

% Omit all events with 10-90% rise times greater than maxPscRiseMs
isSlowRise = pscInfo(:, 7) > maxPscRiseMs/siMs;
if ~isempty(isSlowRise)
    % Store small events in its own matrix eventInfoSlowRise
    eventInfoSlowRise = pscInfo(isSlowRise, :);
    
    % Delete slow events from pscInfo
    pscInfo(isSlowRise, :) = [];
end

% Omit all events with half decay times greater than the following inter-event 
%   intervals
isIncomplete = pscInfo(:, 9) > pscInfo(:, 8);
if ~isempty(isIncomplete)
    % Store incomplete events in its own matrix eventInfoIncomplete
    eventInfoIncomplete = pscInfo(isIncomplete, :);
    
    % Delete incomplete events from pscInfo
    pscInfo(isIncomplete, :) = [];
end

% Omit all events with half decay times greater than maxPscDecayMs
isSlowDecay = pscInfo(:, 9) > maxPscDecayMs/siMs;
if ~isempty(isSlowDecay)
    % Store slow decay events in its own matrix eventInfoSlowDecay
    eventInfoSlowDecay = pscInfo(isSlowDecay, :);
    
    % Delete slow decay events from pscInfo
    pscInfo(isSlowDecay, :) = [];
end

% Omit all events within the seal test window
idxSealTestOn = round(sealTestWindowMs(1)/siMs);
idxSealTestOff = round(sealTestWindowMs(2)/siMs);
isInSealTest = pscInfo(:, 2) >= idxSealTestOn & ...
                    pscInfo(:, 2) <= idxSealTestOff;
if ~isempty(isInSealTest)
    % Store slow decay events in its own matrix eventInfoSlowDecay
    eventInfoInSealTest = pscInfo(isInSealTest, :);
    
    % Delete slow decay events from pscInfo
    pscInfo(isInSealTest, :) = [];
end

eventClass(isPremature & eventClass == 0) = 3;  % CLASS 3 events = Type III PSCs

    % Plot event amplitude threshold (pA) as a cyan line
    line([tVec(1), tVec(end)], ...
            [meanLowpass + eventAmpThreshold, ...
             meanLowpass + eventAmpThreshold] * directionFactor, ...
            'Color', [0 1 1], 'LineStyle', '--', 'LineWidth', 0.5);

%                   - 'BeforeBreakMs': trace length before break point 
%                                       for each PSC (ms)
%                   must be a numeric positive scalar
%                   default == 3
beforeBreakMsDefault = 3;       % default trace length before breakpoint (ms)
                                % TODO: Why 3 ms?
if ~isempty(typeOneInfo)
    [averageTypeOneTrace, allTypeOneTraces] = ...
        compute_average_PSC_trace(typeOneInfo, data, 'SiMs', siMs, ...
                                    'TraceLengthMs', traceLengthMs, ...
                                    'BeforeBreakMs', beforeBreakMs);
end

% Omit all events with half decay times greater than maxPscDecayMs
isSlowDecay = halfDecayTime > maxPscDecayMs/siMs;
eventInfoSlowDecay = eventInfo(isSlowDecay, :); % all events that are slow decay
eventClass(isSlowDecay & eventClass == 0) = 4;  % CLASS 4 events are slow decay

ROI_PLOT = [1000 1500];         % region of interest to plot (ms)
ROI_PLOT = [5200 5300];         % default region of interest to plot 

    plot(tVec(idxEventBreak), ...
         valEventBreak + eventAmpThreshold * directionFactor, ...
         'Color', [0 1 1], 'LineStyle', '--', 'LineWidth', 0.5);

    title('Data Processing with Cumulative Difference', 'FontSize', 14);

% May need to include a "PSC full decay time" that ignores non-PSCs
%       The question is: is a non-PSC event a signal or a noise?

%                   - 'MinPscAmp': minimum allowable PSC amplitude (pA)
%                   must be a nonnegative scalar
%                   default == 10 pA
minPscAmpDefault = 10;          % default minimum allowable PSC amplitude (pA)
                                % TODO: why 10 pA?
addParameter(iP, 'MinPscAmp', minPscAmpDefault, ... 
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
minPscAmp              = iP.Results.MinPscAmp;

% Omit events with amplitude less than minPscAmp
isTooSmall = eventAmplitude < minPscAmp;

%   TODO: may need to compute a "PSC full decay time" 
%           and set to NaN for classes 4~8

            center = idxEventBreak(find(isThisClass) == halfNThisClass) ...
                        * siMs;

if isnan(noiseLevel)
    fprintf('WARNING: noiseLevel is NaN!!\n');
end

fprintf('%s\n', strjoin(message, '\n'));

if any(strfind(message, 'no window'))

if plotEventDetectionFlag
else
    [noiseLevel, ~, message, ~, nKeptWindows] = ...
            compute_rms_Gaussian(dataLowpass, ...
                        'WindowSize', noiseWindowSize, ...
                        'ZSkewnessThres', zSkewnessThres, ...
                        'ZExcessKurtosisThres', zExcessKurtosisThres);
end

    fprintf('Skipped detection!\n');

print_or_show_message(showMessage, message, 'MTitle', mTitle, 'Icon', icon);

%                   - 'ShowMessage': whether to show messages in messages boxes
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
showMessageDefault  = false;        % print to standard output by default
addParameter(iP, 'ShowMessage', showMessageDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
showMessage = iP.Results.ShowMessage;
        if showMessage
            print_or_show_message(message, 'MessageMode', 'show', ...
                                    'MTitle', mTitle, 'Icon', icon);
        else
            print_or_show_message(message, 'MessageMode', 'none', ...
                                    'MTitle', mTitle, 'Icon', icon);
        end
roiPlot = roiPlotDefault; 

    @(x) isempty(x) || (isnumeric(x) && isvector && all(x >= 0)...
                        && all(diff(x) >= 0) && numel(x) == 2);

saveas(fig1, fullfile(outputDirectory, ...
        ['Event_detection_', outputLabel]), figTypes);
saveas(fig1, fullfile(outputDirectory, ...
        ['Event_detection_', outputLabel, '_zoom']), figTypes);

saveas(fig2, fullfile(outputDirectory, ...
        ['Burst_detection_', outputLabel]), figTypes);
saveas(fig2, fullfile(outputDirectory, ...
        ['Burst_detection_', outputLabel, '_zoom']), figTypes);

%}

