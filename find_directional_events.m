function [eventInfo, dataSmooth, dataDirFilt, eventAmpThreshold, ...
                idxBreaks, valBreaks, idxPeaks, valPeaks] = ...
                    find_directional_events (data, direction, varargin)
%% Find all directional events in a data vector
% Usage: [eventInfo, dataSmooth, dataDirFilt, eventAmpThreshold, ...
%               idxBreaks, valBreaks, idxPeaks, valPeaks] = ...
%                   find_directional_events (data, direction, varargin)
%
% Explanation: 
%       Apply a moving average filter to a time series data, then use 
%       the Cohen & Miles method of a cumulative difference function to 
%       obtain the direction-filtered trace for estimating breakpoints and peaks
%       Only leave events that are not "buried in the baseline"
%
% Outputs: 
%       eventInfo   - information for found directional events
%                       Each row corresponds to a directional event
%                       Column assignments for eventInfo:
%                           1 = index at event breakpoint
%                           2 = index at event peak
%                           3 = data value at event breakpoint
%                           4 = data value at event peak
%                           5 = amplitude of the event
%                           6 = 0-100% rise time (samples)
%                           7 = 10-90% rise time (samples)
%                           8 = inter-event interval from this event peak 
%                               to next event peak (samples)
%                           9 = interstimulus interval from this event peak 
%                               to next event breakpoint (samples)
%                           10 = 50% decay time (samples)
%                           11 = "full decay" time (samples):
%                                   time to return within noiseLevel 
%                                   of breakpoint value
%       dataSmooth  - moving-average-filtered version of data
%       dataDirFilt - direction-filtered version of dataSmooth
%       eventAmpThreshold - amplitude threshold (data units) 
%                           used for detecting events
%
% Arguments:
%       data        - vector of current trace data
%                   must be a numeric vector
%       direction   - direction of post-synaptic current to detect
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'Downward'  - downward peaks (e.g., EPSCs)
%                       'Upward'    - upward peaks (e.g., IPSCs)
%       varargin    - 'NoiseLevel': Noise level for the signal to compare with
%                   must be a positive scalar
%                   default == root-mean-square level of Gaussian part of data
%                   - 'Signal2Noise': Signal/noise ratio of events to noise
%                   must be a positive scalar
%                   default == 2
%                   - 'MinAmpThreshold': minimum event amplitude threshold 
%                                       (same units as data vector)
%                   must be a positive scalar
%                   default == 8 data units
%                   - 'BaselineWindow': window for computing baseline (samples)
%                   must be an integer scalar, will be changed to 1 if <= 0
%                   default == 50 samples
%                   - 'MaxBelowBasePerc': maximum below baseline percentage (%)
%                   must be a nonnegative scalar
%                   default == 100 %
%                   - 'SmoothWindow': moving average filter window (samples)
%                   must be an integer scalar, will be changed to 1 if <= 0
%                   default == 5 samples
%
%
% Requires:
%       /home/Matlab/Kojis_Functions/compute_rms_Gaussian.m
%       cd/find_custom.m
%       cd/adjust_peaks.m
%
% Used by:
%       cd/minEASE_detect_gapfree_events.m

% File History:
% ---------- Created by Koji Takahashi & Mark P Beenhakker
% 2017-05-24 AL - Renamed freePSCdetect.m -> find_PSCs.m
% 2017-05-24 AL - Now uses smooth() instead of conv() for smoothing
% 2017-05-24 AL - event breakpoints and peaks are now expressed in samples
% 2017-05-24 AL - Now accepts 'Excitatory' and 'Inhibitory' as directions too
% 2017-05-24 AL - Now uses compute_rms_Gaussian.m to compute default noise level
% 2017-05-25 AL - Added input parser scheme
% 2017-05-25 AL - Organized, simplified and annotated code
% 2017-05-25 AL - Renamed find_PSCs.m -> find_directional_events.m
% 2017-05-24 AL - Generalized to directional events
% 2017-05-24 AL - Direction is now 'Upward' or 'Downward'
% 2017-06-02 AL - Renamed dataFilt -> dataDirFilt
% 2017-06-05 AL - data values and amplitudes are now in terms of the 
%                   original data trace
% 2017-06-05 AL - Now computes the half decay time within the 
%                   inter-stimulus interval instead of the inter-event interval
% 2017-06-05 AL - Added "full decay times"
% 2017-06-07 AL - Fixed bounds for adjust_peaks
% 2017-06-09 AL - Fixed decay times by adding " - 1"
% 2017-06-13 AL - Changed dataDirFilt to data in the computation of 
%                   10-90% rise times
% 2017-06-13 AL - Now uses find_custom in the computation of 
%                   10-90% rise times
% 2017-07-24 AL - 10% rise times are now computed as the 
%                   'last' index 'less' than 10%
% 2017-07-24 AL - Now uses valEventBreaks & eventAmplitudes 
%                   to compute rise10Vals & rise90Vals
% 2017-10-16 AL - Now compares with baseline and allows directional events
%                   to have a maximal below baseline percentage
% 2018-01-28 AL - Now also returns idxBreaks, valBreaks, idxPeaks, valPeaks
% 2018-01-28 AL - Added isdeployed
% 2018-02-08 AL - Added Hard-coded constants section
% 2018-12-18 AL - Now uses arrayfun in place of cellfun for speed
%

%% Hard-coded constants
% For arguments
validDirections = {'Upward', 'Downward'};

%   Used by minEASE.m, minEASE_combine_events.m, minEASE_detect_gapfree_events.m
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

%% Parameters for computing threshold
signal2NoiseDefault = 2;    % default signal/noise ratio of events to noise
                            % TODO: why 2?
minAmpThresholdDefault = 8; % default minimum amplitude threshold (data units)
                            % TODO: 8 pA? why 8?

%% Parameters for comparing with baseline
baselineWindowDefault = 50;   % window for computing baseline (samples)
maxBelowBasePercDefault = 100;  % default maximum below baseline percentage (%)

%% Parameters for smoothing
smoothWindowDefault = 5;    % default moving average filter window (samples)
                            %   Note: the larger the smooth window, 
                            %   the more filtering and degradation of signal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if ~isdeployed
    if exist('/home/Matlab/', 'dir') == 7
        functionsdirectory = '/home/Matlab/';
    elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
        functionsdirectory = '/scratch/al4ng/Matlab/';
    else
        error('Valid functionsdirectory does not exist!');
    end
    addpath(fullfile(functionsdirectory, '/Kojis_Functions/'));    
                                        % for compute_rms_Gaussian.m
    addpath(fullfile(functionsdirectory, '/Adams_Functions/'));    
                                        % for find_custom.m
end

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error('Not enough input arguments, type ''help find_directional_events'' for usage');
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'find_directional_events';

% Add required inputs to an Input Parser
addRequired(iP, 'data', ...                     % vector of current data
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'direction', ...                % event direction to detect
    @(x) any(validatestring(x, validDirections)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'NoiseLevel', [], ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'Signal2Noise', signal2NoiseDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'MinAmpThreshold', minAmpThresholdDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'BaselineWindow', baselineWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
addParameter(iP, 'MaxBelowBasePerc', maxBelowBasePercDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'SmoothWindow', smoothWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));

% Read from the Input Parser
parse(iP, data, direction, varargin{:});
direction = validatestring(direction, validDirections);
noiseLevel       = iP.Results.NoiseLevel;
signal2Noise     = iP.Results.Signal2Noise;
minAmpThreshold  = iP.Results.MinAmpThreshold;
baselineWindow   = iP.Results.BaselineWindow;
maxBelowBasePerc = iP.Results.MaxBelowBasePerc;
smoothWindow     = iP.Results.SmoothWindow;

%% Based on direction provided, find the direction factor
switch direction
case 'Upward'                   % if signal moves upward (e.g., IPSC)
    directionFactor = 1;
case 'Downward'                 % if signal moves downward (e.g., EPSC)
    directionFactor = -1;
otherwise
    error('This event direction is not supported!\n');
end

%% Extract info from data vector
nSamples = length(data);                        % number of sample points

%% Smooth the data with a moving average filter
%   (This will potentially remove tiny local minimums and maximums
%       and generate a more robust baseline)             % TODO: Is this why?

% Make sure the moving average filter is positive
%   This is necessary because of potential rounding issues
if smoothWindow <= 0
    smoothWindow = 1;
    fprintf('Warning: moving average filter window must be positive\n');
    fprintf('         It''s now changed to 1!\n\n');
end

% Apply the default smooth() function, returned vector will be a column
dataSmooth = smooth(data, smoothWindow);        % smoothed data (column)

%% Implement the event detection algorithm used by Cohen & Miles, 2000
% Find the difference vector of the smoothed data
dataSmoothedDiff = diff(dataSmooth);            % difference of smoothed data

% Based on event direction provided, 
%   direction-filter trace to capture downward/upward going signals only
dataDirFilt = zeros(size(dataSmooth));          % direction-filtered data
                                                % Note: the first element is 
                                                %       always zero
for i = 2:nSamples          % start from second sample
    if dataSmoothedDiff(i - 1) * directionFactor > 0  
                            % if signal moved in the given direction
        % Replace that spot in the zero array with the cumulative difference
        dataDirFilt(i) = dataDirFilt(i - 1) + dataSmoothedDiff(i - 1);
    end
end

% Find indices of all potential event breakpoints in the direction-filtered 
%   trace by locating the zero elements that are followed by a non-zero element
%   Note: (1) the very last element will never be a potential event breakpoint
%         (2) the first peak is always preceded by a breakpoint because
%             the first element of dataDirFilt is always zero
idxBreaksFilt = find(~dataDirFilt(1:nSamples-1) & dataDirFilt(2:nSamples));

% Find indices of all potential event peaks in the direction-filtered trace 
%   by locating the nonzero elements that are followed by a zero element
%   Note: (1) the very first element will never be a potential event peak
%             because it's always zero
%         (2) the very last element is a potential event peak as long as
%             it is not zero
idxPeaksFilt = find(dataDirFilt(1:nSamples) & ~[dataDirFilt(2:nSamples); 0]);

% Make sure idxBreaksFilt & idxPeaksFilt have the same length
if length(idxBreaksFilt) ~= length(idxPeaksFilt)
    error(['There''s a problem with the algorithm: ', ...
            'idxBreaksFilt and idxPeaksFilt have different lengths!\n']);
end

% Make sure idxBreaksFilt is never after idxPeaksFilt
if ~isempty(find(idxBreaksFilt > idxPeaksFilt, 1))
    disp(find(idxBreaksFilt > idxPeaksFilt));    
    error(['There''s a problem with the algorithm: ', ...
            'idxBreaksFilt is after idxPeaksFilt at the above indices!\n']);
end

%% Apply a signal-to-noise ratio (amplitude) criterion 
%   on the original data trace to detect events

% If no noiseLevel provided, use the root mean square of the 
%   Gaussian part of the original data trace as the default
if isempty(noiseLevel)
    noiseLevel = compute_rms_Gaussian(data);            % default Gaussian noise level
end

% Compute the amplitude threshold for event detection
%   i.e., this is the amplitude of an event from the last "trough"
eventAmpThreshold = signal2Noise * noiseLevel;  % event amplitude threshold

% Change the threshold to match the required minimum if necessary
if eventAmpThreshold < minAmpThreshold
    eventAmpThreshold = minAmpThreshold;
end

% Find the midpoints between the previous peaks (or 1) and breakpoints
idxMidPrevPeaks2Breaks = ([1; idxPeaksFilt(1:end-1)] + idxBreaksFilt) / 2;

% Find the midpoints between the breakpoints and peaks
idxMidBreaks2Peaks = (idxBreaksFilt + idxPeaksFilt) / 2;

% Find the midpoints between the peaks and next breakpoints (or nSamples)
idxMidPeaks2NextBreaks = (idxPeaksFilt + [idxBreaksFilt(2:end); nSamples]) / 2;

% For all direction-filtered events, 
%   find the corresponding breakpoint values in the original trace
%   Note: (1) the moving average filter window is used as the adjustment window
%         (2) the breakpoint is a peak in the opposite direction
%         (3) the adjustment window cannot go beyond the smoothed midpoint 
%               between the previous peak and breakpoint on the left
%         (4) the adjustment window cannot go beyond the smoothed midpoint 
%               between the breakpoint and the peak on the right
[idxBreaks, valBreaks] = ...
    adjust_peaks(data, idxBreaksFilt, directionFactor * -1, smoothWindow, ...
                 'LeftBounds', idxMidPrevPeaks2Breaks, ...
                 'RightBounds', idxMidBreaks2Peaks);

% For all direction-filtered events, 
%   find the corresponding peak indices and values in the original trace
%   Note: (1) the moving average filter window is used as the adjustment window
%         (2) the adjustment window cannot go beyond the smoothed midpoint 
%               between the breakpoint and the peak on the left
%         (3) the adjustment window cannot go beyond the smoothed midpoint 
%               between the peak and the next breakpoint on the right
[idxPeaks, valPeaks] = ...
    adjust_peaks(data, idxPeaksFilt, directionFactor, smoothWindow, ...
                 'LeftBounds', idxMidBreaks2Peaks, ...
                 'RightBounds', idxMidPeaks2NextBreaks);

% For all direction-filtered events, find the peak amplitudes
peakAmplitudes = abs(valPeaks - valBreaks);

% Find all peaks that satisfies the amplitude threshold
peakNoNotTooSmall = find(peakAmplitudes > eventAmpThreshold);
                                                % peak #s that are not too small

% Make sure the window for computing baseline is positive
%   This is necessary because of potential rounding issues
if baselineWindow <= 0
    baselineWindow = 1;
    fprintf('Warning: Window for computing baseline must be positive\n');
    fprintf('         It''s now changed to 1!\n\n');
end

% Compute baseline value before each breakpoint
baselineStarts = max(1, idxBreaks - baselineWindow + 1);
valBaselines = arrayfun(@(x, y) mean(data(x:y)), baselineStarts, idxBreaks);
% TODO: Try this: compute_means(data, 'Endpoints', transpose(baselineStarts, idxBreaks))
%   from Adams_Functions

% Find all peaks that are not "buried in the baseline"
peakNoNotBuried = find(100 .* (valBaselines - valBreaks) .* directionFactor ...
                            ./ peakAmplitudes <= maxBelowBasePerc);
                                                % peak #s that are not buried

% Find all peaks that are events
peakNoEvent = intersect(peakNoNotTooSmall, peakNoNotBuried);

% If no events are detected, return an empty matrix and exit function
if isempty(peakNoEvent)
    eventInfo = [];
    fprintf('No directional events are detected!!\n')
    return;
end

% Otherwise, get information about the event breakpoints and peaks
idxEventBreaks = idxBreaks(peakNoEvent);        % indices of event breakpoints
idxEventPeaks = idxPeaks(peakNoEvent);          % indices of event peaks
valEventBreaks = valBreaks(peakNoEvent);        % values at breakpoints
valEventPeaks = valPeaks(peakNoEvent);          % values at peaks
eventAmplitudes = peakAmplitudes(peakNoEvent);  % amplitudes of events

% Compute the 0-100% rise times
totalRiseTimes = idxEventPeaks - idxEventBreaks;
                                                % 0-100% rise times (samples)

% Compute the 10-90% rise times
rise10Vals = valEventBreaks + directionFactor * 0.1 * eventAmplitudes;
                                                % 10% of event peak
rise90Vals = valEventBreaks + directionFactor * 0.9 * eventAmplitudes;    
                                                % 90% of event peak
rise10Times = ...
    arrayfun(@(x, y, z) find_custom(data(x:y) * directionFactor < ...
                             z * directionFactor, 1, 'last', ...
                             'ReturnNan', true), ...
            idxEventBreaks, idxEventPeaks, rise10Vals);
rise90Times = ...
    arrayfun(@(x, y, z) find_custom(data(x:y) * directionFactor > ...
                             z * directionFactor, 1, 'first', ...
                             'ReturnNan', true), ...
            idxEventBreaks, idxEventPeaks, rise90Vals);
tenNinetyRiseTimes = rise90Times - rise10Times; % 10-90% rise times (samples)

% Compute the inter-event intervals to the next event or the end of vector
%   Here, "inter-event interval" is defined from peak to peak, 
%   consistent with that used by Molecular Devices
interEventIntervals = [idxEventPeaks(2:end); nSamples] - idxEventPeaks;    
                                            % inter-event intervals (samples)

% Compute the interstimulus intervals to the next event or the end of vector
%   Here, "interstimulus interval" is defined from peak to breakpoint
interStimulusIntervals = [idxEventBreaks(2:end); nSamples] - idxEventPeaks;    
                                            % interstimulus intervals (samples)

% Compute the 50% decay times:
%   The first point where the value returns to 50% of peak amplitude
halfDecayVals = valEventPeaks - directionFactor * 0.5 * eventAmplitudes;
                                            % half decay value (data units)
halfDecayTimes = ...
    arrayfun(@(x, y, z) find_custom(data(x:y) * directionFactor < ...
                             z * directionFactor, 1, 'first', ...
                             'ReturnNan', true), ...
            idxEventPeaks, idxEventPeaks + interStimulusIntervals, ...
            halfDecayVals) - 1;
                                            % half decay times (samples)

% Compute the "full decay times":
%   The first point where the value returns to within noiseLevel 
%       of breakpoint value
fullDecayVals = valEventBreaks + directionFactor * noiseLevel;
                                            % full decay value (data units)
fullDecayTimes = ...
    arrayfun(@(x, y, z) find_custom(data(x:y) * directionFactor < ...
                             z * directionFactor, 1, 'first', ...
                             'ReturnNan', true), ...
            idxEventPeaks, idxEventPeaks + interStimulusIntervals, ...
            fullDecayVals) - 1;
                                            % full decay times (samples)

% Place information about the events into a single matrix
nEvents = length(idxEventPeaks);                % number of events found
eventInfo = zeros(nEvents, 11);
eventInfo(:, IDXBREAK_COLNUM)      = idxEventBreaks;
eventInfo(:, IDXPEAK_COLNUM)       = idxEventPeaks;
eventInfo(:, VALBREAK_COLNUM)      = valEventBreaks;
eventInfo(:, VALPEAK_COLNUM)       = valEventPeaks;
eventInfo(:, EVENTAMP_COLNUM)      = eventAmplitudes;
eventInfo(:, TOTALRISE_COLNUM)     = totalRiseTimes;
eventInfo(:, TENNINETYRISE_COLNUM) = tenNinetyRiseTimes;
eventInfo(:, IEI_COLNUM)           = interEventIntervals;
eventInfo(:, ISI_COLNUM)           = interStimulusIntervals;
eventInfo(:, HALFDECAY_COLNUM)     = halfDecayTimes;
eventInfo(:, FULLDECAY_COLNUM)     = fullDecayTimes;

fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%                           3 = smoothed data value at event breakpoint
%                           4 = smoothed data value at event peak
%                           5 = smoothed amplitude of the event
peakAmplitudes = abs(dataDirFilt(idxPeaks));    % all peak amplitudes
valEventBreaks = dataSmooth(idxEventBreaks);    % smoothed values at breakpoints
valEventPeaks = dataSmooth(idxEventPeaks);      % smoothed values at peaks

halfDecayTimes = ...
    cellfun(@(x, y, z) find_custom(dataSmooth(x:y) * directionFactor < ...
                             z * directionFactor, 1, 'first', ...
                             'ReturnNan', true), ...
            num2cell(idxEventPeaks), ...
            num2cell(idxEventPeaks + interStimulusIntervals), ...
            num2cell(halfDecayVals));

idxBreaks = idxBreaksFilt;                  % use the same breakpoint indices
valBreaks = data(idxBreaksFilt);            % adjust the breakpoint values

[~, rise10Times] = ...
    cellfun(@(x, y, z) min(abs(data(x:y)-z)), ...
            num2cell(idxEventBreaks), num2cell(idxEventPeaks), ...
            num2cell(rise10Vals));
[~, rise90Times] = ...
    cellfun(@(x, y, z) min(abs(data(x:y)-z)), ...
            num2cell(idxEventBreaks), num2cell(idxEventPeaks), ...
            num2cell(rise90Vals));

rise10Vals = 0.1 * data(idxEventPeaks);         % 10% of event peak value
rise90Vals = 0.9 * data(idxEventPeaks);         % 90% of event peak value

rise10Times = ...
    cellfun(@(x, y, z) find_custom(data(x:y) * directionFactor < ...
                             z * directionFactor, 1, 'last', ...
                             'ReturnNan', true), ...
            num2cell(idxEventBreaks), num2cell(idxEventPeaks), ...
            num2cell(rise10Vals)) - 1;
rise90Times = ...
    cellfun(@(x, y, z) find_custom(data(x:y) * directionFactor > ...
                             z * directionFactor, 1, 'first', ...
                             'ReturnNan', true), ...
            num2cell(idxEventBreaks), num2cell(idxEventPeaks), ...
            num2cell(rise90Vals)) - 1;

fullDecayTimes = ...
    cellfun(@(x, y, z) find_custom(data(x:y) * directionFactor < ...
                             z * directionFactor, 1, 'first', ...
                             'ReturnNan', true), ...
            num2cell(idxEventPeaks), ...
            num2cell(idxEventPeaks + interStimulusIntervals), ...
            num2cell(fullDecayVals)) - 1;
                                            % full decay times (samples)
halfDecayTimes = ...
    cellfun(@(x, y, z) find_custom(data(x:y) * directionFactor < ...
                             z * directionFactor, 1, 'first', ...
                             'ReturnNan', true), ...
            num2cell(idxEventPeaks), ...
            num2cell(idxEventPeaks + interStimulusIntervals), ...
            num2cell(halfDecayVals)) - 1;
                                            % half decay times (samples)
rise10Times = ...
    cellfun(@(x, y, z) find_custom(data(x:y) * directionFactor < ...
                             z * directionFactor, 1, 'last', ...
                             'ReturnNan', true), ...
            num2cell(idxEventBreaks), num2cell(idxEventPeaks), ...
            num2cell(rise10Vals));
rise90Times = ...
    cellfun(@(x, y, z) find_custom(data(x:y) * directionFactor > ...
                             z * directionFactor, 1, 'first', ...
                             'ReturnNan', true), ...
            num2cell(idxEventBreaks), num2cell(idxEventPeaks), ...
            num2cell(rise90Vals));
valBaselines = cellfun(@(x, y) mean(data(x:y)), ...
                num2cell(baselineStarts), num2cell(idxBreaks));


%}

