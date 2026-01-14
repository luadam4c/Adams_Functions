function [peakTable, valleyTable, otherResults] = parse_oscillation (vec, varargin)
%% Parses an oscillatory vector to find peaks, valleys, and their amplitudes
% Usage: [peakTable, valleyTable, otherResults] = parse_oscillation (vec, varargin)
% Explanation:
%       This function identifies all significant peaks and valleys (troughs) in a
%       time-series vector. It can optionally filter the vector before analysis.
%       It returns two tables: one for peaks and one for valleys, with detailed metrics.
%
% Example(s):
%       % Create a sample noisy sine wave
%       fs = 1000;
%       t = (0:1/fs:2-1/fs)';
%       vec = sin(2*pi*5*t) + 0.2 * randn(size(t));
%       
%       % Parse the oscillation to get separate peak and valley tables
%       [pT, vT, o] = parse_oscillation(vec, 'TimeVec', t*1000, 'MinPeakProminence', 0.5);
%       
%       % Parse with peak-to-prevalley amplitude definition
%       [pT2, vT2, o2] = parse_oscillation(vec, 'TimeVec', t*1000, 'MinPeakProminence', 0.5, 'AmpMode', 'peak-to-prevalley');
%
%       % Parse the oscillation, applying a bandpass filter internally
%       [pT, vT, o] = parse_oscillation(vec, 'TimeVec', t*1000, 'SamplingIntervalMs', 1, ...
%                                    'FilterCutoffs', [1, 20], 'FilterOrder', 2, ...
%                                    'MinPeakProminence', 0.5);
%
% Outputs:
%       peakTable   - A table for detected peaks with the following columns:
%                       'peakIndex'       - Index of the peak in the vector.
%                       'peakValue'       - Value of the peak.
%                       'preValleyIndex'  - Index of the preceding valley.
%                       'postValleyIndex' - Index of the succeeding valley.
%                       'preValleyValue'  - Value of the preceding valley.
%                       'postValleyValue' - Value of the succeeding valley.
%                       'amplitude'       - Computed amplitude of the peak.
%                       'peakTime'        - (Optional) Time of the peak.
%                       'preValleyTime'   - (Optional) Time of the preceding valley.
%                       'postValleyTime'  - (Optional) Time of the succeeding valley.
%       valleyTable - A table for detected valleys with the following columns:
%                       'valleyIndex'     - Index of the valley in the vector.
%                       'valleyValue'     - Value of the valley.
%                       'prePeakIndex'    - Index of the preceding peak.
%                       'postPeakIndex'   - Index of the succeeding peak.
%                       'prePeakValue'    - Value of the preceding peak.
%                       'postPeakValue'   - Value of the succeeding peak.
%                       'amplitude'       - Computed amplitude of the valley.
%                       'valleyTime'      - (Optional) Time of the valley.
%                       'prePeakTime'     - (Optional) Time of the preceding peak.
%                       'postPeakTime'    - (Optional) Time of the succeeding peak.
%       otherResults- A structure with other results including:
%                       'freqFundamental' - Fundamental frequency (Hz) if computed.
%
% Arguments:
%       vec         - The time-series vector to analyze.
%                   must be a numeric vector
%       varargin    - 'TimeVec': A time vector corresponding to 'vec'
%                   must be a numeric vector of the same length as 'vec'
%                   default == []
%                   - 'TimeUnits': Units for TimeVec
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'ms'  - milliseconds
%                       's'   - seconds
%                   default == 'ms'
%                   - 'SamplingIntervalMs': Sampling interval in ms
%                           Required if filtering is applied.
%                   must be a positive scalar
%                   default == [] computed from time vector that is assumed to be in ms
%                   - 'FundFreqRange': A two-element vector specifying the 
%                                      range [minHz, maxHz] to search for the
%                                      fundamental frequency.
%                   must be a two-element numeric vector or empty
%                   default == [] (no restriction)
%                   - 'FilterCutoffs': Cutoff frequency(ies) for filtering in Hz.
%                                      Can be a 2-element vector for a bandpass
%                                      filter (e.g., [fLow, fHigh]).
%                   must be a numeric vector
%                   default == [] (no filtering or computed from FilterCutoffsRelToFund)
%                   - 'FilterCutoffsRelToFund': Normalized Cutoff frequency(ies) for filtering
%                                      Can be a 2-element vector for a bandpass
%                                      filter (e.g., [fLowRelToFund, fHighRelToFund]).
%                   must be a numeric vector
%                   default == [] (no filtering)
%                   - 'FilterOrder': Order of the Butterworth filter.
%                   must be a positive integer
%                   default == 2
%                   - 'AmpMode': The method for calculating amplitude.
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'peak-to-equilibrium' - Extrenum minus the average of entire vector
%                       'peak-to-avgvalley'   - Extrenum minus the average of surrounding features
%                       'peak-to-prevalley'   - Extrenum minus the preceding feature
%                       'peak-to-postvalley'  - Extrenum minus the succeding feature
%                   default == 'peak-to-equilibrium'
%                   - 'MinPeakDistanceMs': Minimum distance between peaks in ms.
%                   must be a non-negative scalar
%                   default == [] (not provided)
%                   - 'MinPeakDistance': Minimum distance between peaks.
%                   must be empty or a non-negative scalar
%                   default == 0
%                   - 'MinValleyDistanceMs': Minimum distance between valleys in ms.
%                   must be a non-negative scalar
%                   default == [] (not provided)
%                   - 'MinValleyDistance': Minimum distance between valleys.
%                   must be empty or a non-negative scalar
%                   default == 0
%                   - 'PromThresholdPerc': Percentage of full data range for 
%                           the prominence threshold for peak/valley detection. 
%                           Final threshold is maximum of this and 
%                           either MinPeakProminence or MinValleyProminence
%                   must be empty or a non-negative scalar
%                   default == [] (not provided)
%                   - 'MinPeakProminence': Minimum prominence for peak detection.
%                           Final threshold is maximum of this and 
%                           PromThresholdPerc * range(vec)
%                   must be empty or a non-negative scalar
%                   default == 0
%                   - 'MinValleyProminence': Minimum prominence for valley detection.
%                           Final threshold is maximum of this and 
%                           PromThresholdPerc * range(vec)
%                   must be empty or a non-negative scalar
%                   default == 0
%                   - 'MaxDurationForAmplitudeMs': Maximum difference between surrounding features for amplitude calculation
%                   must be empty or a non-negative scalar
%                   default == [] (none)
%                   - Any other parameter-value pair for findpeaks()
%
% Requires:
%       cd/compute_fundamental_frequency.m
%       cd/isnumericvector.m
%       cd/freqfilter.m
%
% Used by:
%       cd/virt_analyze_sniff_whisk.m
%       \Shared\Code\vIRt\virt_moore.m

% File History:
% 2025-09-04 Created by Gemini, based on logic from virt_moore.m
% 2025-09-04 Modified by Gemini to output two separate tables.
% 2025-09-04 Modified by Gemini to include internal filtering capabilities.
% 2025-09-05 Added 'MinPeakDistanceMs', 'MinValleyDistanceMs' and 
%               'PromThresholdPerc', 'FilterCutoffsRelToFund' as optional arguments
% 2025-09-05 Added 'TimeUnits' as an optional argument by Gemini.
% 2025-09-05 Added 'FundFreqRange' as an optional argument
% 2025-09-09 Fixed pre- and post-valley detection
% 2025-09-12 'PromThresholdPerc' and 'MinPeakProminence' are no longer
%               exclusive
% 2025-09-13 Added 'MaxDurationForAmplitudeMs' as an optional argument
% TODO: Make 'ParsePsd' an optional argument and use parse_psd.m and 
%       store results in otherResults

%% Hard-coded parameters
validAmpModes = {'peak-to-equilibrium', 'peak-to-avgvalley', 'peak-to-prevalley', 'peak-to-postvalley'};
validTimeUnits = {'ms', 's'};

%% Default values for optional arguments
timeVecDefault = [];
timeUnitsDefault = 'ms';
siMsDefault = [];
fundFreqRangeDefault = [];          % no frequency range restriction by default
filterCutoffsDefault = [];
filterCutoffsRelToFundDefault = [];
filterOrderDefault = 2;
ampModeDefault = 'peak-to-equilibrium';
minPeakDistanceMsDefault = [];
minPeakDistanceDefault = [];
minValleyDistanceMsDefault = [];
minValleyDistanceDefault = [];
promThresholdPercDefault = [];
minPeakProminenceDefault = [];
minValleyProminenceDefault = [];
maxDurationForAmplitudeMsDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error('A vector must be provided.');
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'vec', @isnumericvector);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TimeVec', timeVecDefault, @isnumericvector);
addParameter(iP, 'TimeUnits', timeUnitsDefault, ...
    @(x) any(validatestring(x, validTimeUnits)));
addParameter(iP, 'SamplingIntervalMs', siMsDefault, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'FundFreqRange', fundFreqRangeDefault, ...         % The frequency range in Hz
    @(x) isempty(x) || (isnumeric(x) && isvector(x) && numel(x) == 2)); % Must be a 2-element vector or empty
addParameter(iP, 'FilterCutoffs', filterCutoffsDefault, @isnumeric);
addParameter(iP, 'FilterCutoffsRelToFund', filterCutoffsRelToFundDefault, @isnumeric);
addParameter(iP, 'FilterOrder', filterOrderDefault, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'AmpMode', ampModeDefault, ...
    @(x) any(validatestring(x, validAmpModes)));
addParameter(iP, 'MinPeakDistanceMs', minPeakDistanceMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'MinPeakDistance', minPeakDistanceDefault, @isnumeric);
addParameter(iP, 'MinValleyDistanceMs', minValleyDistanceMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'MinValleyDistance', minValleyDistanceDefault, @isnumeric);
addParameter(iP, 'PromThresholdPerc', promThresholdPercDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'MinPeakProminence', minPeakProminenceDefault, @isnumeric);
addParameter(iP, 'MinValleyProminence', minValleyProminenceDefault, @isnumeric);
addParameter(iP, 'MaxDurationForAmplitudeMs', maxDurationForAmplitudeMsDefault, @isnumeric);

% Read from the Input Parser
parse(iP, vec, varargin{:});
timeVec = iP.Results.TimeVec;
timeUnits = validatestring(iP.Results.TimeUnits, validTimeUnits);
siMs = iP.Results.SamplingIntervalMs;
fundFreqRange = iP.Results.FundFreqRange;                           % Store the frequency range
filterCutoffs = iP.Results.FilterCutoffs;
filterCutoffsRelToFund = iP.Results.FilterCutoffsRelToFund;
filterOrder = iP.Results.FilterOrder;
ampMode = validatestring(iP.Results.AmpMode, validAmpModes);
minPeakDistanceMs = iP.Results.MinPeakDistanceMs;
minPeakDistance = iP.Results.MinPeakDistance;
minValleyDistanceMs = iP.Results.MinValleyDistanceMs;
minValleyDistance = iP.Results.MinValleyDistance;
promThresholdPerc = iP.Results.PromThresholdPerc;
minPeakProminence = iP.Results.MinPeakProminence;
minValleyProminence = iP.Results.MinValleyProminence;
maxDurationForAmplitudeMs = iP.Results.MaxDurationForAmplitudeMs;

% Keep unmatched arguments for the findpeaks() function
otherArguments = iP.Unmatched;

% Transform to a cell array
otherArguments = struct2arglist(otherArguments);

% Validate TimeVec length if provided
if ~isempty(timeVec) && numel(timeVec) ~= numel(vec)
    error('TimeVec must have the same number of elements as vec.');
end

% Decide whether filtering is requested
if ~isempty(filterCutoffs) || ~isempty(filterCutoffsRelToFund)
    toFilter = true;
    if ~isempty(filterCutoffs) && ~isempty(filterCutoffsRelToFund)
        warning('filterCutoffsRelToFund will be overridden by filterCutoffs!');
    end
else
    toFilter = false;
end    

% Initialize other results structure
otherResults = struct;

%% Decide on default values
% Decide on values for filtering
if toFilter 
    % Decide on the sampling interval
    if isempty(siMs)
        if ~isempty(timeVec)
            % Get the mean sampling interval
            avgSI = mean(diff(timeVec));

            % Convert to ms if time vector is in seconds
            switch timeUnits
            case 's'
                siMs = avgSI * 1000;
            case 'ms'
                siMs = avgSI;
            end
         else
            error('Either Time vector (TimeVec) or Sampling Interval (SamplingIntervalMs) must be provided for filtering.');
        end
    end

    % Compute the fundamental frequency of the oscillation (Hz)
    freqFundamental = ...
        compute_fundamental_frequency(vec, 'SamplingIntervalMs', siMs, ...
                                    'FundFreqRange', fundFreqRange);

    % Define filter cutoffs if fundamental frequency is found
    if isempty(filterCutoffs) && ~isempty(filterCutoffsRelToFund)
        filterCutoffs = filterCutoffsRelToFund * freqFundamental;
    end

    % Save fundamental frequency in otherResults
    otherResults.freqFundamental = freqFundamental;
end


% Decide on minimum peak/valley distance
if isempty(minPeakDistance) && ~isempty(minPeakDistanceMs) && ~isempty(siMs)
    minPeakDistance = round(minPeakDistanceMs / siMs);
else
    minPeakDistance = 0;
end
if isempty(minValleyDistance) && ~isempty(minValleyDistanceMs) && ~isempty(siMs)
    minValleyDistance = round(minValleyDistanceMs / siMs);
else
    minValleyDistance = 0;
end

if isempty(maxDurationForAmplitudeMs)
    maxDurationForAmplitude = [];
else
    maxDurationForAmplitude = round(maxDurationForAmplitudeMs / siMs);
end

% Decide on prominence thresholds
if ~isempty(minPeakProminence) && ~isempty(promThresholdPerc)
    minPeakProminence = max([minPeakProminence, ...
                            (promThresholdPerc / 100) * range(vec)]);
elseif isempty(minPeakProminence) && ~isempty(promThresholdPerc)
    minPeakProminence = (promThresholdPerc / 100) * range(vec);
elseif isempty(minPeakProminence) && isempty(promThresholdPerc)
    minPeakProminence = 0;
else
    % Do nothing
end
if ~isempty(minValleyProminence) && ~isempty(promThresholdPerc)
    minValleyProminence = max([minValleyProminence, ...
                            (promThresholdPerc / 100) * range(vec)]);
elseif isempty(minValleyProminence) && ~isempty(promThresholdPerc)
    minValleyProminence = (promThresholdPerc / 100) * range(vec);
elseif isempty(minValleyProminence) && isempty(promThresholdPerc)
    minValleyProminence = 0;
else
    % Do nothing
end

%% Filter the vector if requested
if ~isempty(filterCutoffs)
    % Apply a Butterworth bandpass filter
    vecFiltered = freqfilter(vec, filterCutoffs, siMs / 1000, ...
                                'FilterType', 'bandpass', ...
                                'FilterOrder', filterOrder);
else
    % Use the original vector if no filtering is specified
    vecFiltered = vec;
end

%% Find all peaks and valleys on the (potentially filtered) vector
% Detect all peaks
[~, peakIndices] = findpeaks(vecFiltered, ...
                            'MinPeakProminence', minPeakProminence, ...
                            'MinPeakDistance',  minPeakDistance, otherArguments{:});


% Detect all valleys (troughs) by inverting the vector
[~, valleyIndices] = findpeaks(-vecFiltered, ...
                            'MinPeakProminence', minValleyProminence, ...
                            'MinPeakDistance', minValleyDistance, otherArguments{:});

% Get values from the ORIGINAL vector
peakValues = vec(peakIndices);
valleyValues = vec(valleyIndices);

%% Create Peak Table
nPeaks = numel(peakIndices);

% Initialize arrays for the peak table
preValleyIndices = nan(nPeaks, 1);
postValleyIndices = nan(nPeaks, 1);
preValleyValues = nan(nPeaks, 1);
postValleyValues = nan(nPeaks, 1);
peakAmplitudes = nan(nPeaks, 1);

% For each peak, find its surrounding valleys and compute amplitude if possible
for i = 1:nPeaks
    idxCurrentPeak = peakIndices(i);

    % Find the closest preceding valley
    precedingValleys = valleyIndices(valleyIndices < idxCurrentPeak);
    if ~isempty(precedingValleys)
        idxPreValley = precedingValleys(end);
        valPreValley = vec(idxPreValley);
        preValleyIndices(i) = idxPreValley;
        preValleyValues(i) = valPreValley;
    end

    % Find the closest succeeding valley
    succeedingValleys = valleyIndices(valleyIndices > idxCurrentPeak);
    if ~isempty(succeedingValleys)
        idxPostValley = succeedingValleys(1);
        valPostValley = vec(idxPostValley);
        postValleyIndices(i) = idxPostValley;
        postValleyValues(i) = valPostValley;
    end

    % Compute the amplitude for the peak based on the selected mode
    %   only if preceding and succeeding valleys are present
    if ~isempty(precedingValleys) && ~isempty(succeedingValleys)
        valPeak = vec(idxCurrentPeak);
        switch ampMode
            case 'peak-to-equilibrium'
                peakAmplitudes(i) = valPeak - mean(vec);
            case 'peak-to-avgvalley'
                peakAmplitudes(i) = valPeak - mean([valPreValley, valPostValley]);
            case 'peak-to-prevalley'
                peakAmplitudes(i) = valPeak - valPreValley;
            case 'peak-to-postvalley'
                peakAmplitudes(i) = valPeak - valPostValley;
        end

        % If amplitude is negative, the data is too noisy so do not 
        %   count this peak
        if peakAmplitudes(i) < 0
            peakAmplitudes(i) = NaN;
        end

        % If maxDurationForAmplitudeMs provided,
        %   do not count this peak if difference between 
        %   preceding and succeeeding valley is greater than maxDurationForAmplitudeMs
        if ~isempty(maxDurationForAmplitude) && ...
                (idxPostValley - idxPreValley) > maxDurationForAmplitude
            peakAmplitudes(i) = NaN;
        end
    end
end

% Create the basic peak table
peakTable = table(peakIndices(:), peakValues(:), ...
                  preValleyIndices, postValleyIndices, ...
                  preValleyValues, postValleyValues, peakAmplitudes, ...
                  'VariableNames', {'peakIndex', 'peakValue', ...
                  'preValleyIndex', 'postValleyIndex', ...
                  'preValleyValue', 'postValleyValue', 'amplitude'});

%% Create Valley Table
nValleys = numel(valleyIndices);

% Initialize arrays for the valley table
prePeakIndices = nan(nValleys, 1);
postPeakIndices = nan(nValleys, 1);
prePeakValues = nan(nValleys, 1);
postPeakValues = nan(nValleys, 1);
valleyAmplitudes = nan(nValleys, 1);

% For each valley, find its surrounding peaks and compute amplitude if possible
for i = 1:nValleys
    idxCurrentValley = valleyIndices(i);

    % Find the closest preceding peak
    precedingPeaks = peakIndices(peakIndices < idxCurrentValley);
    if ~isempty(precedingPeaks)
        idxPrePeak = precedingPeaks(end);
        valPrePeak = vec(idxPrePeak);
        prePeakIndices(i) = idxPrePeak;
        prePeakValues(i) = valPrePeak;
    end

    % Find the closest succeeding peak
    succeedingPeaks = peakIndices(peakIndices > idxCurrentValley);
    if ~isempty(succeedingPeaks)
        idxPostPeak = succeedingPeaks(1);
        valPostPeak = vec(idxPostPeak);
        postPeakIndices(i) = idxPostPeak;
        postPeakValues(i) = valPostPeak;
    end

    % Compute the amplitude for the valley based on the selected mode
    %   only if preceding and succeeding peaks are present
    if ~isempty(precedingPeaks) && ~isempty(succeedingPeaks)
        valValley = vec(idxCurrentValley);
        switch ampMode
            case 'peak-to-equilibrium'
                valleyAmplitudes(i) = mean(vec) - valValley;
            case 'peak-to-avgvalley'
                valleyAmplitudes(i) = mean([valPrePeak, valPostPeak]) - valValley;
            case 'peak-to-prevalley'
                valleyAmplitudes(i) = valPrePeak - valValley;
            case 'peak-to-postvalley'
                valleyAmplitudes(i) = valPostPeak - valValley;
        end
    end
end

% Create the basic valley table
valleyTable = table(valleyIndices(:), valleyValues(:), ...
                    prePeakIndices, postPeakIndices, ...
                    prePeakValues, postPeakValues, valleyAmplitudes, ...
                    'VariableNames', {'valleyIndex', 'valleyValue', ...
                    'prePeakIndex', 'postPeakIndex', ...
                    'prePeakValue', 'postPeakValue', 'amplitude'});

%% Add Time Information if TimeVec is provided
if ~isempty(timeVec)
    % Add time columns to peakTable
    peakTimes = timeVec(peakIndices);
    preValleyTimes = nan(nPeaks, 1);
    postValleyTimes = nan(nPeaks, 1);
    validPreV = ~isnan(preValleyIndices);
    validPostV = ~isnan(postValleyIndices);
    preValleyTimes(validPreV) = timeVec(preValleyIndices(validPreV));
    postValleyTimes(validPostV) = timeVec(postValleyIndices(validPostV));
    peakTable = addvars(peakTable, peakTimes, preValleyTimes, postValleyTimes, ...
                      'NewVariableNames', {'peakTime', 'preValleyTime', 'postValleyTime'});

    % Add time columns to valleyTable
    valleyTimes = timeVec(valleyIndices);
    prePeakTimes = nan(nValleys, 1);
    postPeakTimes = nan(nValleys, 1);
    validPreP = ~isnan(prePeakIndices);
    validPostP = ~isnan(postPeakIndices);
    prePeakTimes(validPreP) = timeVec(prePeakIndices(validPreP));
    postPeakTimes(validPostP) = timeVec(postPeakIndices(validPostP));
    valleyTable = addvars(valleyTable, valleyTimes, prePeakTimes, postPeakTimes, ...
                        'NewVariableNames', {'valleyTime', 'prePeakTime', 'postPeakTime'});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
