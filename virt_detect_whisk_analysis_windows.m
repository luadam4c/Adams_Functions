function [sniffStartWinTable, basalRespCycleTable] = virt_detect_whisk_analysis_windows (peakTable, valleyTable, analysisParams, varargin)
%% Detects whisk analysis windows for sniffing and basal respiration
% Usage: [sniffStartWinTable, basalRespCycleTable] = virt_detect_whisk_analysis_windows (peakTable, valleyTable, analysisParams, varargin)
%
% Explanation:
%       This function identifies time windows corresponding to "sniff starts" and
%       "basal respiration cycles" based on provided whisk and respiratory data.
%       It extracts the relevant whisk peaks, valleys, amplitudes, and computes
%       logarithmic decrements for each identified window. The function can
%       operate in different modes to accommodate data from simulations (driven
%       by pulse cycles) or experimental recordings (driven by detected
%       respiratory events). It also computes phase-response data when in
%       'pulseCycle' or 'respPeak' mode for basal respiration.
%
% Outputs:
%       sniffStartWinTable  - A table where each row represents a detected
%                             sniff-start window and contains associated whisk data.
%                           specified as a table
%       basalRespCycleTable - A table where each row represents a detected
%                             basal respiration cycle and contains associated whisk data.
%                           specified as a table
%
% Arguments:
%       peakTable           - A table of detected whisk peaks from parse_oscillation.m.
%                           must be a table
%       valleyTable         - A table of detected whisk valleys from parse_oscillation.m.
%                           must be a table
%       analysisParams      - A structure containing analysis parameters, such as:
%                               nWhisksSniffStartToAnalyze, sniffFreqThreshold,
%                               minWhisksBasalRespToAnalyze, etc.
%                           must be a structure
%       varargin            - 'SniffDefinitionMode': Method to define sniff periods.
%                           must be one of 'pulseCycle', 'transitionTimes'
%                           default == 'pulseCycle'
%                           - 'BasalDefinitionMode': Method to define basal respiration cycles.
%                           must be one of 'pulseCycle', 'respPeak'
%                           default == 'pulseCycle'
%                           - 'PulseCycleStartTimes': Vector of start times for each pulse cycle.
%                           must be a numeric vector
%                           default == []
%                           - 'PulseCycleEndTimes': Vector of end times for each pulse cycle.
%                           must be a numeric vector
%                           default == []
%                           - 'PulseNumbers': Vector of identifiers for each pulse cycle.
%                           must be a numeric vector
%                           default == []
%                           - 'SeedNumber': The seed number for the simulation.
%                           must be a numeric scalar
%                           default == 0
%                           - 'TrialName': The trial name for the experiment.
%                           must be a string scalar or a character vector
%                           default == ''
%                           - 'TVecMs': The full time vector of the recording in ms (for 'pulseCycle').
%                           must be a numeric vector
%                           default == []
%                           - 'TVecs': Time vectors in s (matrix, one col per sweep, for experimental data).
%                           must be a numeric matrix
%                           default == []
%                           - 'SniffStartTimes': Start times of sniff periods.
%                           must be a numeric vector (for 'pulseCycle') or cell array of vectors (for 'transitionTimes')
%                           default == []
%                           - 'SniffEndTimes': End times of sniff periods.
%                           must be a numeric vector (for 'pulseCycle') or cell array of vectors (for 'transitionTimes')
%                           default == []
%                           - 'NSweeps': Number of sweeps (for experimental data).
%                           must be a numeric scalar
%                           default == 0
%                           - 'WhiskPeakTables': Cell array of whisk peak tables per sweep.
%                           must be a cell array of tables
%                           default == {}
%                           - 'WhiskValleyTables': Cell array of whisk valley tables per sweep.
%                           must be a cell array of tables
%                           default == {}
%                           - 'NPeaksByPulse': Vector of whisk peak counts per pulse cycle.
%                           must be a numeric vector
%                           default == []
%                           - 'PeakAmplitudesByPulse': Cell array of whisk amplitudes per pulse cycle.
%                           must be a cell array of numeric vectors
%                           default == {}
%                           - 'PeakTablesByPulse': Cell array of whisk peak tables per pulse cycle.
%                           must be a cell array of tables
%                           default == {}
%                           - 'BasalRespPeakTable': A table with info about basal respiration peaks.
%                           must be a table
%                           default == table()
%                           - 'FileNumber': The file number for logging purposes.
%                           must be a numeric scalar
%                           default == 0
%                           - 'FileIDDiffs': File handle for logging differences.
%                           must be a numeric scalar
%                           default == 1 (stdout)
%                           - 'ToPostMessage': Whether to display messages.
%                           must be a logical scalar
%                           default == true
%                           - 'MessageHandler': Function handle for displaying messages to a GUI.
%                           must be a function handle
%                           default == []
%
% Requires:
%       cd/compute_phase_response.m
%
% Used by:
%       cd/virt_analyze_sniff_whisk.m
%       \Shared\Code\vIRt-Moore\virt_analyze_whisk.m
%

% File History:
% 2025-10-06 Created by Gemini to modularize window detection logic.
% 2025-10-06 Fixed by Gemini to include messaging, phase response, and annotations.
% 2025-10-06 Updated by Gemini to support experimental data from virt_analyze_sniff_whisk.m
% 2025-10-09 Modified by Gemini to add seedNumber/fileNumber as first column
%               and create mode-specific empty tables to remove unused columns.
% 2025-10-09 Modified by Gemini to also pass in and add trialName for experiments.

%% Hard-coded parameters
validSniffDefinitions = {'pulseCycle', 'transitionTimes'};
validBasalDefinitions = {'pulseCycle', 'respPeak'};

%% Default values for optional arguments
sniffDefinitionModeDefault = 'pulseCycle';
basalDefinitionModeDefault = 'pulseCycle';
pulseCycleStartTimesDefault = [];
pulseCycleEndTimesDefault = [];
pulseNumbersDefault = [];
seedNumberDefault = 0;
trialNameDefault = '';
tVecMsDefault = [];
tVecsDefault = [];
sniffStartTimesDefault = [];
sniffEndTimesDefault = [];
nSweepsDefault = 0;
whiskPeakTablesDefault = {};
whiskValleyTablesDefault = {};
nPeaksByPulseDefault = [];
peakAmplitudesByPulseDefault = {};
peakTablesByPulseDefault = {};
basalRespPeakTableDefault = table();
fileNumberDefault = 0;
fileIDDiffsDefault = 1; % stdout
toPostMessageDefault = true;
messageHandlerDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions
% Function to post messages to command window and/or GUI
function post_message(msg)
    if toPostMessage
        % Print to command window
        fprintf('%s\n', msg);

        % Also send to GUI if handler exists
        if ~isempty(messageHandler)
            messageHandler(msg);
        end
    end
end

function isValid = is_first_n_peaks_valid(peakAmplitudes, N)
%% Checks if the first N elements of a vector are not NaN
    % Compute the maximum index to check
    maxIdxToCheck = min(N, numel(peakAmplitudes));

    % Return whether all of them are not NaNs
    isValid = all(~isnan(peakAmplitudes(1:maxIdxToCheck)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required arguments
addRequired(iP, 'peakTable', @istable);
addRequired(iP, 'valleyTable', @istable);
addRequired(iP, 'analysisParams', @isstruct);

% Add parameter-value pairs
addParameter(iP, 'SniffDefinitionMode', sniffDefinitionModeDefault, ...
    @(x) any(validatestring(x, validSniffDefinitions)));
addParameter(iP, 'BasalDefinitionMode', basalDefinitionModeDefault, ...
    @(x) any(validatestring(x, validBasalDefinitions)));
addParameter(iP, 'PulseCycleStartTimes', pulseCycleStartTimesDefault, @isnumeric);
addParameter(iP, 'PulseCycleEndTimes', pulseCycleEndTimesDefault, @isnumeric);
addParameter(iP, 'PulseNumbers', pulseNumbersDefault, @isnumeric);
addParameter(iP, 'SeedNumber', seedNumberDefault, @isnumeric);
addParameter(iP, 'TrialName', trialNameDefault, @(x) ischar(x) || isstring(x));
addParameter(iP, 'TVecMs', tVecMsDefault, @isnumeric);
addParameter(iP, 'TVecs', tVecsDefault, @isnumeric);
addParameter(iP, 'SniffStartTimes', sniffStartTimesDefault);
addParameter(iP, 'SniffEndTimes', sniffEndTimesDefault);
addParameter(iP, 'NSweeps', nSweepsDefault, @isnumeric);
addParameter(iP, 'WhiskPeakTables', whiskPeakTablesDefault, @iscell);
addParameter(iP, 'WhiskValleyTables', whiskValleyTablesDefault, @iscell);
addParameter(iP, 'NPeaksByPulse', nPeaksByPulseDefault, @isnumeric);
addParameter(iP, 'PeakAmplitudesByPulse', peakAmplitudesByPulseDefault, @iscell);
addParameter(iP, 'PeakTablesByPulse', peakTablesByPulseDefault, @iscell);
addParameter(iP, 'BasalRespPeakTable', basalRespPeakTableDefault, @istable);
addParameter(iP, 'FileNumber', fileNumberDefault, @isnumeric);
addParameter(iP, 'FileIDDiffs', fileIDDiffsDefault);
addParameter(iP, 'ToPostMessage', toPostMessageDefault, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(iP, 'MessageHandler', messageHandlerDefault);

% Parse the inputs
parse(iP, peakTable, valleyTable, analysisParams, varargin{:});
sniffDefinitionMode = iP.Results.SniffDefinitionMode;
basalDefinitionMode = iP.Results.BasalDefinitionMode;
pulseCycleStartTimes = iP.Results.PulseCycleStartTimes;
pulseCycleEndTimes = iP.Results.PulseCycleEndTimes;
pulseNumbers = iP.Results.PulseNumbers;
seedNumber = iP.Results.SeedNumber;
trialName = iP.Results.TrialName;
tVecMs = iP.Results.TVecMs;
tVecs = iP.Results.TVecs;
sniffStartTimes = iP.Results.SniffStartTimes;
sniffEndTimes = iP.Results.SniffEndTimes;
nSweeps = iP.Results.NSweeps;
whiskPeakTables = iP.Results.WhiskPeakTables;
whiskValleyTables = iP.Results.WhiskValleyTables;
nPeaksByPulse = iP.Results.NPeaksByPulse;
peakAmplitudesByPulse = iP.Results.PeakAmplitudesByPulse;
peakTablesByPulse = iP.Results.PeakTablesByPulse;
basalRespPeakTable = iP.Results.BasalRespPeakTable;
fileNumber = iP.Results.FileNumber;
fileIDDiffs = iP.Results.FileIDDiffs;
toPostMessage = iP.Results.ToPostMessage;
messageHandler = iP.Results.MessageHandler;

%% Preparation
% Extract key parameters from the analysisParams struct for easier access
nWhisksSniffStartToAnalyze = analysisParams.nWhisksSniffStartToAnalyze;
sniffFreqThreshold = analysisParams.sniffFreqThreshold;
basalFreqThreshold = analysisParams.basalFreqThreshold;
minWhisksBasalRespToAnalyze = analysisParams.minWhisksBasalRespToAnalyze;
maxWhisksBasalRespToAnalyze = analysisParams.maxWhisksBasalRespToAnalyze;
whiskDirForPhase = analysisParams.whiskDirForPhase;
breathOnsetLatencyMs = analysisParams.breathOnsetLatencyMs;

% Extract all peak and valley times for quick lookups
if ~isempty(peakTable)
    allPeakTimes = peakTable.peakTime;
end

if ~isempty(valleyTable)
    allValleyTimes = valleyTable.valleyTime;
    allValleyValues = valleyTable.valleyValue;
end

%% --- SNIFF START WINDOW DETECTION ---
% Starting message
post_message('Starting Sniff Start Window Detection ...');

switch sniffDefinitionMode
case 'pulseCycle'
    % Define empty table structure for consistent output (Simulation)
    emptySniffTable = table('Size', [0, 15], 'VariableTypes', ...
        {'double', 'double', 'double', 'double', 'double', 'double', 'double', ...
         'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'}, ...
        'VariableNames', {'seedNumber', 'pulseNumber', 'windowNumber', 'sniffStartWinStartTime', ...
                          'sniffStartWinEndTime', 'sniffStartTime', 'sniffEndTime', ...
                          'whiskPeakTimes', 'whiskPeakValues', 'whiskPeakAmplitudes', ...
                          'whiskPreValleyTimes', 'whiskPostValleyTimes', ...
                          'whiskValleyTimes', 'whiskValleyValues', 'whiskLogDecrements'});

    % For simulation data: define sniff periods based on pulse cycle durations
    
    % Compute the sniff cycle duration threshold in ms
    sniffDurationThresholdMs = 1000 / sniffFreqThreshold;

    % Calculate pulse cycle durations
    cycleDurations = pulseCycleEndTimes - pulseCycleStartTimes;

    % Identify sniff cycles (duration shorter than threshold)
    isSniffCycle = cycleDurations < sniffDurationThresholdMs;

    % Find transitions to identify the start and end of sniffing periods
    sniffStartIndices = find(diff([false; isSniffCycle]) == 1);
    sniffEndIndices = find(diff([isSniffCycle; false]) == -1);

    % Get the times for these start and end events
    sniffStartTimes = pulseCycleStartTimes(sniffStartIndices);
    sniffEndTimes = pulseCycleEndTimes(sniffEndIndices);

    % Pair up start and end times correctly, handling edge cases
    % Case 1: An end time occurs before the first start time (e.g., starts mid-sniff)
    if ~isempty(sniffEndTimes) && ~isempty(sniffStartTimes) && sniffEndTimes(1) < sniffStartTimes(1)
        sniffEndTimes(1) = [];
    end
    % Case 2: More start times than end times (e.g., ends mid-sniff)
    if numel(sniffStartTimes) > numel(sniffEndTimes)
        sniffEndTimes(end+1) = tVecMs(end);
    end
    pulseNumbersForSniff = pulseNumbers(sniffStartIndices);

    % Count the number of complete sniff periods
    nSniffPeriods = numel(sniffStartTimes);
    
    % Loop through each complete sniff period
    sniffWinRows = {};
    iSniffStartWin = 0;
    for iSniffPeriod = 1:nSniffPeriods
        currentSniffStart = sniffStartTimes(iSniffPeriod);
        currentSniffEnd = sniffEndTimes(iSniffPeriod);
        
        % Find all peaks that occurred within this entire sniff period
        peaksInPeriod = peakTable(peakTable.peakTime >= currentSniffStart & ...
                                  peakTable.peakTime <= currentSniffEnd, :);

        % Check if there are enough whisks to form a valid window
        if height(peaksInPeriod) >= nWhisksSniffStartToAnalyze
            % Take the first N whisks from the start of the period
            peaksToAnalyze = peaksInPeriod(1:nWhisksSniffStartToAnalyze, :);
            
            % Check if all required peaks have valid amplitudes
            if all(~isnan(peaksToAnalyze.amplitude))
                % Increment window count
                iSniffStartWin = iSniffStartWin + 1;

                % Define window start as the pre-valley of the first whisk
                winStart = peaksToAnalyze.preValleyTime(1);
                if isnan(winStart), winStart = currentSniffStart; end % Fallback

                % Define window end as the post-valley of the last whisk
                winEnd = peaksToAnalyze.postValleyTime(end);
                if isnan(winEnd), winEnd = currentSniffEnd; end % Fallback
                
                % Find all valleys within the sniff start window
                isValleyInWin = allValleyTimes >= winStart & allValleyTimes <= winEnd;
                valleysToAnalyze = valleyTable(isValleyInWin, :);
                
                % Compute the logarithmic decrements of successive whisk peak amplitudes
                peakAmplitudes = peaksToAnalyze.amplitude;
                if numel(peakAmplitudes) > 1
                    logDecrements = log(peakAmplitudes(2:end) ./ peakAmplitudes(1:end-1));
                else
                    logDecrements = nan(nWhisksSniffStartToAnalyze - 1, 1);
                end
                
                % Package the data for this window into a one-row table
                newWindow = table(seedNumber, pulseNumbersForSniff(iSniffPeriod), iSniffStartWin, winStart, winEnd, ...
                            currentSniffStart, currentSniffEnd, ...
                            {peaksToAnalyze.peakTime}, {peaksToAnalyze.peakValue}, {peakAmplitudes}, ...
                            {peaksToAnalyze.preValleyTime}, {peaksToAnalyze.postValleyTime}, ...
                            {valleysToAnalyze.valleyTime}, {valleysToAnalyze.valleyValue}, {logDecrements}, ...
                            'VariableNames', emptySniffTable.Properties.VariableNames);

                % Add to sniff start window rows
                sniffWinRows{end + 1} = newWindow;
            end
        end
    end

case 'transitionTimes'
    % Define empty table structure for consistent output (Experiment)
    emptySniffTable = table('Size', [0, 16], 'VariableTypes', ...
        {'double', 'cellstr', 'double', 'double', 'double', 'double', 'double', 'double', ...
         'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'}, ...
        'VariableNames', {'fileNumber', 'trialName', 'sweepNumber', 'windowNumber', 'sniffStartWinStartTime', ...
                          'sniffStartWinEndTime', 'sniffStartTime', 'sniffEndTime', ...
                          'whiskPeakTimes', 'whiskPeakValues', 'whiskPeakAmplitudes', ...
                          'whiskPreValleyTimes', 'whiskPostValleyTimes', ...
                          'whiskValleyTimes', 'whiskValleyValues', 'whiskLogDecrements'});

    % For experimental data: use pre-detected sniff start/end times per sweep
    sniffWinRows = {};
    breathOnsetLatency = breathOnsetLatencyMs / 1000; % Convert to seconds

    for iSwp = 1:nSweeps
        % Obtain data for this sweep
        whiskPeakTable = whiskPeakTables{iSwp};
        whiskValleyTable = whiskValleyTables{iSwp};
        peakTimes = whiskPeakTable.peakTime;
        preValleyTimes = whiskPeakTable.preValleyTime;
        sniffStartTimesThisSweep = sniffStartTimes{iSwp};
        sniffEndTimesThisSweep = sniffEndTimes{iSwp};
        valleyTimes = whiskValleyTable.valleyTime;

        % If there are no sniff periods or no whisks, skip this sweep
        if isempty(sniffStartTimesThisSweep) || isempty(whiskPeakTable)
            continue;
        end
        
        % If the first sniff period end time is before
        %   the first sniff period start time, remove that first end time
        if ~isempty(sniffEndTimesThisSweep) && sniffEndTimesThisSweep(1) < sniffStartTimesThisSweep(1)
            sniffEndTimesThisSweep(1) = [];
        end

        % If there are more sniff period start times then end times, 
        %   add the end of time vector as an sniff end time
        if numel(sniffStartTimesThisSweep) > numel(sniffEndTimesThisSweep)
            sniffEndTimesThisSweep(end+1) = tVecs(end, iSwp);
        end

        % Count the number of complete sniff periods
        nSniffPeriods = numel(sniffStartTimesThisSweep);

        % Loop through all sniff periods for this sweep
        iSniffStartWin = 0; % Counter for sniff start windows for this sweep
        for iSniffPeriod = 1:nSniffPeriods
            currentSniffStart = sniffStartTimesThisSweep(iSniffPeriod);
            currentSniffEnd = sniffEndTimesThisSweep(iSniffPeriod);
            
            % Find all peaks within the sniff period
            isPeakInSniffPeriod = peakTimes >= currentSniffStart & peakTimes <= currentSniffEnd;
            
            % Moore et al 2013 definition: Find all peaks with prevalley times 
            %   after breathOnsetLatency before sniff start time
            %   but peaks before sniff end time
            isPreValleyInSniffPeriod = preValleyTimes >= currentSniffStart - breathOnsetLatency & peakTimes <= currentSniffEnd;

            % If these definitions do not match, print warning
            peakNumInSniffPeriod = find(isPeakInSniffPeriod);
            preValleyNumInSniffPeriod = find(isPreValleyInSniffPeriod);
            preValleyTooEarly = setdiff(peakNumInSniffPeriod, preValleyNumInSniffPeriod);
            if ~isempty(preValleyTooEarly)
                fprintf(fileIDDiffs, ['Whisk protraction for first whisk peak after sniff start occurred ', ...
                         'earlier than %g ms for sniffing period #%d in sweep %d of file %d ', ...
                         ' for peak numbers (for the sweep): %s!!\n\n'], ...
                         breathOnsetLatencyMs, iSniffPeriod, iSwp, fileNumber, num2str(preValleyTooEarly'));
            end

            % Define whether a whisk is in a sniff period by the peak
            whiskPeaksInSniffPeriod = whiskPeakTable(isPeakInSniffPeriod, :);

            % If there are at least nWhisksSniffStartToAnalyze peaks
            %   within this sniff period and the first of those have all valid amplitudes, 
            %   add to sniff start window with appropriate boundaries
            if height(whiskPeaksInSniffPeriod) >= nWhisksSniffStartToAnalyze
                % Extract whisk peaks to analyze
                whiskPeaksToAnalyze = whiskPeaksInSniffPeriod(1:nWhisksSniffStartToAnalyze, :);

                % Skip this window if some amplitudes not valid
                if any(isnan(whiskPeaksToAnalyze.amplitude))
                    continue;
                end

                % Increment window count
                iSniffStartWin = iSniffStartWin + 1;

                % Start of sniff start window is the pre-valley time of the 
                %   1st peak to analyze, or if doesn't exist, the sniff start time 
                firstPeakToAnalyze = whiskPeaksToAnalyze(1, :);
                firstPreValleyTime = firstPeakToAnalyze.preValleyTime;
                if ~isnan(firstPreValleyTime)
                    sniffStartWinStartTime = firstPreValleyTime;
                else
                    sniffStartWinStartTime = currentSniffStart;
                end

                % End of sniff start window is the post-valley time of the 
                %   last peak to analyze, or if doesn't exist, the sniff end time
                lastPeakToAnalyze = whiskPeaksToAnalyze(end, :);
                lastPostValleyTime = lastPeakToAnalyze.postValleyTime;
                if ~isnan(lastPostValleyTime)
                    sniffStartWinEndTime = lastPostValleyTime;
                else
                    sniffStartWinEndTime = currentSniffEnd;
                end

                % Find all valleys within the sniff start window
                isValleyInWin = valleyTimes >= sniffStartWinStartTime & valleyTimes <= sniffStartWinEndTime;
                whiskValleysToAnalyze = whiskValleyTable(isValleyInWin, :);

                % Compute the logarithmic decrements of successive whisk peak amplitudes
                peakAmplitudes = whiskPeaksToAnalyze.amplitude;
                if numel(peakAmplitudes) > 1
                    logDecrementsSniff = log(peakAmplitudes(2:end) ./ peakAmplitudes(1:end-1));
                else
                    logDecrementsSniff = nan(nWhisksSniffStartToAnalyze - 1, 1);
                end
                                    
                % Create a one-row table for this window
                newWindow = table(fileNumber, {trialName}, iSwp, iSniffStartWin, sniffStartWinStartTime, sniffStartWinEndTime, ...
                    currentSniffStart, currentSniffEnd, ...
                    {whiskPeaksToAnalyze.peakTime}, {whiskPeaksToAnalyze.peakValue}, {peakAmplitudes}, ...
                    {whiskPeaksToAnalyze.preValleyTime}, {whiskPeaksToAnalyze.postValleyTime}, ...
                    {whiskValleysToAnalyze.valleyTime}, {whiskValleysToAnalyze.valleyValue}, {logDecrementsSniff}, ...
                    'VariableNames', emptySniffTable.Properties.VariableNames);

                % Add to sniff start window rows
                sniffWinRows{end + 1} = newWindow;
            end
        end
    end
end

% Vertically concatenate all found windows into a single table
if iSniffStartWin > 0
    sniffStartWinTable = vertcat(sniffWinRows{:});
else
    sniffStartWinTable = emptySniffTable;
end

% Ending message
post_message('Sniff Start Window Detection complete!');


%% --- BASAL RESPIRATION CYCLE DETECTION ---
% Starting message
post_message('Starting Basal Respiration Cycle Detection ...');

switch basalDefinitionMode
case 'pulseCycle'
    % Define empty table structure for consistent output (Simulation)
    emptyBasalTable = table('Size', [0, 22], 'VariableTypes', ...
        {'double', 'double', 'double', 'double', 'double', 'double', 'double', ...
         'double', 'double', 'double', 'double', 'double', 'double', 'double', ...
         'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'}, ...
        'VariableNames', {'seedNumber', 'pulseNumber', 'basalRespCycleNumber', ...
                          'basalRespCycleStartTime', 'basalRespCycleEndTime', ...
                          'breathOnsetTime', 'isStableBasalWhisks', ...
                          'eventTimeWhiskBefore', 'eventTimeTwoWhisksBefore', ...
                          'preIEIWhiskBefore', 'preIEIWhiskAfter', ...
                          'phaseReset', 'phaseChangeWhisk', 'relativeResetTime', ...
                          'whiskPeakTimes', 'whiskPeakValues', 'whiskPeakAmplitudes', ...
                          'whiskPreValleyTimes', 'whiskPostValleyTimes', ...
                          'whiskValleyTimes', 'whiskValleyValues', 'whiskLogDecrements'});

    % For simulation data: define basal respiration cycles based on pulse cycle durations
    
    % Compute the basal respiration cycle duration threshold in ms
    basalDurationThresholdMs = 1000 / basalFreqThreshold;
    
    % Calculate pulse cycle durations
    cycleDurations = pulseCycleEndTimes - pulseCycleStartTimes;
    
    % Check whether all of the first maxWhisksBasalRespToAnalyze peaks
    %   after each pulse are not NaNs
    isFirstNPeaksValid = cellfun(@(x) is_first_n_peaks_valid(x, maxWhisksBasalRespToAnalyze), peakAmplitudesByPulse);
    
    % Filter pulse cycles for basal respiration analysis
    %   A basal respiration cycle must have:
    %   1) Duration greater than basalDurationThresholdMs
    %   2) Contains at least minWhisksBasalRespToAnalyze whisk peaks to be analyzed
    %   3) All whisk peaks to be analyzed have valid amplitudes
    isBasalCycle = cycleDurations >= basalDurationThresholdMs & ...
                   nPeaksByPulse >= minWhisksBasalRespToAnalyze & ...
                   isFirstNPeaksValid;
    
    pulseNumbersBasal = pulseNumbers(isBasalCycle);
    peakTablesForBasal = peakTablesByPulse(isBasalCycle);
    cycleStartsForBasal = pulseCycleStartTimes(isBasalCycle);
    cycleEndsForBasal = pulseCycleEndTimes(isBasalCycle);
    nBasals = sum(isBasalCycle);

    % Decide on the whisk feature for phase definition
    switch whiskDirForPhase
        case 'protraction'
            whiskEventTimes = allValleyTimes;
        case 'retraction'
            whiskEventTimes = allPeakTimes;
    end
    % Get the breath onset times (ms)
    breathOnsetTimes = cycleStartsForBasal + breathOnsetLatencyMs;

    % Get stability for each basal cycle
    if nBasals > 0
        isStableBasalWhisks = [false; diff(pulseNumbersBasal) == 1];
    else
        isStableBasalWhisks = false(size(breathOnsetTimes));
    end

    % Process the identified basal respiration cycles
    basalCycleRows = cell(nBasals, 1);
    for iBasal = 1:nBasals
        % Get peaks for this basal respiration cycle
        peakTableToConsider = peakTablesForBasal{iBasal};
        
        % Truncate to maxWhisksBasalRespToAnalyze
        if height(peakTableToConsider) > maxWhisksBasalRespToAnalyze
            peakTableThis = peakTableToConsider(1:maxWhisksBasalRespToAnalyze, :);
        else
            peakTableThis = peakTableToConsider;
        end
        
        % Define basal respiration cycle start as the preceding valley
        startTime = peakTableThis.preValleyTime(1);
        if isnan(startTime), startTime = cycleStartsForBasal(iBasal); end
        
        % Define basal respiration cycle end as the succeeding valley
        endTime = peakTableThis.postValleyTime(end);
        if isnan(endTime), endTime = cycleEndsForBasal(iBasal); end
        
        % Find all valleys within this basal cycle
        isValleyInCycle = allValleyTimes >= startTime & allValleyTimes <= endTime;
        whiskValleyTimes = allValleyTimes(isValleyInCycle);
        whiskValleyValues = allValleyValues(isValleyInCycle);

        % Extract other whisk data for the table
        whiskPeakTimes = peakTableThis.peakTime;
        whiskPeakValues = peakTableThis.peakValue;
        whiskPeakAmplitudes = peakTableThis.amplitude;
        whiskPreValleyTimes = peakTableThis.preValleyTime;
        whiskPostValleyTimes = peakTableThis.postValleyTime;

        % Compute the logarithmic decrements
        if numel(whiskPeakAmplitudes) > 1
            whiskLogDecrements = log(whiskPeakAmplitudes(2:end) ./ whiskPeakAmplitudes(1:end-1));
        else
            whiskLogDecrements = nan(height(peakTableThis) - 1, 1);
        end

        % Compute phase response for this cycle
        [phaseReset, phaseChangeWhisk, relativeResetTime, preIEIWhiskBefore, ...
         preIEIWhiskAfter, eventTimeWhiskBefore, eventTimeTwoWhisksBefore] = ...
            compute_phase_response(whiskEventTimes, breathOnsetTimes(iBasal), ...
                                   isStableBasalWhisks(iBasal));

        % Package the data for this cycle into a one-row table
        newCycle = table(seedNumber, ...                % seedNumber
                    pulseNumbersBasal(iBasal), ...      % pulseNumber
                    iBasal, ...                         % basalRespCycleNumber
                    startTime, ...                      % basalRespCycleStartTime
                    endTime, ...                        % basalRespCycleEndTime
                    breathOnsetTimes(iBasal), ...       % breathOnsetTime
                    isStableBasalWhisks(iBasal), ...    % isStableBasalWhisks
                    eventTimeWhiskBefore, ...           % eventTimeWhiskBefore
                    eventTimeTwoWhisksBefore, ...       % eventTimeTwoWhisksBefore
                    preIEIWhiskBefore, ...              % preIEIWhiskBefore
                    preIEIWhiskAfter, ...               % preIEIWhiskAfter
                    phaseReset, ...                     % phaseReset
                    phaseChangeWhisk, ...               % phaseChangeWhisk
                    relativeResetTime, ...              % relativeResetTime
                    {whiskPeakTimes}, ...               % whiskPeakTimes
                    {whiskPeakValues}, ...              % whiskPeakValues
                    {whiskPeakAmplitudes}, ...          % whiskPeakAmplitudes
                    {whiskPreValleyTimes}, ...          % whiskPreValleyTimes
                    {whiskPostValleyTimes}, ...         % whiskPostValleyTimes
                    {whiskValleyTimes}, ...             % whiskValleyTimes
                    {whiskValleyValues}, ...            % whiskValleyValues
                    {whiskLogDecrements}, ...           % whiskLogDecrements
                    'VariableNames', emptyBasalTable.Properties.VariableNames);

        basalCycleRows{iBasal} = newCycle;
    end

case 'respPeak'
    % Define empty table structure for consistent output (Experiment)
    emptyBasalTable = table('Size', [0, 33], 'VariableTypes', ...
        {'double', 'cellstr', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', ...
         'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', ...
         'double', 'double', 'double', 'double', 'double', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'}, ...
        'VariableNames', {'fileNumber', 'trialName', 'sweepNumber', 'basalRespCycleNumber', 'basalRespPeakNumber', ...
                          'basalRespCycleStartTime', 'basalRespCycleEndTime', ...
                          'basalRespPeakTime', 'basalRespPeakValue', 'basalRespPreValleyTime', ...
                          'basalRespPostValleyTime', 'basalRespSucceedingIPI', ...
                          'isWithinBasal', 'inspWhiskPeakNum', 'inspWhiskPeakTime', 'inspWhiskPreValleyTime', ...
                          'breathOnsetTime', 'isStableBasalWhisks', 'eventTimeWhiskBefore', 'eventTimeTwoWhisksBefore', 'preIEIWhiskBefore', ...
                          'preIEIWhiskAfter', 'phaseReset', 'phaseChangeWhisk', 'relativeResetTime', ...
                          'whiskPeakTimes', 'whiskPeakValues', 'whiskPeakAmplitudes', ...
                          'whiskPreValleyTimes', 'whiskPostValleyTimes', ...
                          'whiskValleyTimes', 'whiskValleyValues', 'whiskLogDecrements'});

    % For experimental data: use pre-computed basal respiration peak table

    % Convert to seconds
    breathOnsetLatency = breathOnsetLatencyMs / 1000;

    basalCycleRows = {};
    for iSwp = 1:nSweeps
        % Filter the file-level basal peak table for the current sweep
        if ~isempty(basalRespPeakTable)
            basalRespPeaks = basalRespPeakTable(basalRespPeakTable.sweepNumber == iSwp, :);
        else
            basalRespPeaks = table();
        end
        
        % Extract whisk data for this sweep
        whiskPeaks = whiskPeakTables{iSwp};
        whiskValleys = whiskValleyTables{iSwp};
        whiskPeakTimesThisSwp = whiskPeaks.peakTime;
        whiskPreValleyTimesThisSwp = whiskPeaks.preValleyTime;
        whiskValleyTimesThisSwp = whiskValleys.valleyTime;

        % Skip this sweep if no whisks or basal respiration peaks
        if isempty(basalRespPeaks) || isempty(whiskPeaks)
            continue;
        end

        % Count the number of basal respiration peaks for this sweep
        nBasalRespThisSwp = height(basalRespPeaks);

        % 1. Find inspiratory whisks within each basal respiration peak
        whiskPeakNumbers = nan(nBasalRespThisSwp, 1);
        for iResp = 1:nBasalRespThisSwp
            respPeakTime = basalRespPeaks.basalRespPeakTimes(iResp);
            respPreValley = basalRespPeaks.basalRespPreValleyTimes(iResp);
            respPostValley = basalRespPeaks.basalRespPostValleyTimes(iResp);

            % Handle NaN in postValleyTime
            if isnan(respPostValley)
                respPostValley = tVecs(end, iSwp);
            end

            % Moore et al 2013 definition: Find the first whisk peak with pre-valley 
            %   after breathOnsetLatency before respiration pre-valley and
            %   before respiration peak
            firstPreValleyNumber = ...
                find(whiskPreValleyTimesThisSwp >= respPreValley - breathOnsetLatency & ...
                        whiskPreValleyTimesThisSwp < respPeakTime, 1 , 'first');
            
            % Find whisk peaks within the respiration valley-to-valley window
            possibleWhiskPeakNumbers = ...
                find(whiskPeakTimesThisSwp >= respPreValley & ...
                        whiskPeakTimesThisSwp < respPostValley);
        
            % Find the whisk peak closest to the basal respiration peak
            if ~isempty(possibleWhiskPeakNumbers)
                [~, closestIdxInPossible] = min(abs(whiskPeakTimesThisSwp(possibleWhiskPeakNumbers) - respPeakTime));
                closestWhiskPeakNumber = possibleWhiskPeakNumbers(closestIdxInPossible);
            else
                closestWhiskPeakNumber = [];
            end
            

            % If the closest whisk peak is different from the one with 
            %   the first prevalley, print message
            if ~isempty(firstPreValleyNumber) && ~isempty(closestWhiskPeakNumber) && firstPreValleyNumber ~= closestWhiskPeakNumber
                fprintf(fileIDDiffs, ['First whisk protraction (peak number %d) within basal respiration is different ', ...
                         'from closest peak (peak number %d) for basal respiration cycle #%d in sweep %d of file %d!!\n\n'], ...
                         firstPreValleyNumber, closestWhiskPeakNumber, iResp, iSwp, fileNumber);
            end

            % Skip this respiration if no whisks within the basal respiration peak
            if isempty(closestWhiskPeakNumber)
                continue;
            end

            % Set the closest whisk peak as the 'inspiratory whisk'
            whiskPeakNumbers(iResp) = closestWhiskPeakNumber;
        end

        % Add inspiratory whisk peak number to basal respiratory peak table for this sweep
        basalRespPeaks.whiskPeakNumber = whiskPeakNumbers;

        % Check if inspiratory whisk found for each row
        hasInspiratoryWhisk = ~isnan(whiskPeakNumbers);

        % Skip this sweep if no inspiratory whisk found
        if ~any(hasInspiratoryWhisk)
            continue;
        end

        % Get the peak numbers for inspiratory whisks found
        peakNumbersInspWhisk = whiskPeakNumbers(hasInspiratoryWhisk);

        % Remove rows with no inspiratory whisk found)
        inspWhiskTableTemp = basalRespPeaks(hasInspiratoryWhisk, :);

        % Add associated whisk peak information
        inspWhiskTable = horzcat(inspWhiskTableTemp, whiskPeaks(peakNumbersInspWhisk, :));

        % Count the number of inspiratory whisks
        nInspWhisk = height(inspWhiskTable);

        % Get all inspiratory whisk times this sweep
        inspWhiskTimesThisSwp = inspWhiskTable.peakTime;

        % 2. Define basal respiration cycles
        iCycle = 0;             % Counter for basal respiration cycles for this sweep
        lastRespPeakNum = NaN;  % Last respiration peak number for this sweep
        for iInsp = 1:nInspWhisk       
            % Get the current inspiratory whisk info
            currentInspWhisk = inspWhiskTable(iInsp, :);
            inspWhiskPeakNum = currentInspWhisk.whiskPeakNumber;
            inspWhiskPeakTime = currentInspWhisk.peakTime;
            inspWhiskPreValleyTime = currentInspWhisk.preValleyTime;
            respPeakNumber = currentInspWhisk.respPeakNumber;
            basalRespPeakNumber = currentInspWhisk.basalRespPeakNumber;
            basalRespPeakTime = currentInspWhisk.basalRespPeakTimes;
            basalRespPeakValue = currentInspWhisk.basalRespPeakValues;
            basalRespPreValleyTime = currentInspWhisk.basalRespPreValleyTimes;
            basalRespPostValleyTime = currentInspWhisk.basalRespPostValleyTimes;
            basalRespSucceedingIPI = currentInspWhisk.basalRespSucceedingIPIs;
            isWithinBasal = currentInspWhisk.isWithinBasal;

            % Find the next inspiratory whisk time
            if iInsp < nInspWhisk
                nextInspWhiskTime = inspWhiskTimesThisSwp(iInsp + 1);
            else
                % No next inspiratory whisk, set at end of sweep
                nextInspWhiskTime = tVecs(end, iSwp);
            end
                  
            % The search end time is the earliest of the next insp whisk or 
            %   the end of the current respiration (the succeeding valley)
            searchEndTime = min(nextInspWhiskTime, basalRespPostValleyTime);
            
            % Find intervening whisks after the current inspiratory whisk
            %   and before the search end time
            interWhiskPeakNumbers = ...
                find(whiskPeakTimesThisSwp > inspWhiskPeakTime & ...
                    whiskPeakTimesThisSwp < searchEndTime);

            % Define a basal respiration cycle to be analyzed as 
            %   one in which there are enough intervening whisk for a cycle
            %   and that all whisk peaks have valid amplitudes
            if numel(interWhiskPeakNumbers) >= minWhisksBasalRespToAnalyze - 1
                % Collect all whisks for this cycle (inspiratory + all intervening)
                whiskIndicesForCycle = [inspWhiskPeakNum; interWhiskPeakNumbers];
                whiskPeaksToAnalyze = whiskPeaks(whiskIndicesForCycle, :);

                % Skip this inspiratory whisk if some whisk peak
                %   does not have a valid amplitude
                if any(isnan(whiskPeaksToAnalyze.amplitude))
                    continue;
                end

                % Increment cycle number and save 
                iCycle = iCycle + 1;

                % Cycle start time is the preceding valley of the inspiratory whisk
                %   or if not present, the preceding valley of the basal respiration
                if ~isnan(inspWhiskPreValleyTime)
                    basalRespCycleStartTime = inspWhiskPreValleyTime;
                else
                    basalRespCycleStartTime = basalRespPreValleyTime;
                end

                % Cycle end time is the succeeding valley of the last intervening whisk
                %   or if not present, the search end time
                lastWhiskInCycle = whiskPeaksToAnalyze(end, :);
                lastWhiskPostValleyTime = lastWhiskInCycle.postValleyTime;
                if ~isnan(lastWhiskPostValleyTime)
                    basalRespCycleEndTime = lastWhiskPostValleyTime;
                else
                    basalRespCycleEndTime = searchEndTime;
                end
                
                % Find valleys within the cycle
                isValleyInCycle = whiskValleyTimesThisSwp >= basalRespCycleStartTime & ...
                                  whiskValleyTimesThisSwp <= basalRespCycleEndTime;
                whiskValleysToAnalyze = whiskValleys(isValleyInCycle, :);

                % Compute log decrements
                %   diff(log([A1, A2])) == log(A2) - log(A1) == log(A2/A1)
                logDecrements = diff(log(whiskPeaksToAnalyze.amplitude));
                
                % Compute the breath onset time
                breathOnsetTime = basalRespPreValleyTime - breathOnsetLatency;

                % Decide on whisk direction used for phase response
                switch whiskDirForPhase
                case 'protraction'
                    whiskEventTimes = whiskValleyTimesThisSwp;
                case 'retraction'
                    whiskEventTimes = whiskPeakTimesThisSwp;
                end

                % Compute phase response if the current basal respiration cycle
                %       is preceded by a basal respiration cycle (so that both resp cycles and whisk cycles are stable)
                %   phase reset: the phase of breath onset within the whisking cycle
                %   delta phase whisk: the change in whisk phase associated with the breath
                isStableBasalWhisks = respPeakNumber == lastRespPeakNum + 1;
                [phaseReset, phaseChangeWhisk, relativeResetTime, preIEIWhiskBefore, ...
                    preIEIWhiskAfter, eventTimeWhiskBefore, eventTimeTwoWhisksBefore] = ...
                    compute_phase_response(whiskEventTimes, breathOnsetTime, isStableBasalWhisks);

                % Create a one-row table for this cycle
                newCycle = table(fileNumber, {trialName}, iSwp, iCycle, basalRespPeakNumber, ...
                    basalRespCycleStartTime, basalRespCycleEndTime, ...
                    basalRespPeakTime, basalRespPeakValue, basalRespPreValleyTime, basalRespPostValleyTime, basalRespSucceedingIPI, ...
                    isWithinBasal, inspWhiskPeakNum, inspWhiskPeakTime, inspWhiskPreValleyTime, ...
                    breathOnsetTime, isStableBasalWhisks, eventTimeWhiskBefore, eventTimeTwoWhisksBefore, preIEIWhiskBefore, ...
                    preIEIWhiskAfter, phaseReset, phaseChangeWhisk, relativeResetTime, ...
                    {whiskPeaksToAnalyze.peakTime}, {whiskPeaksToAnalyze.peakValue}, ...
                    {whiskPeaksToAnalyze.amplitude}, {whiskPeaksToAnalyze.preValleyTime}, {whiskPeaksToAnalyze.postValleyTime}, ...
                    {whiskValleysToAnalyze.valleyTime}, {whiskValleysToAnalyze.valleyValue}, {logDecrements}, ...
                    'VariableNames', emptyBasalTable.Properties.VariableNames);

                basalCycleRows{end + 1} = newCycle;

                % Save current resp peak number as last basal cycle respiratory peak number
                lastRespPeakNum = respPeakNumber;
            end
        end
    end
end

% Vertically concatenate all found cycles into a single table
if ~isempty(basalCycleRows)
    basalRespCycleTable = vertcat(basalCycleRows{:});
else
    basalRespCycleTable = emptyBasalTable;
end

% Ending message
post_message('Basal Respiration Cycle Detection complete!');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

