function [averagePscTrace, allPscTraces] = ...
                compute_average_PSC_trace (pscInfo, data, varargin)
%% Compute the "average PSC trace"
% Usage: [averagePscTrace, allPscTraces] = ...
%               compute_average_PSC_trace (pscInfo, data, varargin)
% Explanation: 
%       TODO: takes in the sorted epscs array from freEpsc analysis and then returns
%       a data array containing the values for the raw trace.  traces are
%       corrected to baseline based on the first ms or so of data.
%
% Outputs:
%       averagePscTrace     - "averaged PSC trace"
%                           specified as a column vector with 
%                                   a length of TraceLengthSamples
%       allPscTraces        - all PSC traces that were averaged
%                           specified as a matrix with each column a PSC trace
%
% Arguments:    
%       pscInfo     - information for found post-synaptic currents
%                       Each row corresponds to a PSC
%                       Column assignments for pscInfo:
%                           1 = index at PSC breakpoint
%                           2 = index at PSC peak
%                           3 = smoothed data value at PSC breakpoint
%                           4 = smoothed data value at PSC peak
%                           5 = smoothed amplitude of the PSC
%                           6 = 0-100% rise time (samples)
%                           7 = 10-90% rise time (samples)
%                           8 = inter-event interval from this PSC peak 
%                               to next PSC peak (samples)
%                           9 = 50% decay time (samples)
%                   must be a numeric array with at least 8 columns
%       data        - vector of current trace data
%                   must be a numeric vector
%       varargin    - 'TraceLengthSamples': total PSC trace length (samples)
%                   must be a positive integer scalar
%                   default == 500
%                   - 'BeforePeakSamples': trace length before breakpoint 
%                                           for each PSC (samples)
%                   must be a positive integer scalar
%                   default == 30
%                   - 'SiMs': sampling interval in ms/sample
%                   must be a numeric positive scalar
%                   default == []
%                   - 'TraceLengthMs': total trace length for each PSC (ms)
%                   must be a numeric positive scalar
%                   default == []
%                   - 'BeforePeakMs': trace length before breakpoint 
%                                       for each PSC (ms)
%                   must be a numeric positive scalar
%                   default == []
%                   - 'DealWithTooShort': what to do for PSCs that are too short
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'none'      - include original trace and do nothing,
%                                           or omit if at boundaries
%                       'padboth'   - pad short PSCs with NaNs on both sides
%                       'padright'  - pad short PSCs with NaNs on right side only
%                       'omit'      - omit those short PSCs
%                   default == 'none'
%                   - 'MessageMode' - how message boxes are shown
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'wait'  - stops program and waits for the user
%                                   to close the message box
%                       'show'  - does not stop program but still show the
%                                   message box
%                       'none'  - neither stop program nor show a message box
%                   default == 'wait'
%                   - 'Verbose' - whether to print to standard output
%                                   regardless of message mode
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   
%
% Requires:
%       /home/Matlab/Adams_Functions/print_or_show_message.m
%
% Used by:
%       cd/minEASE_detect_gapfree_events.m

% File History:
% ---------- Created by Koji Takahashi & Mark P Beenhakker
% 2017-05-25 AL - Renamed getAverageTrace.m -> compute_average_PSC_trace.m
% 2017-05-25 AL - Organized, simplified and annotated code
% 2017-05-25 AL - Now returns empty matrices if no PSC satisfies criteria
% 2017-05-25 AL - Added input parser scheme
% 2017-05-25 AL - Commented out subtraction of "baseline trace"
% 2017-07-24 AL - Align at peaks instead
% 2017-07-24 AL - Changed criteria for trimming
% 2017-07-24 AL - Add option of padding by NaNs
% 2017-07-24 AL - Added 'DealWithTooShort' as an optional argument
% 2017-07-24 AL - The parameters 'BeforeBreakMs' and 'BeforeBreakSamples' 
%                   are changed to 'BeforePeakMs' and 'BeforePeakSamples'
% 2017-10-23 AL - Fixed traceBeginIdx and traceEndIdx to account for data ends
% 2018-02-02 AL - Added possibleDealWithTooShort
% 2018-02-02 AL - Added showMessage as an optional parameter-value pair argument
% 2018-02-02 AL - Now uses print_or_show_message.m for output
% 2018-02-07 MD - Changed usage of print_or_show_message()
% 2018-02-27 AL - Changed showMessages to messageMode with possible values:
%                   'wait', 'show', 'none'
% 2018-03-02 MD - Defined verbose parameter for print_or_show_message
% TODO for undergrad: make this a more general function 
%                     compute_average_event_trace.m that passes 
%                     idxEventStarts, interEventIntervals, data as required args
%                     and siMs, traceLengthMs, beforeStartMs as optional args
%                     and returns averageEventTrace, allEventTraces

%% Hard-coded parameters
PSCINFO_NCOLS = 11;                 % pscInfo should have this many columns
possibleDealWithTooShort = {'none', 'padboth', 'padright', 'omit'};
validMessageModes = {'wait', 'show', 'none'};

%% Default values for optional arguments
traceLengthSamplesDefault = 500;    % default total PSC trace length (samples)
                                    % TODO: why 500?
beforePeakSamplesDefault = 30;      % default trace length before breakpoint 
                                    %       for each PSC (samples)
                                    % TODO: why 30?
siMsDefault = [];
traceLengthMsDefault = [];
beforePeakMsDefault = [];
dealWithTooShortDefault = 'none';   % do not omit or pad traces by default
messageModeDefault = 'none';        % print to standard output by default

verboseDefault = false;             % default: Program does not print message
                                    %   even if message box is shown

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help compute_average_PSC_trace'' for usage']);
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'compute_average_PSC_trace';

% Add required inputs to an Input Parser
addRequired(iP, 'pscInfo', ...                  % information for found PSCs
    @(x) assert(isnumeric(x) && size(x, 2) >= PSCINFO_NCOLS, ...
            sprintf(['pscInfo should be a numeric array', ...
                    ' with at least %d columns'], PSCINFO_NCOLS)));
addRequired(iP, 'data', ...                     % vector of current data
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TraceLengthSamples', traceLengthSamplesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'BeforePeakSamples', beforePeakSamplesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'SiMs', siMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'TraceLengthMs', traceLengthMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'BeforePeakMs', beforePeakMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'DealWithTooShort', dealWithTooShortDefault, ...
    @(x) any(validatestring(x, possibleDealWithTooShort)));
addParameter(iP, 'MessageMode', messageModeDefault, ...
    @(x) any(validatestring(x, validMessageModes)));
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, pscInfo, data, varargin{:});
traceLengthSamples = iP.Results.TraceLengthSamples;
beforePeakSamples = iP.Results.BeforePeakSamples;
siMs = iP.Results.SiMs;
traceLengthMs = iP.Results.TraceLengthMs;
beforePeakMs = iP.Results.BeforePeakMs;
dealWithTooShort = validatestring(iP.Results.DealWithTooShort, ...
                                    possibleDealWithTooShort);
messageMode = validatestring(iP.Results.MessageMode, validMessageModes);
verbose = iP.Results.Verbose;

%% Convert trace lengths to samples if given in ms
%   Note: sampling interval in ms (siMs) must be provided
if ~isempty(traceLengthMs) && ~isempty(siMs)
    traceLengthSamples = round(traceLengthMs / siMs);
end
if ~isempty(beforePeakMs) && ~isempty(siMs) 
    beforePeakSamples = round(beforePeakMs / siMs);
end

% Compute number of data points
nPoints = length(data);

%% Modify pscInfo according to dealWithTooShort
switch dealWithTooShort
case 'none'
    % Only PSCs that satisfies the following criteria are averaged:
    %   (1) peak index - 1 no less than beforePeakSamples
    %   (2) total indices - peak index + 1 no shorter than 
    %               traceLengthSamples - beforePeakSamples
    pscInfoMod = pscInfo(pscInfo(:, 2) - 1 >= beforePeakSamples & ...
                    nPoints - pscInfo(:, 2) + 1 >= ...
                    traceLengthSamples - beforePeakSamples, :);
case {'padboth', 'padright'}
    % Don't omit any PSC
    pscInfoMod = pscInfo;
case 'omit'
    % Only PSCs that satisfies the following criteria are averaged:
    %   (1) breakpoint to peak not longer than beforePeakSamples
    %   (2) peak to next breakpoint at least as long as 
    %               traceLengthSamples - beforePeakSamples
    pscInfoMod = pscInfo(pscInfo(:, 6) <= beforePeakSamples & ...
                    pscInfo(:, 9) >= traceLengthSamples - beforePeakSamples, :);
otherwise
end

% Extract from pscInfoMod
idxPscBreaks = pscInfoMod(:, 1);            % indices of PSC breakpoints
idxPscPeaks = pscInfoMod(:, 2);             % indices of PSC peaks
interStimulusIntervals = pscInfoMod(:, 9);  % peaks to next breakpoints
nPscs = length(idxPscBreaks);               % number of PSCs left

% If no PSCs left, return empty matrices and exit function
if nPscs == 0
    allPscTraces = [];
    averagePscTrace = [];
    message = sprintf(['No PSCs satisfying criteria ', ...
                        'in the %s mode are found!'], ...
                        dealWithTooShort);
    mTitle = 'Average PSC Trace Warning';
    icon = 'warn';
    print_or_show_message(message, 'MessageMode', messageMode, ...
                            'MTitle', mTitle, 'Icon', icon, 'Verbose', verbose);
    return;
end

%% Compute the "average PSC trace"
% Place all PSC traces together in a matrix (each column is a PSC trace)
allPscTraces = zeros(traceLengthSamples, nPscs);
for i = 1:nPscs
    % Find index for beginning of trace
    traceBeginIdxIdeal = idxPscPeaks(i) - beforePeakSamples;
    traceBeginIdx = max(1, traceBeginIdxIdeal);

    % Find index for end of trace
    traceEndIdxIdeal = traceBeginIdxIdeal + (traceLengthSamples - 1);
    traceEndIdx = min(traceEndIdxIdeal, nPoints);
    
    % Place PSC trace in the ith column of allPscTraces
    switch dealWithTooShort
    case {'none', 'omit'}
        % Take the traces directly from data in these cases
        %   Omit traces that are too short
        if traceEndIdx - traceBeginIdx + 1 == traceLengthSamples
            allPscTraces(:, i) = data(traceBeginIdx:traceEndIdx);
        else
            allPscTraces(:, i) = NaN * ones(traceLengthSamples, 1);                    
            message = sprintf(['The %dth trace was omitted because ', ...
                            'it did not have the correct length!\n'], i);
            mTitle = 'Average PSC Trace Warning';
            icon = 'warn';
            print_or_show_message(message, 'MessageMode', messageMode, ...
                                    'MTitle', mTitle, 'Icon', icon, ...
                                    'Verbose', verbose);
        end
    case {'padboth', 'padright'}
        % If breakpoint to peak shorter than beforePeakSamples, 
        %   pad with NaNs
        pscBeginIdx = max(idxPscBreaks(i), traceBeginIdx);  % PSC start index
        nSamplesToPadLeft = traceLengthSamples - ...
                                (traceEndIdxIdeal - pscBeginIdx + 1);
                                            % # of samples to pad on the left
        if nSamplesToPadLeft > 0
            padLeft = NaN * ones(nSamplesToPadLeft, 1);
        else
            padLeft = [];
        end

        % If peak to next breakpoint shorter than 
        %               traceLengthSamples - beforePeakSamples, pad with NaNs
        pscEndIdx = min(idxPscPeaks(i) + interStimulusIntervals(i) - 1, ...
                            traceEndIdx);               % PSC end index
        nSamplesToPadRight = traceLengthSamples - ...
                                (pscEndIdx - traceBeginIdxIdeal + 1);
                                            % # of samples to pad on the right
        if nSamplesToPadRight > 0
            padRight = NaN * ones(nSamplesToPadRight, 1);
        else
            padRight = [];
        end

        % Combine the padding and the PSC trace and place in the matrix
        if strcmpi(dealWithTooShort, 'padboth') || traceBeginIdx < 1
            allPscTraces(:, i) = [padLeft; ...
                                  data(pscBeginIdx:pscEndIdx); ...
                                  padRight];
        elseif strcmpi(dealWithTooShort, 'padright')
            if traceBeginIdxIdeal < traceBeginIdx
                padLeft = NaN * ones(traceBeginIdx - traceBeginIdxIdeal, 1);
                allPscTraces(:, i) = [padLeft; ...
                                      data(traceBeginIdx:pscEndIdx); ...
                                      padRight];
            else
                allPscTraces(:, i) = [data(traceBeginIdx:pscEndIdx); ...
                                      padRight];
            end
        end
    otherwise
    end
end

%{ 
    % TODO: Why is this necessary? Why 10 to 20?
    baselineTrace = mean(allPscTraces(10:20, i));
    allPscTraces(:, i) = allPscTraces(:, i) - baselineTrace;
%}

% Compute the average trace over all PSC traces
switch dealWithTooShort
case {'none', 'omit'}
    % The averaged trace is a column vector
    averagePscTrace = mean(allPscTraces, 2);
case {'padboth', 'padright'}
    % Ignore NaNs when averaging
    averagePscTrace = nanmean(allPscTraces, 2);
otherwise
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%% Trim pscInfo: Only PSCs that satisfies the following criteria are averaged:
%   (1) peak to next peak at least as long as traceLengthSamples and 
%   (2) breakpoint index lies after beforeBreakSamples
%   TODO: Shouldn't we use the breakpoint to next breakpoint instead? Or line up at peaks instead?
pscInfoTrimmed = pscInfo(pscInfo(:, 8) >= traceLengthSamples & ...
                            pscInfo(:, 1) > beforeBreakSamples, :);

%                   - 'BeforeBreakSamples': trace length before breakpoint 
%                                           for each PSC (samples)
%                   must be a positive integer scalar
%                   default == 30
%                   - 'BeforePeakMs': trace length before breakpoint 
%                                       for each PSC (ms)
%                   must be a numeric positive scalar
%                   default == []
beforeBreakSamplesDefault = 30;     % default trace length before breakpoint 
                                    %        for each PSC (samples)
                                    % TODO: why 30?
beforePeakMsDefault = [];
addParameter(iP, 'BeforeBreakSamples', beforeBreakSamplesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'BeforePeakMs', beforePeakMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
beforeBreakSamples = iP.Results.BeforeBreakSamples;
beforePeakMs = iP.Results.BeforePeakMs;
if ~isempty(beforePeakMs) && ~isempty(siMs) 
    beforeBreakSamples = round(beforePeakMs / siMs);
end

idxPscBreaks = pscInfoTrimmed(:, 1);            % indices of PSC breakpoints

    % Find index for beginning of trace
    %   Note: By criterion (2) this will always be >= 1
    beginIdx = idxPscBreaks(i) - beforeBreakSamples;

    % Find index for end of trace
    %   Note: By criterion (1) this will always be <= length(data)
    endIdx = beginIdx + (traceLengthSamples - 1);

    fprintf(['No PSCs with peak to next peak greater than', ...
                ' or equal to the trace length are found!\n\n']);

        nSamplesToPadLeft = max(pscBeginIdx - traceBeginIdx, ...
                                traceLengthSamples - ...
                                (traceEndIdx - pscBeginIdx + 1));
        nSamplesToPadRight = max(traceEndIdx - pscEndIdx, ...
                                traceLengthSamples - ...
                                (pscEndIdx - traceBeginIdx + 1));

fprintf('No PSCs satisfying criteria in the %s mode are found!\n\n', ...
            dealWithTooShort);
fprintf(['The %dth trace is omitted because ', ...
    'it does not have the correct length!\n'], i);

            print_or_show_message(showMessage, message, ...
                            'MTitle', mTitle, 'Icon', icon);
%                   - 'ShowMessage': whether to show messages in messages boxes
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
showMessageDefault  = false;        % print to standard output by default
addParameter(iP, 'ShowMessage', showMessageDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
showMessage = iP.Results.ShowMessage;
            if%                   - 'ShowMessage': whether to show messages in messages boxes
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

%}
