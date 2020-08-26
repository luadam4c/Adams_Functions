function [eventInfoExpanded, varargout] = ...
                identify_bursts (eventInfo, varargin)
%% Identify bursts from a list of directional events
% Usage: [eventInfoExpanded, varargout] = ...
%               identify_bursts (eventInfo, varargin)
% Explanation: 
%       Use baseline differences between adjacent events to detect regions
%       that includes a burst. 
%       First, use a baseline difference criterion to identify
%       a large "crude burst region" of a given minimum region size
%       Then use the maximum interstimulus interval and minimum
%       number of spikes per burst to identify "bursts"
%       TODO: Problem 1: If two burst windows overlap, only one burst will
%                           be identifies instead of 2
% Outputs:
%       eventInfoExpanded       - expanded information for directional events
%                               Each row corresponds to a directional event
%                               Column assignments for eventInfoExpanded:
%                                   1 = index at event breakpoint
%                                   2 = index at event peak
%                                   3 = smoothed data value at event breakpoint
%                                   4 = smoothed data value at event peak
%                                   5 = smoothed amplitude of the event
%                                   6 = 0-100% rise time (samples)
%                                   7 = 10-90% rise time (samples)
%                                   8 = inter-event interval from this event 
%                                       peak to next event peak (samples)
%                                   9 = interstimulus interval from this event 
%                                       peak to next event breakpoint (samples)
%                                   10 = 50% decay time (samples)
%                                   11 = "full decay" time (samples):
%                                       time to return within noiseLevel 
%                                       of breakpoint value
%                                   12 = whether in a burst
%                                   13 = whether in a "crude burst region"
%                                   14 = burst region number
%                               specified as a numeric array with 13 columns
%       varargout               {1} = absDiffBreaks
%                               {2} = willJump
%                               {3} = inCrude
%                               {4} = RegionNo
%                               {5} = inBurst
% Arguments:    
%       eventInfo       - information for found directional events
%                       Each row corresponds to a directional event
%                       Column assignments for eventInfo:
%                           1 = index at event breakpoint
%                           2 = index at event peak
%                           3 = smoothed data value at event breakpoint
%                           4 = smoothed data value at event peak
%                           5 = smoothed amplitude of the event
%                           6 = 0-100% rise time (samples)
%                           7 = 10-90% rise time (samples)
%                           8 = inter-event interval from this event 
%                               peak to next event peak (samples)
%                           9 = interstimulus interval from this event 
%                               peak to next event breakpoint (samples)
%                           10 = 50% decay time (samples)
%                           11 = "full decay" time (samples):
%                                   time to return within noiseLevel 
%                                   of breakpoint value
%                       must be a numeric array with at least 9 columns
%       varargin    - 'MinBaselineDiff': minimum baseline difference
%                           (data units) for identifying a "crude burst region"
%                   must be a positive scalar
%                   default == 40
%                   - 'CrudeRegionSize': crude burst region size (events)
%                   must be a positive integer scalar
%                   default == 50
%                   - 'MinSpikesPerBurst': minimum number of spikes in a burst
%                   must be a positive integer scalar
%                   default == 8
%                   - 'MaxIsi': maximum interstimulus interval (samples) 
%                   must be a positive integer scalar
%                   default == 500
%
% Used by:
%       cd/minEASE_detect_gapfree_events.m

% File History:
% ---------- Created by Koji Takahashi & Mark P Beenhakker
% 2017-05-24 AL - Renamed burstFinder.m -> identify_bursts.m
% 2017-05-24 AL - Added input parser scheme
% 2017-05-26 AL - Simplified code (by a lot)
% 2017-05-26 AL - Now includes the last spike in the burst
% 2017-06-02 AL - Removed refined burst region
% 2017-06-02 AL - Identify all bursts that satisfy the two criteria within 
%                   a burst region
% 2017-06-05 AL - Changed inter-event intervals to interstimulus intervals
% TODO: Plot the burst analysis here
% TODO: Implement plotflag

%% Constants
EVENTINFO_NCOLS = 11;        % eventInfo should have this many columns

%% Default parameters
minBaselineDiffDefault = 40;    % Default minimum baseline difference (d. u.)
                                %   for identifying a crude burst region
                                % TODO: Why 40 pA?
crudeRegionSizeDefault = 50;    % Default burst region size (events)
                                % TODO: Why 50 spikes?
minSpikesPerBurstDefault = 8;   % Default minimum number of spikes in a burst
                                % TODO: Why 8 spikes?
maxIsiDefault = 500;            % Default maximum interstimulus interval (samples) 
                                %   within a burst
                                % TODO: Why 50 ms?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add required inputs to an Input Parser
addRequired(iP, 'eventInfo', ...                  % information for found PSCs
    @(x) assert(isnumeric(x) && size(x, 2) >= EVENTINFO_NCOLS, ...
            sprintf(['eventInfo should be a numeric array', ...
                    ' with at least %d columns'], EVENTINFO_NCOLS)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MinBaselineDiff', minBaselineDiffDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'CrudeRegionSize', crudeRegionSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'MinSpikesPerBurst', minSpikesPerBurstDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'MaxIsi', maxIsiDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Read from the Input Parser
parse(iP, eventInfo, varargin{:});
minBaselineDiff   = iP.Results.MinBaselineDiff;
crudeRegionSize   = iP.Results.CrudeRegionSize;
minSpikesPerBurst = iP.Results.MinSpikesPerBurst;
maxIsi            = iP.Results.MaxIsi;

%% Extract from eventInfo
valEventBreaks = eventInfo(:, 3);           % values at event breakpoints (pA)
interStimulusIntervals = eventInfo(:, 9);   % interstimulus intervals (samples)
nEvents = length(valEventBreaks);           % number of events

%% Initialize output vectors
willJump = zeros(nEvents, 1);   % whether the event will 
                                %   "jump" to the next by minBaselineDiff
inCrude = zeros(nEvents, 1);    % whether the event is in a crude burst region
RegionNo = zeros(nEvents, 1);   % stores the burst region number for each event;
                                %   0 if not in a crude burst region
inBurst = zeros(nEvents, 1);    % whether the event is in a real burst

%% Identify "crude burst regions"
%   An event is within a "crude burst region" if one of the following is true: 
%       (1) absolute difference between the breakpoint baseline
%           of this event and of the next one is at least minBaselineDiff
% TODO: Why?
%       (2) it is within crudeRegionSize/2 of an event satisfying (1)
% TODO: Why?
%   Note: the last event will never be fulfill (1) but might fulfill 2

% Calculate absolute differences in the event breakpoint baselines 
absDiffBreaks = abs(diff(valEventBreaks));

% In a vector called willJump, 
%   mark events satisfying criterion (1) by 1 and others by 0
willJump(absDiffBreaks >= minBaselineDiff) = 1;

% Apply a moving average filter to willJump to identify all events
%   satisfying criterion (2) as well and mark with 1 in a vector inCrude
willJumpSmoothed = smooth(willJump, crudeRegionSize);
inCrude(willJumpSmoothed > 0) = 1;

% Tally results
rowsInCrude = find(inCrude > 0);     % rows of events in crude burst regions
nEventsInCrude = length(rowsInCrude);% number of events in crude burst regions

% If no bursts are detected, return empty columns and exit function
if isempty(rowsInCrude)             % if at least one burst region is detected
    eventInfoExpanded = zeros(nEvents, 14);
    eventInfoExpanded(:, 1:11) = eventInfo;
    
    varargout{1} = absDiffBreaks;
    varargout{2} = willJump;
    varargout{3} = inCrude;
    varargout{4} = RegionNo;
    varargout{5} = inBurst;

    fprintf('No bursts are detected!\n\n');
    return;
end

% Identify events in a contiguous crude burst region with the same number
RegionNo(rowsInCrude(1)) = 1;           % initialize the first region
for i = 2:nEventsInCrude
    if rowsInCrude(i) - rowsInCrude(i-1) == 1       % still in same region
        RegionNo(rowsInCrude(i)) = RegionNo(rowsInCrude(i-1));
    else                                            % in a different region
        RegionNo(rowsInCrude(i)) = RegionNo(rowsInCrude(i-1)) + 1;
    end
end

% Tally the total number of crude burst regions
nRegions = max(RegionNo);               % number of crude burst regions

%% Identify "real bursts" within the crude burst region
%   A "burst" is the largest group of contiguous events in each refined burst
%   region that satisfies both of:
%       (1) each interstimulus interval to the next event is not more than maxIsi,
%               except for the last event
%       (2) total number of events is at least minSpikesPerBurst
% TODO: Is (2) necessary?
for iRegion = 1:nRegions        % for each burst region
    % Whether the event satisfies criterion (1)
    satisfyIsi = (interStimulusIntervals <= maxIsi & RegionNo == iRegion);

    % Construct a truthStatement to find the first 
    %   and last+1 rows of each run of events satisfying criterion (1)
    truthStatement = diff([0; satisfyIsi; 0]);
    firstRows = find(truthStatement == 1);
    lastPlusOneRows = find(truthStatement == -1);

    % Find the lengths of each run of events
    runLengths = lastPlusOneRows - firstRows;
    
    % Only mark as a burst if the run length is at least minSpikesPerBurst 
    %   (criterion 2)
    toMark = find(runLengths >= minSpikesPerBurst);
    for iMark = 1:length(toMark)
        inBurst(firstRows(toMark(iMark)):lastPlusOneRows(toMark(iMark))) = 1;
    end
end

% Store all event information in eventInfoExpanded
% The first 10 columns are the same as eventInfo
eventInfoExpanded = zeros(nEvents, 14);
eventInfoExpanded(:, 1:11) = eventInfo;
eventInfoExpanded(:, 12) = inBurst;
eventInfoExpanded(:, 13) = inCrude;
eventInfoExpanded(:, 14) = RegionNo;

varargout{1} = absDiffBreaks;
varargout{2} = willJump;
varargout{3} = inCrude;
varargout{4} = RegionNo;
varargout{5} = inBurst;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%       First, use a stringent baseline difference criterion to identify
%       a large "crude burst region" of a given length
%       Next, use a looser baseline difference criterion to identify a 
%       "refined burst region."
%       Finally, use the maximum interstimulus interval and minimum
%       number of spikes per burst to identify the "burst"

%                   - 'MinDiffRefined': minimum baseline difference (data units)
%                                       for identifying a "refined burst region"
%                   must be a positive scalar
%                   default == 8
minDiffRefinedDefault = 8;      % Default minimum baseline difference (d. u.)
                                %   for identifying a refined burst region
                                % TODO: Why 8 pA?
addParameter(iP, 'MinDiffRefined', minDiffRefinedDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
minDiffRefined    = iP.Results.MinDiffRefined;
inRefined = zeros(nEvents, 1);  % whether the event is in a refined burst region

%                                   11 = whether in a "refined burst region"
%                                   12 = whether in a "crude burst region"
%                                   13 = burst region number
%                               {5} = inRefined
%                               {6} = inBurst
    varargout{5} = inRefined;
    varargout{6} = inBurst;
eventInfoExpanded(:, 11) = inRefined;
eventInfoExpanded(:, 12) = inCrude;
eventInfoExpanded(:, 13) = RegionNo;
varargout{5} = inRefined;
varargout{6} = inBurst;

%% Identify a "refined burst region" within each crude burst region
%   A "refined burst region" is the largest group of contiguous events
%   such that the absolute difference between the breakpoint baseline of the 
%   first/last event and of the following event is at least minDiffRefined

% In a vector called inRefined, 
%   mark all events in a refined burst region with 1
for iRegion = 1:nRegions        % for each burst region
    % Find the first index in the region that satisfies jump criterion
    firstJumpIdx = find(absDiffBreaks >= minDiffRefined & ...
                        RegionNo == iRegion, 1, 'first');

    % Find the last index in the region that satisfies jump criterion
    lastJumpIdx = find(absDiffBreaks >= minDiffRefined & ...
                        RegionNo == iRegion, 1, 'last');

    % Mark the events between these rows
    inRefined(firstJumpIdx:lastJumpIdx) = 1;
end

%   A "burst" is the largest group of contiguous events in each refined burst
%   region that satisfies both of:
%       (1) each inter-event interval to the next event is not more than maxIsi,
%               except for the last event
%       (2) total number of events is at least minSpikesPerBurst

    % Find the longest run of events
    [runLength, runNo] = max(lastPlusOneRows - firstRows);
    
    % Only make this a burst if the runLength is at least minSpikesPerBurst 
    %   (criterion 2)
    if runLength >= minSpikesPerBurst
        inBurst(firstRows(runNo):lastPlusOneRows(runNo)) = 1;
    end

%}

