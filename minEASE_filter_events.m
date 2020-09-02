function [eventInfoFiltered, eventClassFiltered, ...
                    isCheckedFiltered, interEventIntervals] = ...
                minEASE_filter_events (eventInfo, eventClass, isChecked, ...
                                nSamples, siMs, timeStart, varargin)
%% Filter events and compute inter-event intervals based on criteria
% Usage: [eventInfoFiltered, eventClassFiltered, ...
%                   isCheckedFiltered, interEventIntervals] = ...
%               minEASE_filter_events (eventInfo, eventClass, isChecked, ...
%                               nSamples, siMs, timeStart, varargin)
% Outputs:
%       eventInfoFiltered - filtered eventInfo
%                   specified as a numeric 2d array
%       eventClassFiltered - filtered eventClass
%                   specified as a numeric vector
%       isCheckedFiltered - filtered isChecked
%                   specified as a numeric vector
%       interEventIntervals - inter-event intervals
%                   specified as a numeric vector
%
% Arguments:
%       eventInfo   - information about events (cf. find_directional_events.m)
%                   must be a numeric 2d array
%       eventClass  - class of events (cf. minEASE_detect_gapfree_events.m)
%                   must be a numeric 2d array
%       isChecked   - whether the event is checked
%                   must be a numeric 2d array
%       nSamples    - number of total samples in data
%                   must be a positive integer scalar
%       siMs        - sampling interval in milliseconds
%                   must be a positive scalar
%       timeStart   - start of time in milliseconds
%                   must be a numeric scalar
%       varargin    - 'ClassesToInclude': classes of events to include
%                   must be a positive integer vector
%                   default == 1:3 (Type I~III PSCs)
%                   - 'TimeWindow': time window (in ms) to look for events
%                   must be a numeric vector
%                   default == [min(timeVector), max(timeVector)]
%                   - 'EventInfoTimeUnits':units used for time in eventInfo
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'ms'        - milliseconds
%                       'samples'   - samples
%                   default == 'samples'
%
% Used by:
%       cd/minEASE_filter_output.m

% File History:
%   2018-07-26 Modified from /home/Matlab/Marks_Functions/paula/Oct2017/loadAllTheData.m
%   2018-07-29 Renamed compute_IEIs.m -> minEASE_filter_events.m
%   2018-07-29 Now returns eventInfoFiltered, eventClassFiltered
%   2018-08-03 Added isCheckedFiltered
%   2018-08-12 Added 'EventInfoTimeUnits' and implement the case when it is 'ms'
%   TODO: Make it an option to use breakpoint indices instead of peak indices for IEI
% 

%% Hard-coded constants
validTimeUnits = {'samples', 'ms'};

%   Constants to be consistent with minEASE_detect_gapfree_events.m
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

%% Default values for optional arguments
classesToIncludeDefault = 1:3;
timeWindowDefault = [];
eventInfoTimeUnitsDefault = 'samples';
                            % eventInfo times is in samples by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 6
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'eventInfo', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addRequired(iP, 'eventClass', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addRequired(iP, 'isChecked', ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'2d'}));
addRequired(iP, 'nSamples', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
addRequired(iP, 'siMs', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'scalar'}));
addRequired(iP, 'timeStart', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ClassesToInclude', classesToIncludeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'vector'}));
addParameter(iP, 'TimeWindow', timeWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'EventInfoTimeUnits', eventInfoTimeUnitsDefault, ...
    @(x) any(validatestring(x, validTimeUnits)));

% Read from the Input Parser
parse(iP, eventInfo, eventClass, isChecked, ...
          nSamples, siMs, timeStart, varargin{:});
classesToInclude = iP.Results.ClassesToInclude;
timeWindow = iP.Results.TimeWindow;
eventInfoTimeUnits = validatestring(iP.Results.EventInfoTimeUnits, ...
                                    validTimeUnits);

%% Preparation
% Construct a time vector if needed
if strcmp(eventInfoTimeUnits, 'samples')
    % Construct a time vector
    timeVector = (timeStart - siMs) + (1:nSamples)' * siMs;
end

% Create a default time window if not provided
if isempty(timeWindow)
    timeWindow = [timeStart, timeStart + nSamples * siMs];
end

%% Filter the events
% If there are no events, return an empty vector
if isempty(eventInfo)
    eventInfoFiltered = [];
    eventClassFiltered = [];
    isCheckedFiltered = [];
    interEventIntervals = [];
    return;
end

% Get the minimum and maximum time to look for events
minTime = timeWindow(1);
maxTime = timeWindow(2);

% Extract the event times
if strcmp(eventInfoTimeUnits, 'samples')
    % Get the peak indices
    indPeaks = eventInfo(:, IDXPEAK_COLNUM);

    % Make sure indPeaks are integers
    if ~isempty(indPeaks)
        validateattributes(indPeaks, {'numeric'}, {'positive', 'integer'});
    end

    % Use the peak indices for the event times
    eventTimes = timeVector(indPeaks);
elseif strcmp(eventInfoTimeUnits, 'ms')
    % eventInfo already contains actual times in ms
    eventTimes = eventInfo(:, IDXPEAK_COLNUM);
else
    error('The time units %s is unrecognized!', eventInfoTimeUnits);
end

% Get the indices to include:
%   only events of specific classes and within specific time windows
indToInclude = ismember(eventClass, classesToInclude) & ...
                eventTimes >= minTime & eventTimes <= maxTime;

% Filter eventInfo, eventClass and isChecked accordingly
eventInfoFiltered = eventInfo(indToInclude, :);
eventClassFiltered = eventClass(indToInclude, :);
isCheckedFiltered = isChecked(indToInclude, :);

%% Compute the inter-event intervals
% If there are no more events, return an empty vector
if isempty(indToInclude)
    eventInfoFiltered = [];
    eventClassFiltered = [];
    isCheckedFiltered = [];
    interEventIntervals = [];
    return;
end

% Restrict the event times to the indices to include
eventTimesToInclude = eventTimes(indToInclude);

% Compute the inter-event intervals
interEventIntervals = diff(eventTimesToInclude);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%