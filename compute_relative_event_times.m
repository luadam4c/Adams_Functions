function [relEventTimes, relativeTimeWindow, origInd] = ...
                compute_relative_event_times (eventTimes, stimTimes, varargin)
%% Computes the relative event times from event times and stimulus times
% Usage: [relEventTimes, relativeTimeWindow, origInd] = ...
%               compute_relative_event_times (eventTimes, stimTimes, varargin)
% Explanation:
%       Computes the relative event times around stimulus times
%       For each stimulus time, a window around it is used
%       Note: Output is:
%               - a numeric vector if eventTimes is a vector 
%                       and stimTimes is a scalar
%               - a cell array of numeric vectors if eventTimes is a vector 
%                       and stimTimes is a vector
%               - a cell array of cell arrays of numeric vectors 
%                       if either eventTimes or stimTimes is a cell array
%
% Example(s):
%       eventTimes = randi(100, 100, 1);
%       stimTimes = 50;
%       relEventTimes = compute_relative_event_times(eventTimes, stimTimes)
%       eventTimes = randi(100, 100, 1);
%       stimTimes = 10:10:80;
%       relEventTimes = compute_relative_event_times(eventTimes, stimTimes)
%       eventTimes = {randi(100, 100, 1); randi(100, 100, 1) + 100};
%       stimTimes = {50; 150};
%       relEventTimes = compute_relative_event_times(eventTimes, stimTimes)
%       eventTimes = {randi(100, 100, 1); randi(100, 100, 1) + 100};
%       stimTimes = {10:10:80; 110:10:200};
%       relEventTimes = compute_relative_event_times(eventTimes, stimTimes)
%
% Outputs:
%       relEventTimes   - relative event times
%                       specified as a numeric vector,
%                           a cell array of numeric vectors or
%                           a cell array of cell arrays of numeric vectors
%       relativeTimeWindow - relative time window used
%                       specified as a numeric vector
%       origInd         - original indices of the extracted relative event times
%                       specified as a numeric vector,
%                           a cell array of numeric vectors or
%                           a cell array of cell arrays of numeric vectors
%
% Arguments:
%       eventTimes  - event times
%                       Note: If a cell array, each cell 
%                               will be matched with stimTimes
%                   must be a numeric vector or a cell array of numeric vectors
%       stimTimes   - stimulus times
%                       Note: If a cell array, each cell 
%                               will be matched with eventTimes
%                   must be a numeric vector or a cell array of numeric vectors
%       varargin    - 'RelativeTimeWindow': relative time window
%                   must be a 2-element numeric vector
%                   default == interStimInterval * 0.5 * [-1, 1]
%                   - 'StimIndices': stimulation indices to restrict to
%                   must be a positive integer array or string recognized by
%                        the 'Pattern' option of extract_subvectors.m
%                           'odd'   - odd indices
%                           'even'  - even indices
%                   default == no restrictions
%                   - 'ForceMatrixOutput': whether to force output as a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/argfun.m
%       cd/compare_events_pre_post_stim.m
%       cd/create_error_for_nargin.m
%       cd/create_shifted_vectors.m
%       cd/extract_subvectors.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
%       cd/force_row_vector.m
%       cd/iscellnumeric.m
%       cd/match_format_vector_sets.m
%       cd/plot_relative_events.m
%
% Used by:
%       cd/compute_psth.m

% File History:
% 2019-09-11 Moved from compute_psth.m
% 2019-09-15 Now returns relative event time window used
% 2019-10-10 Added 'StimIndices' as an optional argument
% 2019-10-10 Added 'ForceMatrixOutput' as an optional argument
% 2020-06-26 Added 'origInd' as an output
% TODO: Add option to shift relative event times by stimDelay
% 

%% Hard-coded parameters

%% Default values for optional arguments
stimIndicesDefault = [];        % take all stims by default
relativeTimeWindowDefault = [];
forceMatrixOutputDefault = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'eventTimes', ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['eventTimes must be either empty or a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'stimTimes', ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['stimTimes must be either empty or a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'RelativeTimeWindow', relativeTimeWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'StimIndices', stimIndicesDefault, ...
    @(x) validateattributes(x, {'numeric', 'char'}, {'2d'}));
addParameter(iP, 'ForceMatrixOutput', forceMatrixOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, eventTimes, stimTimes, varargin{:});
relativeTimeWindow = iP.Results.RelativeTimeWindow;
stimIndices = iP.Results.StimIndices;
forceMatrixOutput = iP.Results.ForceMatrixOutput;

%% Preparation
% Force as column vectors
[eventTimes, stimTimes] = ...
    argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', false), ...
            eventTimes, stimTimes);

% Make sure the vector numbers are identical and force as column cell arrays
[eventTimesCell, stimTimesCell] = ...
    match_format_vector_sets(eventTimes, stimTimes, 'ForceCellOutputs', true);

% Sort the times in ascending order
[eventTimesCell, stimTimesCell] = ...
    argfun(@(x) cellfun(@sort, x, 'UniformOutput', false), ...
            eventTimesCell, stimTimesCell);

% Restrict to certain stimulation windows if requested
if isempty(stimIndices)
    % Do nothing
elseif isnumeric(stimIndices)
    stimTimesCell = extract_subvectors(stimTimesCell, 'Indices', stimIndices);
elseif ischar(stimIndices)
    stimTimesCell = extract_subvectors(stimTimesCell, 'Pattern', stimIndices);
end

% Compute the default relative time window
if isempty(relativeTimeWindow)
    % Compute the average inter-stimulus interval
    interStimInterval = compute_average_interval(stimTimesCell);

    if isnan(interStimInterval)
        % Compute the minimum interval that contains 
        %   all event times for each set
        minHalfWidthEachCell = ...
            cellfun(@(x, y) max(abs([max(x) - y, min(x) - y])), ...
                    eventTimesCell, stimTimesCell);

        % Use an interval that works for all sets
        windowHalfWidth = max(minHalfWidthEachCell);
    else
        % Use half of the average inter-stimulus interval on each side
        windowHalfWidth = interStimInterval * 0.5;
    end

    relativeTimeWindow = windowHalfWidth * [-1, 1];
end

%% Do the job
% Extract relative event times for each window
[relEventTimesCellCell, origIndCellCell] = ...
    cellfun(@(x, y) compute_relative_event_times_helper(x, y, ...
                                                    relativeTimeWindow), ...
            eventTimesCell, stimTimesCell, 'UniformOutput', false);

%% Output results
if ~iscell(eventTimes) && ~iscell(stimTimes)
    % Extract from the cell array if there is only one set of event times
    %   and one set of stim times
    if numel(relEventTimesCellCell) == 1
        relEventTimesCell = relEventTimesCellCell{1};
        origIndCell = origIndCellCell{1};

        if numel(relEventTimesCell) == 1
            relEventTimes = relEventTimesCell{1};
            origInd = origIndCell{1};
        else
            relEventTimes = relEventTimesCell;
            origInd = origIndCell;
        end
    else
        error('Not implemented yet!')
    end
else
    if forceMatrixOutput
        % Put the event time arrays in a cell matrix
        %   Note: Each column is a file
        %         Each row is a stim
        [relEventTimes, origInd] = ...
            argfun(@(x) force_matrix(x, 'TreatCellNumAsArray', true), ...
                    relEventTimesCellCell, origIndCellCell);
    else
        relEventTimes = relEventTimesCellCell;
        origInd = origIndCellCell;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function meanInterval = compute_average_interval(eventTimes)
%% Compute the average inter-event interval
% TODO: Pull out as its own function

% Construct a function that calculates the mean over a matrix
nanmeansq = @(x) nanmean(nanmean(x));

% Compute the average inter-event interval
if iscell(eventTimes)
    meanInterval = nanmeansq(cellfun(@(x) nanmeansq(diff(x)), eventTimes));
else
    meanInterval = nanmeansq(diff(eventTimes));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [relEventTimes, origInd] = ...
                compute_relative_event_times_helper(eventTimes, stimTimes, ...
                                                            relativeTimeWindow)
%% Computes the relative event times from event times and stimulus times

% Extract a time window for each stimulus time
windows = create_shifted_vectors(relativeTimeWindow, stimTimes);

% Extract the event times corresponding to each time window
[eventTimesEachWindow, origInd] = ...
    extract_subvectors(eventTimes, 'Windows', windows, 'ForceCellOutput', true);

% Compute the relative event times for each time window
relEventTimes = cellfun(@(x, y) subtract(x, y), eventTimesEachWindow, ...
                        num2cell(stimTimes), 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function difference = subtract(x, y)

if isempty(x) || isempty(y)
    difference = x;
else
    difference = x - y;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vecsNew = create_shifted_vectors (vecOrig, shiftValues)
%% Creates shifted vectors based on values to shift
% TODO: Pull out as its own function

% Force the original vector as a column vector
vecOrig = force_column_vector(vecOrig);

% Force the values to shift as a row vector
shiftValues = force_row_vector(shiftValues);

% Create shifted vectors
%   Note: Each column corresponds to a time window
vecsNew = repmat(shiftValues, size(vecOrig)) + ...
            repmat(vecOrig, size(shiftValues));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
