function parsedParams = parse_peaks (vec, varargin)
%% Parses peaks for one vector
% Usage: parsedParams = parse_peaks (vec, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load sunspot.dat
%       year = sunspot(:, 1);
%       avSpots = sunspot(:, 2);
%       parsedParams = parse_peaks(avSpots)
%       parsedParams = parse_peaks(avSpots, 'ParseMode', 'max')
%       parsedParams = parse_peaks(avSpots, 'ParseMode', 'first')
%
% Outputs:
%       parsedParams    - parsed scalars
%                       specified as a table or struct array
%
% Arguments:
%       vec         - vector to parse
%                   must be a numeric vector
%       varargin    - 'ParseMode': parse mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'all'   - all peaks
%                       'max'   - the maximum peak
%                       'maxOfAll'  - the maximum peak, considering all peaks
%                       'first' - the first peak
%                   default == 'all'
%                   - 'OutputMode': output mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'table' - table
%                       'struct'- structure array
%                       'auto'  - 'table' if parseMode is 'all';
%                                   'struct' otherwise
%                   default == 'auto'
%                   - 'PeakLowerBound': peak lower bound
%                   must be a numeric scalar
%                   default == none
%                   - Any other parameter-value pair for findpeaks()
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/find_custom.m
%       cd/find_troughs_from_peaks.m
%       cd/force_column_vector.m
%       cd/isnumericvector.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/m3ha_plot_simulated_traces.m

% File History:
% 2020-04-20 Created by Adam Lu
% 2020-04-22 Added 'maxOfAll' as a valid parseMode

%% Hard-coded parameters
validParseModes = {'all', 'max', 'maxOfAll', 'first'};
validOutputModes = {'auto', 'table', 'struct'};
parsedParamsVariableNames = {'peakNum'; 'idxPeak'; 'peakAmp'; ...
                            'peakWidth'; 'peakProm'; ...
                            'idxPeakStart'; 'idxPeakEnd'};
parsedParamsVariableTypes = {'double'; 'double'; 'double'; ...
                            'double'; 'double'; 'double'; 'double'};

%% Default values for optional arguments
parseModeDefault = 'all';           % set later
outputModeDefault = 'auto';         % set later
peakLowerBoundDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'vec', ...
    @(x) assert(isnumericvector(x), ...
                'vecs must be a numeric vector!'));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ParseMode', parseModeDefault, ...
    @(x) any(validatestring(x, validParseModes)));
addParameter(iP, 'OutputMode', outputModeDefault, ...
    @(x) any(validatestring(x, validOutputModes)));
addParameter(iP, 'PeakLowerBound', peakLowerBoundDefault, ...
    @(x) assert(isnumericvector(x), ...
                'vecs must be a numeric vector!'));

% Read from the Input Parser
parse(iP, vec, varargin{:});
parseMode = validatestring(iP.Results.ParseMode, validParseModes);
outputMode = validatestring(iP.Results.OutputMode, validOutputModes);
peakLowerBound = iP.Results.PeakLowerBound;

% Keep unmatched arguments for the findpeaks() function
otherArguments = iP.Unmatched;

% Update 'NPeaks' if necessary
switch parseMode
    case {'max', 'first'}
        otherArguments.NPeaks = 1;
    otherwise
        % Do nothing
end

% Update 'MinPeakHeight' if necessary
if strcmp(parseMode, 'max')
    % Find the second highest value
    secondHighestValue = find_second_highest_value(vec);

    % Set as minimum peak height
    if ~isnan(secondHighestValue)
        otherArguments.MinPeakHeight = secondHighestValue;
    end
end

% Update 'SortStr' if necessary
if strcmp(parseMode, 'maxOfAll')
    otherArguments.SortStr = 'descend';
end

% Transform to a cell array
otherArguments = struct2arglist(otherArguments);

%% Preparation
% Decide on the output mode
if strcmp(outputMode, 'auto')
    switch parseMode
        case 'all'
            outputMode = 'table';
        case {'max', 'maxOfAll', 'first'}
            outputMode = 'struct';
        otherwise
            error('parseMode unrecognized!');
    end
end

% Initialize output
switch outputMode
     case 'table'
        % TODO: create_empty_table.m
        parsedParams = table('Size', [0, numel(parsedParamsVariableNames)], ...
                            'VariableNames', parsedParamsVariableNames, ...
                            'VariableTypes', parsedParamsVariableTypes);
     case 'struct'
    parsedParams.peakNum = NaN;
    parsedParams.idxPeak = NaN;
    parsedParams.peakAmp = NaN;
    parsedParams.peakWidth = NaN;
    parsedParams.peakProm = NaN;
    parsedParams.idxPeakStart = NaN;
    parsedParams.idxPeakEnd = NaN;
     otherwise
         % Do nothing
end

%% Find all peaks
% Count the number of samples
nSamples = numel(vec);

% Return if empty
if nSamples == 0
    return
end

% Find all peaks
[peakAmp, idxPeak, peakWidth, peakProm] = findpeaks(vec, otherArguments{:});

% Return if empty
if isempty(peakAmp)
    return
end

% Force as column vectors
[peakAmp, idxPeak, peakWidth, peakProm] = ...
    argfun(@force_column_vector, peakAmp, idxPeak, peakWidth, peakProm);

% Count the number of peaks
nPeaks = numel(peakAmp);

% Create peak numbers
peakNum = transpose(1:nPeaks);

% Place in a table
peakTable = table(peakNum, idxPeak, peakAmp, peakWidth, peakProm);

%% Find peak regions
% Sort the peaks
[idxPeakSorted, origInd] = sort(idxPeak, 'ascend');
[~, sortedInd] = sort(origInd);

% Find the minimums between each peak and the bounds
[~, indTroughsSorted] = ...
    find_troughs_from_peaks(vec, [1; idxPeakSorted; nSamples]);

% Set the peak starts to be the trough on the left
idxPeakStartSorted = indTroughsSorted(1:end-1);

% Set the peak ends to be the trough on the right
idxPeakEndSorted = indTroughsSorted(2:end);

% Make sure peak start and end values are at least lower bound
if ~isempty(peakLowerBound)
    [idxPeakStartSorted, idxPeakEndSorted] = ...
        arrayfun(@(x, y) modify_peak_endpoints(vec, x, y, peakLowerBound), ...
                idxPeakStartSorted, idxPeakEndSorted);
end

% Reorder to original
idxPeakStart = idxPeakStartSorted(sortedInd);
idxPeakEnd = idxPeakEndSorted(sortedInd);

% Add to table
peakTable = addvars(peakTable, idxPeakStart, idxPeakEnd);

%% Output results
% Restrict to the maximum peak
if strcmp(parseMode, 'maxOfAll')
    peakTable = peakTable(1, :);
end

% Reorganize outputs if necessary
switch outputMode
    case 'table'
        parsedParams = peakTable;
    case 'struct'
        parsedParams = table2struct(peakTable, 'ToScalar', false);
    otherwise
        error('outputMode unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function secondHighestValue = find_second_highest_value (vec)
% TODO: Pull out to its own function

% Find the maximum value
maxValue = max(vec);

% Make all maximum values NaN
vec(vec == maxValue) = NaN;

% Find the next highest value
secondHighestValue = max(vec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idxPeakStart, idxPeakEnd] = ...
                modify_peak_endpoints(vec, idxPeakStart, idxPeakEnd, ...
                                        peakLowerBound)
%% Modify peak endpoints

% Restrict vector to old peak region
vecPeakOld = vec(idxPeakStart:idxPeakEnd);
idxVecPeakOldStart = idxPeakStart;

% Find the minimum peak region value
minPeakRegionValue = min(vecPeakOld);

% If the minimum peak region value is at least the lower bound, we are done
if minPeakRegionValue >= peakLowerBound
    return
end

% Find the peak index relative to peak region
[maxPeakRegionValue, idxRelPeak] = max(vecPeakOld);

% If the peak value is less than the lower bound, remove the peak
if maxPeakRegionValue < peakLowerBound
    idxPeakStart = NaN;
    idxPeakEnd = NaN;
    return
end

% Restrict old peak region to left and right parts
vecPeakLeftOld = vecPeakOld(1:idxRelPeak);
idxVecPeakLeftOldStart = idxVecPeakOldStart;
vecPeakRightOld = vecPeakOld(idxRelPeak:end);
idxVecPeakRightOldStart = idxVecPeakOldStart - 1 + idxRelPeak;

% Update the index peak start if necessary
if min(vecPeakLeftOld) < peakLowerBound
    % Find the relative index right before new peak start
    idxRelPeakBeforeStart = find(vecPeakLeftOld < peakLowerBound, 1, 'last');

    % Update the index peak start
    idxPeakStart = (idxVecPeakLeftOldStart - 1) + idxRelPeakBeforeStart + 1;
end

% Update the index peak end if necessary
if min(vecPeakRightOld) < peakLowerBound
    % Find the relative index of new peak end
    idxRelPeakAfterEnd = find(vecPeakRightOld < peakLowerBound, 1, 'first');

    % Update the index peak end
    idxPeakEnd = (idxVecPeakRightOldStart - 1) + idxRelPeakAfterEnd - 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%