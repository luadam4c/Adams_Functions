function varargout = parse_repetitive_pulses (vectors, siMs, varargin)
%% Parses repetitive pulses
% Usage: [parsedParams, parsedData] = parse_repetitive_pulses (vectors, siMs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       spike2MatPath = '/media/shareX/2019octoberR01/Pleth/Data_Not_Used/code_testing/test2AtNight_200Hz.mat';
%       spike2Table = parse_spike2_mat(spike2MatPath);
%       channelValues = spike2Table.channelValues;
%       channelNames = spike2Table.channelNames;
%       gasVec = channelValues{strcmp(channelNames, 'O2')};
%       siMs = spike2Table{strcmp(channelNames, 'O2'), 'siSeconds'} * 1000;
%       [parsedParams, parsedData] = parse_repetitive_pulses(gasVec, siMs, 'PulseDirection', 'downward', 'TraceFileName', spike2MatPath);
%
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
%
% Arguments:
%       vectors     - vectors containing many pulse responses
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%       siMs        - sampling interval in ms
%                   must be a positive vector
%       varargin    - 'TraceFileName': Name of the corresponding trace file(s)
%                   must be empty, a character vector, a string array 
%                       or a cell array of character arrays
%                   default == extracted from the .atf file
%                   - 'PulseTableSuffix': Suffix for the saved pulse table
%                   must be a character vector or a string scalar 
%                   default == '_pulses'
%                   - 'PulseShape': shape of expected pulse
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'square'        - square pulse
%                       'first-order'   - first-order responses
%                   default = 'square'
%                   - 'PulseDirection': direction of expected pulse
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'downward'  - downward peaks (e.g., hypoxia for O2)
%                       'upward'    - upward peaks (e.g., hypercapnia for CO2)
%                       'auto'      - no preference (whichever is largest)
%                   default = 'auto'
%                   - 'MinInterPulseIntervalMs' - minimum inter-pulse interval
%                                                   in ms
%                   must be a numeric scalar
%                   default = 0 ms
%                   - 'MinPulseAmplitude' - minimum pulse amplitude
%                   must be a numeric vector
%                   default = use 2 * compute_rms_error(vectors) by default
%                   - 'FileBase': file base for output files
%                   must be a string scalar or a character vector
%                   default == match the trace file name
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/compute_rms_error.m
%       cd/construct_and_check_fullpath.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/create_time_vectors.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/force_column_vector.m
%       cd/force_column_cell.m
%       cd/force_string_end.m
%       cd/isemptycell.m
%       cd/match_format_vectors.m
%       cd/match_row_count.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/parse_gas_trace.m
%       cd/parse_laser_trace.m

% File History:
% 2019-09-09 Created by Adam Lu
% 2019-09-10 Added 'PulseDirection' as an optional argument
% 2019-09-13 Now uses parse_repetitive_pulses.m
% 2019-09-13 Added 'PulseShape' as an optional argument
% 2019-09-14 Now applies detection of square pulses with faster algorithm
% 2019-09-15 Added 'FileBase' as an optional argument
% 2020-07-16 Added 'MinPulseAmplitude' as an optional argument
% TODO: Allow different pulse directions for different vectors
% 

%% Hard-coded parameters
MS_PER_S = 1000;
validPulseShapes = {'square', 'first-order'};
validPulseDirections = {'upward', 'downward', 'auto'};

%% Default values for optional arguments
traceFileNameDefault = '';      % set later
pulseTableSuffixDefault = '_pulses';
pulseShapeDefault = 'square';   % detects square pulses by default
pulseDirectionDefault = 'auto'; % automatically detect largest peak by default
minInterPulseIntervalMsDefault = 0;
minPulseAmplitudeDefault = [];  % use compute_rms_error.m by default
fileBaseDefault = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'vectors', ...                   % vectors
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'siMs', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TraceFileName', traceFileNameDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'PulseTableSuffix', pulseTableSuffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PulseShape', pulseShapeDefault, ...
    @(x) any(validatestring(x, validPulseShapes)));
addParameter(iP, 'PulseDirection', pulseDirectionDefault, ...
    @(x) any(validatestring(x, validPulseDirections)));
addParameter(iP, 'MinInterPulseIntervalMs', minInterPulseIntervalMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'MinPulseAmplitude', minPulseAmplitudeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'FileBase', fileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, vectors, siMs, varargin{:});
traceFileName = iP.Results.TraceFileName;
pulseTableSuffix = iP.Results.PulseTableSuffix;
pulseShape = validatestring(iP.Results.PulseShape, validPulseShapes);
pulseDirection = validatestring(iP.Results.PulseDirection, validPulseDirections);
minInterPulseIntervalMs = iP.Results.MinInterPulseIntervalMs;
minPulseAmplitude = iP.Results.MinPulseAmplitude;
fileBase = iP.Results.FileBase;

% Keep unmatched arguments for the TODO() function
% otherArguments = iP.Unmatched;
% otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% TODO: Condense to a function?
% Force vectors to be a column vectors
vectors = force_column_vector(vectors);

% Force vectors to be a column cell array of column vectors
vectors = force_column_cell(vectors);

% Count the number of samples
nSamples = count_samples(vectors);

% Count the number of vectors
nVectors = count_vectors(vectors);

% Create time vectors in seconds
timeVec = create_time_vectors(nSamples, 'TimeUnits', 's', ...
                                'SamplingIntervalMs', siMs);

% Match the format
timeVec = match_format_vector_sets(timeVec, vectors);

% Convert sampling interval to seconds
siSeconds = siMs / MS_PER_S;

% Match the row count
siSeconds = match_row_count(siSeconds, nVectors);

% Match the number of file names to the number of vectors
[traceFileName, vectors] = match_format_vector_sets(traceFileName, vectors);

% Make sure fileBase is in agreement with nVectors
if isempty(fileBase)
    if isemptycell(traceFileName)
        fileBase = create_labels_from_numbers(1:nVectors, ...
                                'Prefix', strcat(create_time_stamp, '_'));
    else
        % Extract from the trace file name
        fileBase = extract_fileparts(traceFileName, 'pathbase');

        % Force as a cell array
        fileBase = force_column_cell(fileBase);
    end
else
    % Force as a cell array
    fileBase = force_column_cell(fileBase);

    if numel(fileBase) ~= nVectors
        fprintf('Number of file bases must match the number of vectors!\n');
        varargout{1} = table.empty;
        varargout{2} = table.empty;
        return
    end
end

%% Parse the amplitude of any change
% Compute the extreme values
maxValue = extract_elements(vectors, 'max');
minValue = extract_elements(vectors, 'min');

% Compute the baseline and steady-state values
%   Note: this is either maxValue or minValue, whichever is closest to
%           the first value
firstValue = extract_elements(vectors, 'first');
[baseValue, steadyValue] = ...
    arrayfun(@(x, y, z) decide_on_baseline_steadystate(x, y, z, ...
                                                        pulseDirection), ...
            firstValue, maxValue, minValue);

% Compute the amplitudes
amplitude = steadyValue - baseValue;

% Compute maximum noise amplitudes if not provided
if isempty(minPulseAmplitude)
    minPulseAmplitude = 2 * compute_rms_error(vectors);
else
    minPulseAmplitude = match_format_vectors(minPulseAmplitude, amplitude);
end

%% Parse stimulus protocol
% Create pulse tables
[pulseTables, secEndPoints, pulseEndPoints, pulseTablePaths] = ...
    cellfun(@(a, b, c, d, e, f, g) ...
            create_pulse_tables(a, b, c, d, e, f, g, ...
                        pulseTableSuffix, pulseShape, ...
                        pulseDirection, minInterPulseIntervalMs), ...
            timeVec, vectors, num2cell(baseValue), ...
            num2cell(amplitude), num2cell(minPulseAmplitude), ...
            traceFileName, fileBase, 'UniformOutput', false);

%% Output results
if nVectors > 1
    % Put parameters in a table
    parsedParams = table(nSamples, siSeconds, maxValue, minValue, ...
                        baseValue, steadyValue, amplitude, pulseTablePaths);

    % Put data in a table
    parsedData = table(timeVec, sections, secEndPoints, ...
                        pulseEndPoints, pulseWindows);
else
    % Put parameters in a structure
    parsedParams.nSamples = nSamples;
    parsedParams.siSeconds = siSeconds;
    parsedParams.maxValue = maxValue;
    parsedParams.minValue = minValue;
    parsedParams.baseValue = baseValue;
    parsedParams.steadyValue = steadyValue;
    parsedParams.amplitude = amplitude;
    parsedParams.pulseTablePaths = pulseTablePaths;

    % Put data in a structure
    parsedData.timeVec = timeVec;
    parsedData.pulseTables = pulseTables;
    parsedData.secEndPoints = secEndPoints;
    parsedData.pulseEndPoints = pulseEndPoints;
end
% Output variably
varargout{1} = parsedParams;
if nargout > 1
    varargout{2} = parsedData;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [baseValue, steadyValue] = ...
                decide_on_baseline_steadystate (firstValue, maxValue, ...
                                                minValue, pulseDirection)

% Hard-coded parameters
maxPulseAmpRatio = 10;

% Compute the distance from the first value to the extreme values
distFirstToMax = abs(firstValue - maxValue);
distFirstToMin = abs(firstValue - minValue);

% Decide which values are baseline and steady state
if distFirstToMax / distFirstToMin > maxPulseAmpRatio
    % Baseline is close to minimum
    if strcmp(pulseDirection, 'downward')
        % Force going down
        baseValue = firstValue;
        steadyValue = minValue;
    else
        baseValue = mean([firstValue, minValue]);
        steadyValue = maxValue;
    end
elseif distFirstToMin / distFirstToMax > maxPulseAmpRatio
    % Baseline is close to maximum
    if strcmp(pulseDirection, 'upward')
        % Force going up
        baseValue = firstValue;
        steadyValue = maxValue;
    else
        % Use the average of firstValue and maxValue as baseline
        baseValue = mean([firstValue, maxValue]);
        steadyValue = minValue;
    end
else
    % Baseline is in between minimum and maximum
    baseValue = firstValue;

    % Steady state depends on the pulse direction
    switch pulseDirection
        case 'auto'
            % Use whatever is larger
            if distFirstToMax > distFirstToMin
                steadyValue = maxValue;
            else
                steadyValue = minValue;
            end
        case 'upward'
            steadyValue = maxValue;
        case 'downward'
            steadyValue = minValue;
        otherwise
            error('pulseDirection unrecognized!');
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pulseTable, secEndPoints, pulseEndPoints, ...
                    pulseTablePath, handles] = ...
                create_pulse_tables(timeVec, vector, baseValue, amplitude, ...
                            minPulseAmplitude, traceFileName, ...
                            fileBase, pulseTableSuffix, ...
                            pulseShape, pulseDirection, minInterPulseIntervalMs)
%% Create pulse tables
% TODO: Pull out as its own function

%% Hard-coded parameters
MS_PER_S = 1000;

%% Preparation
% Create a figure suffix
figSuffix = strcat(pulseTableSuffix, '_detection');

% Make sure the sheet base ends with pulseTableSuffix
pulseTableBase = force_string_end(fileBase, pulseTableSuffix);
figBase = force_string_end(fileBase, figSuffix);

% Create paths for the pulse table
pulseTablePath = force_string_end(pulseTableBase, '.csv');

% Compute the sampling interval in ms
siMs = (timeVec(2) - timeVec(1)) * MS_PER_S;

% Compute the minimum inter-pulse interval in samples
minInterPulseIntervalSamples = floor(minInterPulseIntervalMs / siMs);

%% Make sure a pulse exists
% A pulse exists only if the pulse amplitude is greater than a threshold
if abs(amplitude) > minPulseAmplitude
    pulseExists = true;
else
    pulseExists = false;

    % Display warning
    fprintf('Warning: There is either too much noise or no pulse in %s\n', ...
            traceFileName);

    % Initialize outputs
    secEndPoints = [];
    pulseEndPoints = [];
    startTime = [];
    endTime = [];
    duration = [];
    tracePath = {};
    pathExists = [];
    idxStart = [];
    idxEnd = [];
    pulseWindows = [];
end

%% Find pulse endpoints
if pulseExists
    switch pulseShape
    case 'square'
        disp('Finding endpoints for square pulses ...');
        pulseEndPoints = find_all_pulse_endpoints(vector, baseValue, ...
                                                amplitude, pulseDirection, ...
                                                minInterPulseIntervalSamples);
        secEndPoints = pulseEndPoints;
    case 'first-order'
        % TODO: Use minInterPulseIntervalSamples

        % Find approximate windows for first-order responses
        disp(['Finding approximate sections that ', ...
                'contain first-order responses ...']);
        secEndPoints = find_section_endpoints(vector, baseValue, ...
                                                amplitude, pulseDirection);

        % Remove sections with less than 4 data points
        %   Note: The fit() function requires 4 data points
        %           to fit a first order response
        nSamples = diff(secEndPoints, 1, 1) + 1;
        toRemove = nSamples <= 4;
        secEndPoints(:, toRemove) = [];

        % Extract subvectors
        sections = extract_subvectors(vector, 'EndPoints', secEndPoints);

        % Find first-order response end points relative to each section
        disp('Finding first-order response end points ...');
        [idxPulseStartsRel, idxPulseEndsRel] = ...
            find_pulse_response_endpoints(sections, siMs, ...
                                            'ResponseLengthMs', 0);

        % Find the relative first-order response end points
        % TODO: Use fit_first_order_response.m directly and get amplitude and tau info
        %       and put in pulse table
        % Put together as end points
        pulseEndPointsRel = transpose([idxPulseStartsRel, idxPulseEndsRel]);

        % Shift the relative first-order response end points
        %   by the section starting indices
        pulseEndPoints = pulseEndPointsRel + ...
                            repmat(secEndPoints(1, :), 2, 1) - 1;
    otherwise
        error('pulseShape unrecognized!');
    end

    % Convert to pulse windows
    pulseWindows = timeVec(pulseEndPoints);

    % Extract the index
    idxStart = transpose(pulseEndPoints(1, :));
    idxEnd = transpose(pulseEndPoints(2, :));

    % Extract the time
    startTime = transpose(pulseWindows(1, :));
    endTime = transpose(pulseWindows(2, :));

    % Compute the duration
    duration = endTime - startTime;

    % Force file name as a cell array
    traceFileName = force_column_cell(traceFileName);

    % Match number of file names with the number of startTime
    traceFileName = match_row_count(traceFileName, size(startTime, 1));

    % Determine if the trace file exists
    [tracePath, pathExists] = construct_and_check_fullpath(traceFileName);
end

%% Outputs
% Create a pulse table
pulseTable = table(startTime, endTime, duration, ...
                    tracePath, pathExists, idxStart, idxEnd);

% Write to spreadsheet files
writetable(pulseTable, pulseTablePath);

% Plot and save a figure that verifies the pulse window detection
handles = plot_pulse_detection(timeVec, vector, pulseWindows, figBase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function endPoints = find_all_pulse_endpoints (vector, baseValue, ...
                                                amplitude, pulseDirection, ...
                                                minInterPulseIntervalSamples)
%% Find all pulse endpoints
% TODO: Improve the performance with a better algorithm

% Standardize the pulse vector
[vecPos, ampPos] = ...
    standardize_pulse(vector, baseValue, amplitude, pulseDirection);

% Compute the amplitude threshold
thresAmp = ampPos / 2;

% Get all differences
diffVecPos = diff(vecPos);

% Compute the maximum positive jump
maxJump = max(diffVecPos);

% Compute the minimum negative fall
minFall = min(diffVecPos);

% Compute the thresholds for minimum jump and minimum fall
thresJump = maxJump * 1 / 4;
thresFall = minFall * 1 / 4;

% Find the indices before jumps and before falls 
indBeforeJumps = find(diffVecPos > thresJump & vecPos(1:end-1) < thresAmp);
indBeforeFalls = find(diffVecPos < thresFall & vecPos(1:end-1) > thresAmp);

% Compute the differences between the indices before jumps
diffIndBeforeJumps = diff(indBeforeJumps);

% Decide whether each index is at least minInterPulseIntervalSamples after
%   the previous index
isIdxStart = [true; diffIndBeforeJumps >= minInterPulseIntervalSamples];

% Get all starting indices
indStart = indBeforeJumps(isIdxStart);

% For each starting index after the first, 
%   find the most recent index before fall
if numel(indStart) > 1
    indEnd = nan(size(indStart));
    for iStart = 2:numel(indStart)
        iTemp = find(indBeforeFalls < indStart(iStart), 1, 'last');
        if isempty(iTemp)
            error('Code logic error!');
        else
            indEnd(iStart - 1, 1) = indBeforeFalls(iTemp);
        end
    end
end

% The last index before fall is the last ending index
indEnd(numel(indStart), 1) = indBeforeFalls(end);

% Put together as end points
endPoints = transpose([indStart, indEnd]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function endPoints = find_section_endpoints(vector, baseValue, ...
                                            amplitude, pulseDirection)


% Standardize the pulse vector
[vecPos, ampPos] = ...
    standardize_pulse(vector, baseValue, amplitude, pulseDirection);

% Rename thresholds
oneFourthAmp = ampPos / 4;
threeFourthsAmp = ampPos * 3 / 4;
oneTenthAmp = ampPos / 10;

% Find the first points that reaches 1/4 and 3/4 of the amplitude
idxRef = 0;
idx25Rel = find(vecPos > oneFourthAmp, 1, 'first');
idx75Rel = find(vecPos > threeFourthsAmp, 1, 'first');
endPoints = [];
while ~isempty(idx25Rel) && ~isempty(idx75Rel)
    % Compute the half rise time in samples
    halfRiseTimeSamples = idx75Rel - idx25Rel;

    % Compute the starting index (relative)
    idxStartRel = max(1, idx25Rel - halfRiseTimeSamples);

    % Compute the starting index (original)
    idxStart = idxRef + idxStartRel;

    % Remove up to idx75Rel and update the reference index
    vecPos(1:idx75Rel) = [];
    idxRef = idxRef + idx75Rel;

    % Find the first point that reaches 1/10 of the amplitude
    idxEndRel = find(vecPos < oneTenthAmp, 1, 'first');

    % If not found, we are done
    if isempty(idxEndRel)
        break
    end

    % Compute the ending index (original)
    idxEnd = idxRef + idxEndRel;

    % Put end points together
    endPointsThis = [idxStart; idxEnd];

    % Add to all end points
    endPoints = [endPoints, endPointsThis];

    % Remove up to idxEndRel and update the reference index
    vecPos(1:idxEndRel) = [];
    idxRef = idxRef + idxEndRel;

    % Find the first points that reaches 1/4 and 3/4 of the amplitude again
    idx25Rel = find(vecPos > oneFourthAmp, 1, 'first');
    idx75Rel = find(vecPos > threeFourthsAmp, 1, 'first');
end

% If not found, use the first and last point
if isempty(endPoints)
    endPoints = [1; numel(vecPos)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vecPos, ampPos] = standardize_pulse (vector, baseValue, ...
                                                amplitude, pulseDirection)
%% Standardizes a pulse vector by making the values of interest 
%   start at zero and positive

% Subtract the vector by the baseline value 
vecShifted = vector - baseValue;

% Make amplitude and vector positive
switch pulseDirection
    case 'auto'
        % Detect all pulses
        ampPos = abs(amplitude);
        vecPos = abs(vecShifted);
    case 'upward'
        % Detect only upward pulses
        ampPos = amplitude;
        vecPos = vecShifted;
    case 'downward'
        % Detect only downward pulses
        ampPos = -amplitude;
        vecPos = -vecShifted;
    otherwise
        error('pulseDirection unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = plot_pulse_detection (timeVec, gasVec, pulseWindows, figBase)
%% Plots gas pulse window detection

% Display message
disp('Plotting pulse detection ...');

% Create a new figure
fig = figure;

% Plot the original gas trace
p = plot(timeVec, gasVec);

% Hold on
hold on;

% Plot the window boundaries
if ~isempty(pulseWindows)
    shades = plot_window_boundaries(pulseWindows(:), ...
                                    'BoundaryType', 'verticalShade');
else
    shades = gobjects(1);
end

% Save figure if requested                            
if ~isempty(figBase)
    saveas(fig, figBase, 'png');
end

% Return handles
handles.fig = fig;
handles.p = p;
handles.shades = shades;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Note: The following assumes that there are only two possible values
candValue = [minValue, maxValue];
firstValue = extract_elements(vectors, 'first');
[~, iTemp1] = min(abs([firstValue, firstValue] - candValue), [], 2);
baseValue = candValue(iTemp1);

% Note: the following is very slow:
% Rename thresholds
oneFourthAmp = ampPos / 4;
threeFourthsAmp = ampPos * 3 / 4;

% Find the index before first points that reaches 1/4 and 3/4 of the amplitude
idxRef = 0;
idxFirstAbove = find(vecPos > oneFourthAmp, 1, 'first');
endPoints = [];
while ~isempty(idxFirstAbove)
    % Compute the relative starting index
    idxStartRel = idxFirstAbove - 1;

    % Compute the original starting index
    idxStart = idxRef + idxStartRel;

    % Remove up to idxStartRel and update the reference index
    vecPos(1:idxStartRel) = [];
    idxRef = idxRef + idxStartRel;

    % Find the index before first point that reaches below 3/4 of the amplitude
    idxFirstBelow = find(vecPos < threeFourthsAmp, 1, 'first');

    % If not found, we are done
    if isempty(idxFirstBelow)
        break
    end

    % Compute the relative ending index
    idxEndRel = idxFirstBelow - 1;

    % Compute the original ending index
    idxEnd = idxRef + idxEndRel;

    % Put end points together
    endPointsThis = [idxStart; idxEnd];

    % Decide what to do to the new end points
    if ~isempty(endPoints) && ...
            idxStart - endPoints(2, end) < minInterPulseIntervalSamples
        % The pulse is part of the previous pulse: replace the last end point
        endPoints(2, end) = idxEnd;
    else
        % The pulse is new: add to all previous endpoints
        endPoints = [endPoints, endPointsThis];
    end

    % Remove up to idxEndRel and update the reference index
    vecPos(1:idxEndRel) = [];
    idxRef = idxRef + idxEndRel;

    % Find the index before first point that reaches 1/4 of the amplitude again
    idxFirstAbove = find(vecPos > oneFourthAmp, 1, 'first');
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
