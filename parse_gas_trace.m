function varargout = parse_gas_trace (vectors, siMs, varargin)
%% Parses gas traces
% Usage: [parsedParams, parsedData] = parse_gas_trace (vectors, siMs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       spike2MatPath = '/media/shareX/2019octoberR01/Pleth/Data/new_pleth_data/20190904T185956_Animal1_Norm-Hypo-Norm-HypoCO2-Norm_5-10-10-10-5.mat';
%       spike2Table = parse_spike2_mat(spike2MatPath);
%       channelValues = spike2Table.channelValues;
%       channelNames = spike2Table.channelNames;
%       gasVec = channelValues{strcmp(channelNames, 'O2')};
%       siMs = spike2Table{4, 'siSeconds'} * 1000;
%       [parsedParams, parsedData] = parse_gas_trace(gasVec, siMs, 'PulseDirection', 'downward', 'TraceFileName', spike2MatPath);
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
%                   must be empty, a characeter vector, a string array 
%                       or a cell array of character arrays
%                   default == extracted from the .atf file
%                   - 'PulseDirection': direction of expected pulse
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'downward'  - downward peaks (e.g., hypoxia for O2)
%                       'upward'    - upward peaks (e.g., hypercapnia for CO2)
%                       'auto'      - no preference (whichever is largest)
%                   default = 'auto'
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/construct_and_check_fullpath.m
%       cd/create_error_for_nargin.m
%       cd/create_time_vectors.m
%       cd/struct2arglist.m
% TODO:
% create_labels_from_numbers
% extract_elements
% extract_subvectors
% force_column_vector
% force_column_cell
% force_string_end
% count_samples
% isemptycell
% match_row_count
%
% Used by:
%       cd/parse_spike2_mat.m

% File History:
% 2019-09-09 Created by Adam Lu
% 2018-09-10 Added 'PulseDirection' as an optional argument
% TODO: Allow different pulse directions for different vectors
% 

%% Hard-coded parameters
MS_PER_S = 1000;
validPulseDirections = {'upward', 'downward', 'auto'};

% TODO: Make optional argument
fileBase = '';                  % file base for output files

%% Default values for optional arguments
traceFileNameDefault = '';      % set later
pulseDirectionDefault = 'auto';  % automatically detect largest peak by default

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
addParameter(iP, 'PulseDirection', pulseDirectionDefault, ...
    @(x) any(validatestring(x, validPulseDirections)));

% Read from the Input Parser
parse(iP, vectors, siMs, varargin{:});
traceFileName = iP.Results.TraceFileName;
pulseDirection = validatestring(iP.Results.PulseDirection, validPulseDirections);

% Keep unmatched arguments for the TODO() function
% otherArguments = iP.Unmatched;
% otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
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
    arrayfun(@(x, y, z) decide_on_baseline_steadystate(x, y, z, pulseDirection), ...
                    firstValue, maxValue, minValue);

% Compute the amplitudes
amplitude = steadyValue - baseValue;

%% Parse stimulus protocol
% Create pulse tables
[pulseTables, secEndPoints, pulseEndPoints, pulseTablePaths] = ...
    cellfun(@(x, y, z, u, v, w) ...
            create_pulse_tables(x, y, z, u, v, w, pulseDirection), ...
            timeVec, vectors, num2cell(baseValue), num2cell(amplitude), ...
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

function [pulseTable, secEndPoints, pulseEndPoints, pulseTablePath, handles] = ...
                create_pulse_tables(timeVec, vector, baseValue, amplitude, ...
                                    traceFileName, fileBase, pulseDirection)
%% Create pulse tables
% TODO: Pull out as its own function

% Hard-coded parameters
MS_PER_S = 1000;

% Note: Must be consistent with parse_iox.m
pulseTableSuffix = '_gas_pulses';
figSuffix = '_gas_pulse_detection';

% Make sure the sheet base ends with pulseTableSuffix
pulseTableBase = force_string_end(fileBase, pulseTableSuffix);
figBase = force_string_end(fileBase, figSuffix);

% Create paths for the pulse table
pulseTablePath = force_string_end(pulseTableBase, '.csv');

% Find approximate windows for first-order responses
disp('Finding approximate gas pulse sections ...');
secEndPoints = find_section_endpoints(vector, baseValue, ...
                                        amplitude, pulseDirection);

% Extract subvectors
sections = extract_subvectors(vector, 'EndPoints', secEndPoints);

% Compute the sampling interval in ms
siMs = (timeVec(2) - timeVec(1)) * MS_PER_S;

% Find the protocol end points
% TODO: Use fit_first_order_response.m directly and get amplitude and tau info
%       and put in pulse table
% Find pulse response end points relative to each section
disp('Finding gas pulse end points ...');
[idxPulseStartsRel, idxPulseEndsRel] = ...
    find_pulse_response_endpoints(sections, siMs, 'ResponseLengthMs', 0);

% Put together as end points
pulseEndPointsRel = transpose([idxPulseStartsRel, idxPulseEndsRel]);

% Fix the pulse end points
pulseEndPoints = pulseEndPointsRel + repmat(secEndPoints(1, :), 2, 1);

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

% Match traceFileName with the number of sections
traceFileName = match_format_vector_sets(traceFileName, sections);

% Determine if the trace file exists
[tracePath, pathExists] = construct_and_check_fullpath(traceFileName);

% Create a pulse table
pulseTable = table(startTime, endTime, duration, ...
                    tracePath, pathExists, idxStart, idxEnd);

% Write to spreadsheet files
writetable(pulseTable, pulseTablePath);

% Plot and save a figure that verifies the pulse window detection
handles = plot_gas_pulse_detection(timeVec, vector, pulseWindows, figBase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = plot_gas_pulse_detection (timeVec, gasVec, pulseWindows, figBase)
%% Plots gas pulse window detection

% Display message
disp('Plotting gas pulse detection ...');

% Create a new figure
fig = figure;

% Plot the original gas trace
p = plot(timeVec, gasVec);

% Hold on
hold on;

% Plot the window boundaries
shades = plot_window_boundaries(pulseWindows(:), ...
                                'BoundaryType', 'verticalShade');

% Save figure if requested                            
if ~isempty(figBase)
    saveas(fig, figBase, 'png');
end

% Return handles
handles.fig = fig;
handles.p = p;
handles.shades = shades;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function endPoints = find_section_endpoints(vector, baseValue, ...
                                            amplitude, pulseDirection)

% Subtract the vector by the baseline value 
vecShifted = vector - baseValue;

% Make amplitude and vector positive
switch pulseDirection
    case 'auto'
        % Detect all pulses
        absAmp = abs(amplitude);
        vecPos = abs(vecShifted);
    case 'upward'
        % Detect only upward pulses
        absAmp = amplitude;
        vecPos = vecShifted;
    case 'downward'
        % Detect only downward pulses
        absAmp = -amplitude;
        vecPos = -vecShifted;
    otherwise
        error('pulseDirection unrecognized!');
end

% Find the first points that reaches 1/4 and 3/4 of the amplitude
idxRef = 0;
idx25Rel = find(vecPos > absAmp / 4, 1, 'first');
idx75Rel = find(vecPos > absAmp * 3 / 4, 1, 'first');
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
    idxEndRel = find(vecPos < absAmp / 10, 1, 'first');

    % If not found, we are done
    if isempty(idxEndRel)
        break;
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
    idx25Rel = find(vecPos > absAmp / 4, 1, 'first');
    idx75Rel = find(vecPos > absAmp * 3 / 4, 1, 'first');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

candValue = [minValue, maxValue];
firstValue = extract_elements(vectors, 'first');
[~, iTemp1] = min(abs([firstValue, firstValue] - candValue), [], 2);
baseValue = candValue(iTemp1);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%