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
%       [parsedParams, parsedData] = parse_gas_trace(gasVec, siMs);
%       timeVec = parsedData.timeVec;
%       pulseWindows = parsedData.pulseWindows;
%       pulseWindowBoundaries = pulseWindows(:)
%       plot(timeVec, gasVec); hold on;
%       plot_window_boundaries(pulseWindowBoundaries, 'BoundaryType', 'verticalShade')
%
%       parse_gas_trace(vectors, siMs, 'FileBases', 
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
%       varargin    - 'FileBases': Base name of the corresponding trace file(s)
%                   must be empty, a characeter vector, a string array 
%                       or a cell array of character arrays
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/create_time_vectors.m
%       cd/struct2arglist.m
% TODO:
% create_labels_from_numbers
% extract_elements
% extract_subvectors
% force_column_vector
% force_column_cell
% count_samples
% match_row_count
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-09 Created by Adam Lu
% 

%% Hard-coded parameters
MS_PER_S = 1000;

%% Default values for optional arguments
fileBasesDefault = '';      % set later

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
addParameter(iP, 'FileBases', fileBasesDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));

% Read from the Input Parser
parse(iP, vectors, siMs, varargin{:});
fileBases = iP.Results.FileBases;

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

% Make sure fileBases is in agreement with nVectors
if isempty(fileBases)
    fileBases = create_labels_from_numbers(1:nVectors, ...
                                'Prefix', strcat(create_time_stamp, '_'));
else
    fileBases = force_column_cell(fileBases);
    
    if numel(fileBases) ~= nVectors
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

% Compute the baseline values
%   Note: this is either maxValue or minValue, whichever is closest to
%           the first value
candValue = [minValue, maxValue];
firstValue = extract_elements(vectors, 'first');
[~, iTemp1] = min(abs([firstValue, firstValue] - candValue), [], 2);
baseValue = candValue(iTemp1);

% Compute the steady-state values
steadyValue = candValue(3 - iTemp1);

% Compute the amplitudes
amplitude = steadyValue - baseValue;

%% Parse stimulus protocol
%% TODO: Pull out as a function
% Hard-coded parameters
pulseTableSuffix = 'gas_pulses';

% Create paths for the pulse table
pulseTablePaths = strcat(fileBases, '_', pulseTableSuffix, '.csv');

% Find approximate windows for first-order responses
secEndPoints = cellfun(@(x, y, z) find_section_endpoints(x, y, z), ...
                        vectors, num2cell(baseValue), num2cell(amplitude), ...
                        'UniformOutput', false);

% Extract subvectors
sections = cellfun(@(x, y) extract_subvectors(x, 'EndPoints', y), ...
                        vectors, secEndPoints, 'UniformOutput', false);

% Find the protocol end points
% TODO: Use fit_first_order_response.m directly and get amplitude and tau info
%       and put in stim table
pulseEndPointsRel = ...
    cellfun(@(x) find_pulse_endpoints_from_response(x, siMs), ...
            sections, 'UniformOutput', false);

% Fix the pulse end points
pulseEndPoints = cellfun(@(x, y) x + repmat(y(1, :), 2, 1), ...
                    pulseEndPointsRel, secEndPoints, 'UniformOutput', false);

% Convert to pulse windows
pulseWindows = cellfun(@(x, y) x(y), ...
                    timeVec, pulseEndPoints, 'UniformOutput', false);

% Create stim tables
pulseTables = cellfun(@(x, y) separate_and_create_table(x, y), ...
                    pulseEndPoints, pulseWindows, 'UniformOutput', false);

% Write to spreadsheet files
cellfun(@(x, y) writetable(x, y), pulseTables, pulseTablePaths);

%% Output results
if nVectors > 1
    % Put parameters in a table
    parsedParams = table(nSamples, siSeconds, maxValue, minValue, ...
                        baseValue, steadyValue, amplitude);

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

    % Put data in a structure
    parsedData.timeVec = timeVec;
    parsedData.sections = sections;
    parsedData.secEndPoints = secEndPoints;
    parsedData.pulseEndPoints = pulseEndPoints;
    parsedData.pulseWindows = pulseWindows;
end
% Output variably
varargout{1} = parsedParams;
if nargout > 1
    varargout{2} = parsedData;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function endPoints = find_section_endpoints(vector, baseValue, amplitude)

% Subtract the vector by the baseline value 
vecShifted = vector - baseValue;

% Make amplitude and vector positive
absAmp = abs(amplitude);
vecPos = abs(vecShifted);

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

function endPoints = find_pulse_endpoints_from_response (vec, siMs)
% Finds pulse endpoints from shape of pulse response

% Find pulse response end points
[idxPulseStartsRel, idxPulseEndsRel] = ...
    find_pulse_response_endpoints(vec, siMs, 'ResponseLengthMs', 0);

% Put together as end points
endPoints = transpose([idxPulseStartsRel, idxPulseEndsRel]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pulseTable = separate_and_create_table (pulseEndPoints, pulseWindows)
%% Create a pulse table

% Extract the index
idxStart = transpose(pulseEndPoints(1, :));
idxEnd = transpose(pulseEndPoints(2, :));

% Extract the time
timeStart = transpose(pulseWindows(1, :));
timeEnd = transpose(pulseWindows(2, :));

% Create a table
pulseTable = table(idxStart, idxEnd, timeStart, timeEnd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%