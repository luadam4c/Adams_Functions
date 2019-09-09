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
%       timeVecs = parsedData{1, 'timeVecs'};
%       timeVec = timeVecs{1};
%       pulseWindows = parsedData{1, 'pulseWindows'};
%       pulseWindowBoundaries = pulseWindows{1}(:)
%       plot(timeVec, gasVec); hold on;
%       plot_window_boundaries(pulseWindowBoundaries, 'BoundaryType', 'verticalShade')
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
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/create_time_vectors.m
%       cd/struct2arglist.m
% extract_elements
% extract_subvectors
% force_column_vector
% force_column_cell
% count_samples
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-09 Created by Adam Lu
% 

%% Hard-coded parameters
MS_PER_S = 1000;

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
% addParameter(iP, 'param1', param1Default, ...
    % TODO: validation function %);

% Read from the Input Parser
parse(iP, vectors, siMs, varargin{:});
% Force vectors to be a column cell array
vectors = force_column_cell(vectors);
% param1 = iP.Results.param1;

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
timeVecs = ...
    create_time_vectors(nSamples, 'SamplingIntervalMs', siMs, 'TimeUnits', 's');

% Match the format
timeVecs = match_format_vector_sets(timeVecs, vectors);

% Convert sampling interval to seconds
siSeconds = siMs / MS_PER_S;

% Match the row count
siSeconds = match_row_count(siSeconds, nVectors);

%% Do the job
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
                    timeVecs, pulseEndPoints, 'UniformOutput', false);

% Create stim tables
stimTables = cellfun()

%% Output results
% Put parameters in a table
parsedParams = table(nSamples, siSeconds, maxValue, minValue, ...
                    baseValue, steadyValue, amplitude);

% Put data in a table
parsedData = table(timeVecs, sections, secEndPoints, ...
                    pulseEndPoints, pulseWindows);

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

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%