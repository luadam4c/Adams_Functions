function varargout = parse_pleth_trace (vectors, siSeconds, varargin)
%% Parses pleth traces
% Usage: [parsedParams, parsedData] = parse_pleth_trace (vectors, siSeconds, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       spike2MatPath = '/media/shareX/Data_for_test_analysis/parse_pleth_trace/20200811T175223_000.mat';
%       spike2Table = parse_spike2_mat(spike2MatPath);
%       channelValues = spike2Table.channelValues;
%       channelNames = spike2Table.channelNames;
%       plethVec = channelValues{strcmp(channelNames, 'Pleth 2')};
%       siSeconds = spike2Table{strcmp(channelNames, 'Pleth 2'), 'siSeconds'};
%       timeVecs = create_time_vectors('TimeUnits', 's', 'SamplingIntervalSeconds', siSeconds, 'Vectors', plethVec);
%       [parsedParams, parsedData] = parse_pleth_trace(plethVec, siSeconds, 'TraceFileName', spike2MatPath);
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
%       siSeconds   - sampling interval in seconds
%                   must be a positive vector
%       varargin    - 'TraceFileName': Name of the corresponding trace file(s)
%                   must be empty, a character vector, a string array 
%                       or a cell array of character arrays
%                   default == extracted from the .atf file
%                   - 'FileBase': file base for output files
%                   must be a string scalar or a character vector
%                   default == match the trace file name
%                   - 'PlethWindowSeconds' - pleth window for computing
%                                               breathing rate in seconds
%                   must be a numeric scalar
%                   default = 5 seconds
%                   - Any other parameter-value pair for TODO
%
% Requires:
%       cd/compute_running_windows.m
%       cd/create_error_for_nargin.m
% TODO
% count_samples
% count_vectors
% create_labels_from_numbers
% create_time_vectors
% extract_fileparts
% force_column_vector
% force_column_cell
% isemptycell
% match_row_count
% match_format_vector_sets
%
% Used by:
%       cd/parse_spike2_mat.m

% File History:
% 2020-08-12 Modified from parse_laser_trace.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
traceFileNameDefault = '';      % set later
plethWindowSecondsDefault = 5;
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
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'vectors', ...                   % vectors
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'siSeconds', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TraceFileName', traceFileNameDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'FileBase', fileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PlethWindowSeconds', plethWindowSecondsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Read from the Input Parser
parse(iP, vectors, siSeconds, varargin{:});
traceFileName = iP.Results.TraceFileName;
fileBase = iP.Results.FileBase;
plethWindowSeconds = iP.Results.PlethWindowSeconds;

% Keep unmatched arguments for the TODO function
otherArguments = iP.Unmatched;

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
timeVecs = create_time_vectors(nSamples, 'TimeUnits', 's', ...
                    'SamplingIntervalSeconds', siSeconds, 'Vectors', vectors);

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

%% Do the job
% Compute running windows
endPointsWindows = compute_running_windows(timeVecs, plethWindowSeconds);

% Extract pieces
pieces = extract_subvectors(vectors, endPointsWindows);

%% Output results
% parsedParams

parsedData.endPointsWindows = endPointsWindows;

% Output variably
varargout{1} = parsedParams;
if nargout > 1
    varargout{2} = parsedData;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
