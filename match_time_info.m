function varargout = match_time_info (tVecs, varargin)
%% Match time vector(s) with sampling interval(s) and number(s) of samples
% Usage: [tVecs, siMs, nSamples] = match_time_info (tVecs, siMs, nSamples, varargin)
% Explanation:
%       TODO
% Example(s):
%       tVecs = match_time_info([], 0.1, 10)
%       [~, siMs, nSamples] = match_time_info(1:10)
% Outputs:
%       tVecs       - time vector(s)
%                   specified as a numeric array 
%                       or a cell array of numeric arrays
%       siMs        - sampling interval(s) in milliseconds
%                   specified as a positive vector
%       nSamples    - number of samples
%                   specified as a positive integer vector
% Arguments:
%       tVecs       - time vector(s)
%                   must be empty or 
%                       a numeric array or a cell array of numeric arrays
%       siMs        - (opt) sampling interval(s) in milliseconds
%                   must be empty or a positive vector
%       nSamples    - (opt) number of samples for each vector
%                   must be empty or a positive integer vector
%       varargin    - 'TimeUnits': output time units
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'min'   - minutes
%                       's'     - seconds
%                       'ms'    - milliseconds
%                       'us'    - microseconds
%                   default == 's'
%
% Requires:
%       cd/compute_sampling_interval.m
%       cd/create_error_for_nargin.m
%       cd/create_time_vectors.m
%       cd/ispositivevector.m
%       cd/ispositiveintegervector.m
%       cd/match_row_count.m
%       cd/match_vector_count.m
%
% Used by:
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/parse_current_family.m
%       cd/parse_multiunit.m

% File History:
% 2019-02-20 Created by Adam Lu
% 2019-10-15 Added 'TimeUnits' as optional arguments
% 

%% Hard-coded parameters
validTimeUnits = {'min', 's', 'ms', 'us'};

%% Default values for optional arguments
siMsDefault = [];
nSamplesDefault = [];
timeUnitsDefault = 'ms';           	% return time vector in ms by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'tVecs', ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['tVecs must be either empty, a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add optional inputs to the Input Parser
addOptional(iP, 'siMs', siMsDefault, ...
    @(x) assert(isempty(x) || ispositivevector(x), ...
                'siMs must be either empty or a positive vector!'));
addOptional(iP, 'nSamples', nSamplesDefault, ...
    @(x) assert(isempty(x) || ispositiveintegervector(x), ...
                'nSamples must be either empty or a positive integer vector!'));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TimeUnits', timeUnitsDefault, ...
    @(x) any(validatestring(x, validTimeUnits)));

% Read from the Input Parser
parse(iP, tVecs, varargin{:});
siMs = iP.Results.siMs;
nSamples = iP.Results.nSamples;
timeUnits = validatestring(iP.Results.TimeUnits, validTimeUnits);

%% Do the job
% Compute sampling interval(s) and create time vector(s)
if ~isempty(tVecs)
    % Compute sampling interval(s) in ms
    if isempty(siMs)
        siMs = compute_sampling_interval(tVecs);
    end

    % Compute number of samples
    if nargout >= 3 && isempty(nSamples)
        nSamples = count_samples(tVecs);
    end
elseif isempty(tVecs) && ~isempty(siMs) && ~isempty(nSamples)
    % Create time vector(s)
    tVecs = create_time_vectors(nSamples, 'SamplingIntervalMs', siMs, ...
                                    'TimeUnits', timeUnits);
elseif isempty(tVecs) && (isempty(siMs) || isempty(nSamples))
    error('One of tVecs or both siMs and nSamples must be provided!');
end

% Count the maximum number of entries
if nargout >= 3
    nEntries = max([count_vectors(tVecs), numel(siMs), numel(nSamples)]);
else
    nEntries = max([count_vectors(tVecs), numel(siMs)]);
end

% Match the number of entries
tVecs = match_vector_count(tVecs, nEntries);
siMs = match_row_count(siMs, nEntries);
if nargout >= 3
    nSamples = match_row_count(nSamples, nEntries);
end

%% Outputs
varargout{1} = tVecs;
varargout{2} = siMs;
if nargout >= 3
    varargout{3} = nSamples;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%