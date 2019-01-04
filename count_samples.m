function nSamples = count_samples (vectors, varargin)
%% Counts the number of samples whether given an array or a cell array
% Usage: nSamples = count_samples (vectors, varargin)
% Explanation:
%       Uses either numel(x) for vectors, 
%           size(x, 1) for arrays,
%           or cellfun(@numel, x) for cell arrays
% Example(s):
%       nSamples = count_samples(data)
% Outputs:
%       nSamples    - number of samples for each vector
%                   specified as a column vector 
%                       or a cell array of column vectors
%
% Arguments:
%       vectors     - vectors to count samples from
%                   Note: If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array
%       varargin    - 'ForceColumnOutput': whether to force output as a column
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/iscellnumericvector.m
%       cd/force_column_numeric.m
%       cd/match_row_count.m
%
% Used by:    
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/create_average_time_vector.m
%       cd/create_indices.m
%       cd/extract_columns.m
%       cd/extract_subvectors.m
%       cd/find_window_endpoints.m
%       cd/m3ha_xolotl_plot.m
%       cd/parse_pulse.m
%       cd/parse_pulse_response.m

% File History:
% 2018-10-10 Created by Adam Lu
% 2018-12-18 Now not longers forces nSamples to be a column vector
% 2019-01-03 Now accepts cell arrays of non-vector arrays
% 

%% Default values for optional arguments
forceColumnOutputDefault = true;  % force output as a column by default

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
addRequired(iP, 'vectors', ...                   % vectors
    @(x) assert(isnumeric(x) || iscell(x), ...
                'vectors must be either a numeric array or a cell array!'));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ForceColumnOutput', forceColumnOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vectors, varargin{:});
forceColumnOutput = iP.Results.ForceColumnOutput;

%% Do the job
% Decide based on input type
if iscellnumericvector(vectors)
    % Count the number of elements for each vector
    nSamples = cellfun(@numel, vectors);
elseif iscell(vectors)
    % Count the number of elements for each array
    nSamples = cellfun(@(x) count_samples(x, ...
                        'ForceColumnOutput', forceColumnOutput), ...
                        vectors, 'UniformOutput', false);
elseif isnumeric(vectors)
    if isvector(vectors)
        % Count the number of elements
        nSamples = numel(vectors);
    else
        % All the vectors have the same number of samples
        nSamplesScalar = size(vectors, 1);

        % Count the number of vectors
        nVectors = size(vectors, 2);

        % Repeat to make a column vector
        nSamples = match_row_count(nSamplesScalar, nVectors);
    end
else
    error('vectors is not the right type!');
end

% Force nSamples to be a column vector unless requested not to
if forceColumnOutput
    nSamples = force_column_numeric(nSamples);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if iscell(vectors) && ~iscolumn(vectors)
    vectors = vectors(:);
end

%       cd/force_column_cell.m
% If vectors is a cell array, force vectors to be a column cell array
%   Note: this will make nSamples a column vector
if iscell(vectors)
    vectors = force_column_cell(vectors);
end

nSamples = cellfun(@length, vectors);
nSamples = length(vectors);
nSamples = ones(nVectors, 1) * nSamplesScalar;

%                   must be a numeric array or a cell array of numeric vectors

@(x) assert(isnumeric(x) || iscellnumericvector(x), ...
            ['vectors must be either a numeric array', ...
                'or a cell array of numeric vectors!']));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%