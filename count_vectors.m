function nVectors = count_vectors (vectors, varargin)
%% Counts the number of vectors whether given an array or a cell array
% Usage: nVectors = count_vectors (vectors, varargin)
% Explanation:
%       Uses either 1 for vectors,
%               size(x, 2) for arrays 
%               or numel() for cell arrays
% Example(s):
%       nVectors = count_vectors(data)
% Outputs:
%       nVectors    - number of vectors
%                   specified as a nonnegative integer vector
%
% Arguments:
%       vectors     - vectors to count
%                   Note: If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array
%       varargin    - 'TreatArrayAsVector': whether to treat a non-vector array 
%                                           as a single vector
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/iscellnumericvector.m
%
% Used by:
%       cd/collapse_identical_vectors.m
%       cd/compute_all_pulse_responses.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/extract_channel.m
%       cd/extract_elements.m
%       cd/identify_repetitive_pulses.m
%       cd/parse_pulse_response.m
%       cd/plot_protocols.m
%       cd/plot_traces.m

% File History:
% 2018-10-10 Created by Adam Lu
% 2019-01-03 Now returns a vector if input is a cell array of non-vectors
% 

%% Default values for optional arguments
treatArrayAsVectorDefault = false;  % treat arrays as many vectors by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vectors', ...                   % vectors
    @(x) assert(isnumeric(x) || iscell(x), ...
                'vectors must be either a numeric array or a cell array!'));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TreatArrayAsVector', treatArrayAsVectorDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vectors, varargin{:});
treatArrayAsVector = iP.Results.TreatArrayAsVector;

%% Do the job
if iscell(vectors) && (treatArrayAsVector || iscellnumericvector(vectors))
    nVectors = numel(vectors);
elseif iscell(vectors)
    nVectors = cellfun(@count_vectors, vectors);
elseif isnumeric(vectors)    
    if isvector(vectors)
        nVectors = 1;
    else
        nVectors = size(vectors, 2);
    end
else
    error('vectors is not the right type!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors

    @(x) assert(isnumeric(x) || iscellnumericvector(x), ...
                ['vectors must be either a numeric array', ...
                    'or a cell array of numeric vectors!']));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%