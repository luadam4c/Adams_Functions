function nVectors = count_vectors (vectors)
%% Counts the number of vectors whether given an array or a cell array
% Usage: nVectors = count_vectors (vectors)
% Explanation:
%       Uses either size(x, 2) for arrays and numel() for cell arrays
% Example(s):
%       nVectors = count_vectors(data)
% Outputs:
%       nVectors    - number of vectors
%                   specified as a nonnegative integer scalar
%
% Arguments:
%       vectors     - vectors to count
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%
% Requires:
%       cd/iscellnumeric.m
%
% Used by:
%       cd/compute_single_neuron_errors.m
%       cd/find_pulse_endpoints.m
%       cd/match_vector_counts.m

% File History:
% 2018-10-10 Created by Adam Lu
% 

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
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Read from the Input Parser
parse(iP, vectors);

%% Do the job
if iscell(vectors)
    nVectors = numel(vectors);
elseif isnumeric(vectors)    
    nVectors = size(vectors, 2);
else
    error('vectors is not the right type!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%