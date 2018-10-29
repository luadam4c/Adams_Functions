function nVectors = count_vectors (vectors)
%% Counts the number of vectors whether given an array or a cell array
% Usage: nVectors = count_vectors (vectors)
% Explanation:
%       Uses either 1 for vectors,
%               size(x, 2) for arrays 
%               or numel() for cell arrays
% Example(s):
%       nVectors = count_vectors(data)
% Outputs:
%       nVectors    - number of vectors
%                   specified as a nonnegative integer scalar
%
% Arguments:
%       vectors     - vectors to count
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%
% Requires:
%       cd/iscellnumericvector.m
%
% Used by:
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/find_pulse_endpoints.m
%       cd/force_column_cell.m
%       cd/plot_traces.m

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
    @(x) assert(isnumeric(x) || iscellnumericvector(x), ...
                ['vectors must be either a numeric array', ...
                    'or a cell array of numeric vectors!']));

% Read from the Input Parser
parse(iP, vectors);

%% Do the job
if iscell(vectors)
    nVectors = numel(vectors);
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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%