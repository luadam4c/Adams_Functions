function nSamples = count_samples (vectors)
%% Counts the number of samples whether given an array or a cell array
% Usage: nSamples = count_samples (vectors)
% Explanation:
%       Uses either size(x, 1) for arrays and 
%           cellfun(@length, x) for cell arrays
% Example(s):
%       nSamples = count_samples(data)
% Outputs:
%       nSamples    - number of samples for each vector
%                   specified as a column vector
%
% Arguments:
%       vectors     - vectors to count samples from
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%
% Requires:
%       cd/iscellnumeric.m
%
% Used by:    
%       cd/compute_single_neuron_errors.m
%       cd/force_column_cell.m
%       cd/parse_pulse.m
%       cd/parse_pulse_response.m

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

% Force vectors to be a column cell array
%   Note: this will make nSamples a column vector
vectors = force_column_cell(vectors);

%% Do the job
if iscell(vectors)
    % Count for each vector
    nSamples = cellfun(@length, vectors);
elseif isnumeric(vectors)
    % All the vectors have the same number of samples
    nSamplesScalar = size(vectors, 1);

    % Repeat for the number of vectors
    nVectors = size(vectors, 2);
    nSamples = ones(nVectors, 1) * nSamplesScalar;
else
    error('vectors is not the right type!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if iscell(vectors) && ~iscolumn(vectors)
    vectors = vectors(:);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%