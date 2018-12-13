function nSamples = count_samples (vectors)
%% Counts the number of samples whether given an array or a cell array
% Usage: nSamples = count_samples (vectors)
% Explanation:
%       Uses either length(x) for vectors, 
%           size(x, 1) for arrays,
%           or cellfun(@length, x) for cell arrays
% Example(s):
%       nSamples = count_samples(data)
% Outputs:
%       nSamples    - number of samples for each vector
%                   specified as a column vector
%
% Arguments:
%       vectors     - vectors to count samples from
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%
% Requires:
%       cd/iscellnumericvector.m
%       cd/force_column_numeric.m
%
% Used by:    
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/m3ha_xolotl_plot.m
%       cd/parse_pulse.m
%       cd/parse_pulse_response.m

% File History:
% 2018-10-10 Created by Adam Lu
% 2018-10-27 
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
    % Count for each vector
    nSamples = cellfun(@length, vectors);
elseif isnumeric(vectors)
    if isvector(vectors)
        nSamples = length(vectors);
    else
        % All the vectors have the same number of samples
        nSamplesScalar = size(vectors, 1);

        % Repeat for the number of vectors
        nVectors = size(vectors, 2);
        nSamples = ones(nVectors, 1) * nSamplesScalar;
    end
else
    error('vectors is not the right type!');
end

% Force nSamples to be a column vector
nSamples = force_column_numeric(nSamples);

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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%