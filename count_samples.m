function nSamples = count_samples (vectors, varargin)
%% Counts the number of samples whether given an array or a cell array
% Usage: nSamples = count_samples (vectors, varargin)
% Explanation:
%       Uses either numel(x) for vectors, 
%           size(x, 1) for arrays,
%           or cellfun(@numel, x) for cell arrays
% Example(s):
%       nSamples = count_samples(data)
%       count_samples(repmat({repmat({'sdf'}, 3, 1)}, 3, 1))
%       count_samples(repmat({'sdf'}, 3, 4))
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
%                   - 'TreatMatrixAsVector': whether to treat a non-vector array 
%                                           as a single vector
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatRowAsMatrix': whether to treat a row vector
%                                           as many one-element vectors
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatCellAsArray': whether to treat a cell array
%                                           as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatCellStrAsArray': whether to treat a cell array
%                                       of character arrays as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/iscellnumericvector.m
%       cd/force_column_vector.m
%       cd/match_row_count.m
%
% Used by:    
%       cd/compute_combined_trace.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/create_average_time_vector.m
%       cd/create_indices.m
%       cd/extract_columns.m
%       cd/extract_subvectors.m
%       cd/find_window_endpoints.m
%       cd/m3ha_xolotl_plot.m
%       cd/parse_lts.m
%       cd/parse_pulse.m
%       cd/parse_pulse_response.m

% File History:
% 2018-10-10 Created by Adam Lu
% 2018-12-18 Now no longers forces nSamples to be a column vector
% 2019-01-03 Now accepts cell arrays of non-vector arrays
% 2019-01-03 Added 'TreatMatrixAsVector' as an optional argument
% 2019-01-03 Added 'TreatRowAsMatrix' as an optional argument
% 2019-01-04 Now uses isnum.m
% 2019-01-04 Added 'TreatCellAsArray' (default == 'false')
% 2019-01-04 Added 'TreatCellStrAsArray' (default == 'true')
% 2019-01-22 Now uses iscellnumericvector instead of iscellvector
% 2019-01-22 Now returns 0 if vectors is empty
% 2019-01-23 Now maintains uniform output if possible
% 

%% Default values for optional arguments
forceColumnOutputDefault = true;    % force output as a column by default
treatMatrixAsVectorDefault = false; % treat a matrix as many vectors by default
treatRowAsMatrixDefault = false;    % treat a row vector as a vector by default
treatCellAsArrayDefault = false;    % treat cell arrays as many arrays by default
treatCellStrAsArrayDefault = true;  % treat cell arrays of character arrays
                                    %   as an array by default

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
addRequired(iP, 'vectors');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ForceColumnOutput', forceColumnOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatMatrixAsVector', treatMatrixAsVectorDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatRowAsMatrix', treatRowAsMatrixDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellAsArray', treatCellAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellStrAsArray', treatCellStrAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vectors, varargin{:});
forceColumnOutput = iP.Results.ForceColumnOutput;
treatMatrixAsVector = iP.Results.TreatMatrixAsVector;
treatRowAsMatrix = iP.Results.TreatRowAsMatrix;
treatCellAsArray = iP.Results.TreatCellAsArray;
treatCellStrAsArray = iP.Results.TreatCellStrAsArray;

%% Do the job
% Decide based on input type
if isempty(vectors)
    nSamples = 0;
elseif iscell(vectors) && ~treatCellAsArray && ...
        ~(iscellstr(vectors) && treatCellStrAsArray)
    % Count samples in each cell
    if iscellnumericvector(vectors) || treatMatrixAsVector || ...
            iscellstr(vectors) && ~treatCellStrAsArray
        % Count the number of elements for each vector in each cell
        nSamples = cellfun(@numel, vectors);
    else
        % Count the number of elements for each array in each cell,
        %   maintaining uniform output if possible
        try
            nSamples = cellfun(@(x) count_samples(x, ...
                                'ForceColumnOutput', forceColumnOutput, ...
                                'TreatCellAsArray', treatCellAsArray, ...
                                'TreatCellStrAsArray', treatCellStrAsArray), ...
                                vectors, 'UniformOutput', true);
        catch
            nSamples = cellfun(@(x) count_samples(x, ...
                                'ForceColumnOutput', forceColumnOutput, ...
                                'TreatCellAsArray', treatCellAsArray, ...
                                'TreatCellStrAsArray', treatCellStrAsArray), ...
                                vectors, 'UniformOutput', false);
        end
    end
else
    if treatMatrixAsVector || isvector(vectors) && ~treatRowAsMatrix
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
end

% Force nSamples to be a column vector unless requested not to
if forceColumnOutput
    nSamples = force_column_vector(nSamples);
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

elseif isnum(vectors) || iscell(vectors) && treatCellAsArray || ...
        iscellstr(vectors) && treatCellStrAsArray
else
    error('vectors is not the right type!');

@(x) assert(isnum(x) || iscell(x), ...
            'vectors must be either a numeric array or a cell array!'));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
