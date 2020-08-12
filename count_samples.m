function nSamples = count_samples (vectors, varargin)
%% Counts the number of samples whether given an array or a cell array
% Usage: nSamples = count_samples (vectors, varargin)
% Explanation:
%       Uses either numel(x) for vectors, 
%           size(x, 1) for arrays,
%           or cellfun(@numel, x) for cell arrays
%
% Example(s):
%       nSamples = count_samples(data)
%       count_samples(magic(3))
%       count_samples(magic(3), 'TreatMatrixAsVector', true)
%       count_samples(magic(3), 'CountMethod', 'nrows')
%       count_samples(rand(2, 3), 'CountMethod', 'length')
%       count_samples(repmat({repmat({'sdf'}, 3, 1)}, 3, 1))
%       count_samples(repmat({1:3}, 3, 4))
%       count_samples(repmat({1:3}, 3, 4), 'TreatCellNumAsArray', true)
%       count_samples(repmat({'sdf'}, 3, 4))
%       count_samples(repmat({'sdf'}, 3, 4), 'TreatCellStrAsArray', false)
%
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
%                   - 'TreatCellNumAsArray': whether to treat a cell array
%                                       of numeric arrays as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatCellStrAsArray': whether to treat a cell array
%                                       of character arrays as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'CountMethod': method for counting samples
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'veclength' - length of each vector
%                       'length'    - largest dimension of each vector
%                       'numel'     - number of elements
%                       'nrows'     - number of rows
%                       'ncols'     - number of columns
%                       'nRowsEachCol' - number of rows for each column
%                       'nColsEachRow' - number of columns for each row
%                   default == 'nRowsEachCol'
%
% Requires:
%       cd/array_fun.m
%       cd/create_error_for_nargin.m
%       cd/iscellnumeric.m
%       cd/force_column_vector.m
%       cd/match_row_count.m
%
% Used by:    
%       cd/combine_abf_data.m
%       cd/compute_combined_trace.m
%       cd/compute_lts_errors.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sampsizepwr.m
%       cd/compute_sweep_errors.m
%       cd/create_average_time_vector.m
%       cd/create_indices.m
%       cd/create_time_vectors.m
%       cd/extract_columns.m
%       cd/extract_subvectors.m
%       cd/find_window_endpoints.m
%       cd/m3ha_find_decision_point.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_xolotl_plot.m
%       cd/parse_lts.m
%       cd/parse_multiunit.m
%       cd/parse_phase_info.m
%       cd/parse_pulse.m
%       cd/parse_pulse_response.m
%       cd/parse_repetitive_pulses.m
%       cd/parse_spike2_mat.m
%       cd/plot_tuning_curve.m
%       cd/select_similar_values.m
%       cd/test_difference.m

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
% 2019-04-24 Added 'CountMethod' as an optional argument
% 2019-10-03 Added 'TreatCellNumAsArray' as an optional argument
% 

%% Hard-coded parameters
validCountMethods = {'veclength', 'length', 'numel', 'nrows', 'ncols', ...
                        'nRowsEachCol', 'nColsEachRow'};

%% Default values for optional arguments
forceColumnOutputDefault = true;    % force output as a column by default
treatMatrixAsVectorDefault = false; % treat a matrix as many vectors by default
treatRowAsMatrixDefault = false;    % treat a row vector as a vector by default
treatCellAsArrayDefault = false;    % treat cell arrays as many arrays by default
treatCellNumAsArrayDefault = false; % treat cell arrays of numeric arrays
                                    %   as many arrays by default
treatCellStrAsArrayDefault = true;  % treat cell arrays of character arrays
                                    %   as an array by default
countMethodDefault = 'veclength';   % count the length of each vector by default

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
addParameter(iP, 'TreatCellNumAsArray', treatCellNumAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellStrAsArray', treatCellStrAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'CountMethod', countMethodDefault, ...
    @(x) any(validatestring(x, validCountMethods)));

% Read from the Input Parser
parse(iP, vectors, varargin{:});
forceColumnOutput = iP.Results.ForceColumnOutput;
treatMatrixAsVector = iP.Results.TreatMatrixAsVector;
treatRowAsMatrix = iP.Results.TreatRowAsMatrix;
treatCellAsArray = iP.Results.TreatCellAsArray;
treatCellNumAsArray = iP.Results.TreatCellNumAsArray;
treatCellStrAsArray = iP.Results.TreatCellStrAsArray;
countMethod = validatestring(iP.Results.CountMethod, validCountMethods);

%% Do the job
% Decide based on input type
if isempty(vectors)
    nSamples = 0;
elseif iscell(vectors) && ~treatCellAsArray && ...
        ~(iscellnumeric(vectors) && treatCellNumAsArray) && ...
        ~(iscellstr(vectors) && treatCellStrAsArray)
    % Count the number of elements for each array in each cell,
    %   maintaining uniform output if possible
    try
        nSamples = array_fun(@(x) count_samples(x, ...
                            'ForceColumnOutput', forceColumnOutput, ...
                            'TreatCellAsArray', treatCellAsArray, ...
                            'TreatCellNumAsArray', treatCellNumAsArray, ...
                            'TreatCellStrAsArray', treatCellStrAsArray, ...
                            'CountMethod', countMethod), ...
                            vectors, 'UniformOutput', true);
    catch
        nSamples = array_fun(@(x) count_samples(x, ...
                            'ForceColumnOutput', forceColumnOutput, ...
                            'TreatCellAsArray', treatCellAsArray, ...
                            'TreatCellNumAsArray', treatCellNumAsArray, ...
                            'TreatCellStrAsArray', treatCellStrAsArray, ...
                            'CountMethod', countMethod), ...
                            vectors, 'UniformOutput', false);
    end
else
    % Either a non-cell array or a cell array treated as an array
    %   or a cell array of strings treated as an array
    switch countMethod
        case 'veclength'
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
        case 'length'
            % Compute the largest dimension of each vector
            nSamples = length(vectors);
        case 'numel'
            % Count the number of elements
            nSamples = numel(vectors);
        case 'nrows'
            % Count the number of rows
            nSamples = size(vectors, 1);
        case 'ncols'
            % Count the number of columns
            nSamples = size(vectors, 2);
        case 'nRowsEachCol'
            % All the vectors have the same number of samples
            nSamplesScalar = size(vectors, 1);

            % Count the number of vectors
            nVectors = size(vectors, 2);

            % Repeat to make a column vector
            nSamples = match_row_count(nSamplesScalar, nVectors);
        case 'nColsEachRow'
            % All the vectors have the same number of samples
            nSamplesScalar = size(vectors, 2);

            % Count the number of vectors
            nVectors = size(vectors, 1);

            % Repeat to make a column vector
            nSamples = match_row_count(nSamplesScalar, nVectors);
        otherwise
            error('countMethod unrecognized!');
    end
end

% Force nSamples to be a column vector unless requested not to
if forceColumnOutput
    nSamples = force_column_vector(nSamples);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if iscellnumericvector(vectors) || treatMatrixAsVector || ...
        iscellstr(vectors) && ~treatCellStrAsArray
    % Count the number of elements for each vector in each cell
    nSamples = cellfun(@numel, vectors);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
