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
%       cd/iscellvector.m
%       cd/force_column_vector.m
%
% Used by:
%       cd/collapse_identical_vectors.m
%       cd/compute_all_pulse_responses.m
%       cd/compute_combined_trace.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/extract_channel.m
%       cd/extract_elements.m
%       cd/identify_repetitive_pulses.m
%       cd/m3ha_average_by_group.m
%       cd/m3ha_plot_individual_traces.m
%       cd/parse_pulse_response.m
%       cd/plot_protocols.m
%       cd/plot_traces.m

% File History:
% 2018-10-10 Created by Adam Lu
% 2019-01-03 Now returns a vector if input is a cell array of non-vectors
% 2019-01-03 Added 'TreatMatrixAsVector' as an optional argument
% 2019-01-03 Added 'TreatRowAsMatrix' as an optional argument
% 2018-01-03 Added 'ForceColumnOutput' as an optional argument with default true
% 2019-01-04 Added 'TreatCellAsArray' (default == 'false')
% 2019-01-04 Added 'TreatCellStrAsArray' (default == 'true')
% 

%% Default values for optional arguments
forceColumnOutputDefault = true;    % force output as a column by default
treatMatrixAsVectorDefault = false; % treat a matrix as many vectors by default
treatRowAsMatrixDefault = false;    % treat a row vector as a vector by default
treatCellAsArrayDefault = false;% treat cell arrays as many arrays by default
treatCellStrAsArrayDefault = true;  % treat cell arrays of character arrays
                                    %   as an array by default

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
if iscell(vectors) && ~treatCellAsArray && ...
        ~(iscellstr(vectors) && treatCellStrAsArray)
    % Count number of vectors
    if iscellvector(vectors) || treatMatrixAsVector || ...
            iscellstr(vectors) && ~treatCellStrAsArray
        % Count number of cells
        nVectors = numel(vectors);
    else
        % TODO: the case when not uniform output
        nVectors = cellfun(@count_vectors, vectors);
    end
else
    if treatMatrixAsVector || isvector(vectors) && ~treatRowAsMatrix
        nVectors = 1;
    else
        nVectors = size(vectors, 2);
    end
end

% Force nVectors to be a column vector unless requested not to
if forceColumnOutput
    nVectors = force_column_vector(nVectors);
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

if isnum(vectors) || iscell(vectors) && treatCellAsArray

%                   must be a numeric array or a cell array
addRequired(iP, 'vectors', ...                   % vectors
    @(x) assert(isnum(x) || iscell(x), ...
                'vectors must be either a numeric array or a cell array!'));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
