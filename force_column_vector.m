function vectors = force_column_vector (vectors, varargin)
%% Transform row vector(s) or array(s) to column vector(s)
% Usage: vectors = force_column_vector (vectors, varargin)
% Explanation:
%       Starting with a cell array of vectors, 
%           this function makes sure each vector is a column vector.
%       If a single vector (may be empty) is provided, 
%           the function makes sure it's a column vector.
%       If a non-vector array is provided, the function 
%           force_column_cell.m is applied 
%           unless 'IgnoreNonVectors' is set to true
%
% Example(s):
%       vector = force_column_vector(vector);
%       vectors = force_column_vector(vectors);
%       force_column_vector({[3, 4], [5; 6], magic(3)})
%
% Outputs:
%       vectors     - vectors transformed
%
% Arguments:
%       vectors     - original vectors
%       varargin    - 'IgnoreNonVectors': whether to ignore non-vectors
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ForceCellOutput': whether to force output as 
%                                           a cell array
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
%                   - 'TreatCharAsScalar': whether to treat character arrays 
%                                           as scalars
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'ToLinearize': whether to linearize a non-vector array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RowInstead': whether to force as row vector instead
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%       cd/isnum.m
%
% Used by:    
%   TODO: Check if some of these can use 
%           match_format_cell or force_column_cell instead
%       cd/annotation_in_plot.m
%       cd/compute_centers_from_edges.m
%       cd/compute_combined_trace.m
%       cd/compute_bins.m
%       cd/compute_default_sweep_info.m
%       cd/compute_peak_decay.m
%       cd/compute_peak_halfwidth.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/create_default_grouping.m
%       cd/create_indices.m
%       cd/extract_columns.m
%       cd/extract_subvectors.m
%       cd/fit_2exp.m
%       cd/force_column_cell.m
%       cd/force_row_vector.m
%       cd/collapse_identical_vectors.m
%       cd/match_format_vectors.m
%       cd/match_format_vector_sets.m
%       cd/match_reciprocals.m
%       cd/match_vector_count.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_neuron_create_simulation_params.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_plot_individual_traces.m
%       cd/plot_bar.m
%       cd/plot_cfit_pulse_response.m
%       cd/plot_histogram.m
%       cd/plot_raster.m
%       cd/plot_window_boundaries.m
%       cd/remove_outliers.m
%       cd/xolotl_set_simparams.m

% File History:
% 2018-10-12 Created by Adam Lu
% 2018-10-27 Added 'IgnoreNonVectors' as an optional argument
% 2018-12-11 Now accepts logical arrays
% 2018-12-28 Now accepts all types of arrays
% 2019-01-03 Now adds ignoreNonVectors to recursive part
% 2019-01-04 Added 'TreatCellAsArray' (default == 'false')
% 2019-01-04 Added 'TreatCellStrAsArray' (default == 'true')
% 2019-01-04 Added 'TreatCharAsScalar' (default == 'true')
% 2019-01-04 Fixed bugs for cellstrs
% 2019-01-08 Added 'ForceCellOutput' as an optional argument
% 2019-01-09 Added 'ToLinearize' as an optional argument
% 2019-01-13 Added 'RowInstead' as an optional argument
% TODO: Deal with 3D arrays
% 

%% Default values for optional arguments
ignoreNonVectorsDefault = false;    % don't ignore non-vectors by default
forceCellOutputDefault = false;     % don't force as cell array by default
treatCellAsArrayDefault = false;% treat cell arrays as many arrays by default
treatCellStrAsArrayDefault = true;  % treat cell arrays of character arrays
                                    %   as an array by default
treatCharAsScalarDefault = true;% treat character arrays as scalars by default
toLinearizeDefault = false;     % whether to linearize a nonvector array
rowInsteadDefault = false;      % whether to force as row vector instead

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
addParameter(iP, 'IgnoreNonVectors', ignoreNonVectorsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ForceCellOutput', forceCellOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellAsArray', treatCellAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellStrAsArray', treatCellStrAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCharAsScalar', treatCharAsScalarDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ToLinearize', toLinearizeDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RowInstead', rowInsteadDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vectors, varargin{:});
ignoreNonVectors = iP.Results.IgnoreNonVectors;
forceCellOutput = iP.Results.ForceCellOutput;
treatCellAsArray = iP.Results.TreatCellAsArray;
treatCellStrAsArray = iP.Results.TreatCellStrAsArray;
treatCharAsScalar = iP.Results.TreatCharAsScalar;
toLinearize = iP.Results.ToLinearize;
rowInstead = iP.Results.RowInstead;

%% Do the job
if iscell(vectors) && ~treatCellAsArray && ...
        ~(iscellstr(vectors) && treatCellStrAsArray)
    % Extract as a cell array
    %   Note: this will have a recursive effect
    vectors = cellfun(@(x) force_column_vector(x, ...
                            'IgnoreNonVectors', ignoreNonVectors, ...
                            'ToLinearize', toLinearize, ...
                            'RowInstead', rowInstead), ...
                    vectors, 'UniformOutput', false);
elseif rowInstead && ~isrow(vectors) || ~rowInstead && ~iscolumn(vectors)
    if isempty(vectors) || ischar(vectors) && treatCharAsScalar
        % Do nothing
    elseif isvector(vectors)
        % Transpose it
        vectors = transpose(vectors);
    else
        % Must be a non-vector
        if ~ignoreNonVectors
            if toLinearize
                % Linearize as a column vector
                vectors = vectors(:);

                % Transpose to a row vector if requested 
                if rowInstead
                    vectors = transpose(vectors);
                end
            else
                % Transform to a column cell array of column vectors
                %   or a row cell array of row vectors
                vectors = force_column_cell(vectors, 'ToLinearize', false, ...
                                            'RowInstead', rowInstead);
            end
        end
    end
else
    % Do nothing
end

%% Force output as a cell array if requested
if forceCellOutput
    % Reassign as a column cell array of column vectors
    %   or a row cell array of row vectors
    vectors = force_column_cell(vectors, 'ToLinearize', false, ...
                                'RowInstead', rowInstead);
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

addRequired(iP, 'vectors', ...
    @(x) isnumeric(x) && isvector(x) || ...
        iscell(x) && all(cellfun(@(x) isnumeric(x) && isvector(x), x)) );

%       vectors     - original vectors
%                   must be a numeric vector or a cell array of numeric vectors
addRequired(iP, 'vectors', ...
    @(x) isnumeric(x) || ...
        iscell(x) && all(cellfun(@(x) isnumeric(x) && isvector(x), x)) );

@(x) assert(isempty(x) || isnumeric(x) || islogical(x) || iscellnumeric(x), ...
            ['vectors must be either empty or a numeric array ', ...
                'or a cell array of numeric arrays!']));

%                   must be a numeric array or a cell array of numeric arrays

addRequired(iP, 'vectors', ...                   % vectors
    @(x) assert(isempty(x) || isnumeric(x) || islogical(x) || ...
                isdatetime(x) || isduration(x) || iscell(x), ...
                ['vectors must be either empty or a numeric array ', ...
                    'or a cell array!']));

if (isnumeric(vectors) || islogical(vectors) || ...
        isdatetime(vectors) || isduration(vectors)) && ~iscolumn(vectors)

%                   specified as a numeric array 
%                       or a cell array of numeric vectors
%                   must be a numeric array or a cell array

% Reassign as a column
vectors = vectors(:);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
