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
%       force_column_vector({1:3, 4:7})
%       force_column_vector({1:3, 4:7}, 'TreatCellNumAsArray', true)
%       force_column_vector({1:3, 4:7}, 'CombineAcrossCells', true)
%       force_column_vector({ones(2, 1), magic(3)}, 'ToLinearize', true)
%       force_column_vector({ones(2, 1), magic(3)}, 'CombineAcrossCells', true)
%       force_column_vector({ones(2, 1), magic(3)}, 'ToLinearize', true, 'CombineAcrossCells', true)
%       force_column_vector({{ones(1, 2), ones(1, 2)}, {[], []}}, 'CombineAcrossCells', true)
%       force_column_vector({{ones(2, 1), ones(2, 1)}, {[], []}}, 'ToLinearize', true, 'CombineAcrossCells', true)
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
%                   - 'TreatCellNumAsArray': whether to treat a cell array
%                                       of numeric arrays as a single array
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
%                   - 'CombineAcrossCells': whether to combine across cells
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RowInstead': whether to force as row vector instead
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_cell.m
%       cd/force_matrix.m
%       cd/iscellnumeric.m
%       cd/isnum.m
%
% Used by:
%       cd/compute_running_windows.m
%       cd/compute_sampling_interval.m
%       cd/create_shifted_vectors.m
%       cd/plot_correlation_coefficient.m
%   TODO: Check if some of these can use 
%           match_format_cell or force_column_cell instead
%       cd/adjust_edges.m
%       cd/align_vectors_by_index.m
%       cd/annotation_in_plot.m
%       cd/combine_variables_across_tables.m
%       cd/compute_activation_profile.m
%       cd/compute_centers_from_edges.m
%       cd/compute_combined_trace.m
%       cd/create_default_grouping.m
%       cd/compute_derivative_trace.m
%       cd/compute_bins.m
%       cd/compute_default_sweep_info.m
%       cd/compute_gabab_conductance.m
%       cd/compute_grouped_histcounts.m
%       cd/compute_lts_errors.m
%       cd/compute_peak_decay.m
%       cd/compute_peak_halfwidth.m
%       cd/compute_relative_event_times.m
%       cd/compute_rms_error.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/compute_weighted_average.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/create_indices.m
%       cd/create_looped_params.m
%       cd/detect_spikes_current_clamp.m
%       cd/extract_columns.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/fit_2exp.m
%       cd/find_closest.m
%       cd/find_window_endpoints.m
%       cd/force_column_cell.m
%       cd/force_data_as_matrix.m
%       cd/force_matrix.m
%       cd/force_row_vector.m
%       cd/collapse_identical_vectors.m
%       cd/m3ha_compute_statistics.m
%       cd/match_format_vectors.m
%       cd/match_format_vector_sets.m
%       cd/match_reciprocals.m
%       cd/match_vector_count.m
%       cd/m3ha_network_change_params.m
%       cd/m3ha_network_launch.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_neuron_create_sim_params.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_simulate_population.m
%       cd/movingaveragefilter.m
%       cd/nan_except.m
%       cd/parse_peaks.m
%       cd/parse_pleth_trace.m
%       cd/parse_repetitive_pulses.m
%       cd/plot_fitted_traces.m
%       cd/plot_bar.m
%       cd/plot_cfit_pulse_response.m
%       cd/plot_error_bar.m
%       cd/plot_histogram.m
%       cd/plot_table_parallel.m
%       cd/plot_raster.m
%       cd/plot_traces.m
%       cd/plot_window_boundaries.m
%       cd/remove_outliers.m
%       cd/select_similar_values.m
%       cd/unique_groups.m
%       cd/vecfun.m
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
% 2019-04-24 Added 'CombineAcrossCells' as an optional argument
% 2019-04-27 Fixed the case when 'ToLinearize' and 'CombineAcrossCells' are
%               both true
% 2019-10-03 Added 'TreatCellNumAsArray' as an optional argument
% 
% TODO: Need to debug 'CombineAcrossCells' and 'ToLinearize'
% TODO: Deal with 3D arrays
% 

%% Default values for optional arguments
ignoreNonVectorsDefault = false;    % don't ignore non-vectors by default
forceCellOutputDefault = false;     % don't force as cell array by default
treatCellAsArrayDefault = false;% treat cell arrays as many arrays by default
treatCellNumAsArrayDefault = false; % treat cell arrays of numeric arrays
                                    %   as many arrays by default
treatCellStrAsArrayDefault = true;  % treat cell arrays of character arrays
                                    %   as an array by default
treatCharAsScalarDefault = true;% treat character arrays as scalars by default
toLinearizeDefault = false;     % whether to linearize a nonvector array
combineAcrossCellsDefault = false;  % whether to combine across cells
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
addParameter(iP, 'TreatCellNumAsArray', treatCellNumAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellStrAsArray', treatCellStrAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCharAsScalar', treatCharAsScalarDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ToLinearize', toLinearizeDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'CombineAcrossCells', combineAcrossCellsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RowInstead', rowInsteadDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vectors, varargin{:});
ignoreNonVectors = iP.Results.IgnoreNonVectors;
forceCellOutput = iP.Results.ForceCellOutput;
treatCellAsArray = iP.Results.TreatCellAsArray;
treatCellNumAsArray = iP.Results.TreatCellNumAsArray;
treatCellStrAsArray = iP.Results.TreatCellStrAsArray;
treatCharAsScalar = iP.Results.TreatCharAsScalar;
toLinearize = iP.Results.ToLinearize;
combineAcrossCells = iP.Results.CombineAcrossCells;
rowInstead = iP.Results.RowInstead;

%% Do the job
if iscell(vectors) && ~treatCellAsArray && ...
        ~(iscellnumeric(vectors) && treatCellNumAsArray) && ...
        ~(iscellstr(vectors) && treatCellStrAsArray)
    if combineAcrossCells
        % Break up contents in cells
        while iscell(vectors) && (toLinearize || ~iscellnumeric(vectors))
            % Linearize the cell array
            vectors = force_column_cell(vectors, 'ToLinearize', true);

            % Force as a matrix (horizontally-concatenated vectors)
            if toLinearize
                vectors = force_matrix(vectors, 'Verbose', false, ...
                                        'AlignMethod', 'leftAdjustPad');
            else
                vectors = force_matrix(vectors, 'AlignMethod', 'none');
            end
            
            % Apply the function recursively on the horizontally-concatenated
            %   vectors
            vectors = force_column_vector(vectors, ...
                                    'IgnoreNonVectors', ignoreNonVectors, ...
                                    'ToLinearize', toLinearize, ...
                                    'RowInstead', rowInstead);
        end

        % If linearized, remove padded NaNs when forced as a matrix
        if toLinearize
            vectors(isnan(vectors)) = [];
        end
    else
        % Apply the function recursively on each cell
        vectors = cellfun(@(x) force_column_vector(x, ...
                                'IgnoreNonVectors', ignoreNonVectors, ...
                                'ToLinearize', toLinearize, ...
                                'RowInstead', rowInstead), ...
                        vectors, 'UniformOutput', false);
    end
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
