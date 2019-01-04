function vectors = force_column_numeric (vectors, varargin)
%% Transform row numeric vector(s) or numeric array(s) to column numeric vector(s)
% Usage: vectors = force_column_numeric (vectors, varargin)
% Explanation:
%       Starting with a cell array of numeric vectors, 
%           this function makes sure each vector is a column vector.
%       If a single numeric vector (may be empty) is provided, 
%           the function makes sure it's a column vector.
%       If a numeric non-vector array is provided, the function 
%           force_column_cell.m is applied 
%           unless 'IgnoreNonVectors' is set to true
%
% Example(s):
%       vector = force_column_numeric(vector);
%       vectors = force_column_numeric(vectors);
%       force_column_numeric({[3, 4], [5; 6], magic(3)})
%
% Outputs:
%       vectors     - vectors transformed
%
% Arguments:
%       vectors     - original vectors
%       varargin    - 'IgnoreNonVectors': whether to ignore non-vectors
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%       cd/isnum.m
%
% Used by:    
%   TODO: Check if some of these can use 
%           match_format_cell or force_column_cell instead
%       cd/annotation_in_plot.m
%       cd/compute_average_trace.m
%       cd/compute_bins.m
%       cd/compute_default_sweep_info.m
%       cd/compute_peak_decay.m
%       cd/compute_peak_halfwidth.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/create_indices.m
%       cd/extract_subvectors.m
%       cd/fit_2exp.m
%       cd/force_column_cell.m
%       cd/collapse_identical_vectors.m
%       cd/match_format_vectors.m
%       cd/match_format_vector_sets.m
%       cd/match_reciprocals.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_neuron_create_simulation_params.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_plot_individual_traces.m
%       cd/plot_cfit_pulse_response.m
%       cd/plot_raster.m
%       cd/plot_window_boundaries.m
%       cd/xolotl_set_simparams.m

% File History:
% 2018-10-12 Created by Adam Lu
% 2018-10-27 Added 'IgnoreNonVectors' as an optional argument
% 2018-12-11 Now accepts logical arrays
% 2018-12-28 Now accepts all types of arrays
% 2019-01-03 Now adds ignoreNonVectors to recursive part
% TODO: Deal with 3D arrays
% 

%% Default values for optional arguments
ignoreNonVectorsDefault = false;    % don't ignore non-vectors by default

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
addParameter(iP, 'IgnoreNonVectors', ignoreNonVectorsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vectors, varargin{:});
ignoreNonVectors = iP.Results.IgnoreNonVectors;

%% Do the job
if ~iscell(vectors) && ~iscolumn(vectors)
    if isempty(vectors)
        % Do nothing
    elseif isvector(vectors)
        % Reassign as a column
        vectors = vectors(:);
    else
        % Must be a non-vector
        if ~ignoreNonVectors
            % Reassign as a column cell array of column vectors
            vectors = force_column_cell(vectors);
        end
    end
elseif iscell(vectors)
    % Extract as a cell array
    %   Note: this will have a recursive effect
    vectors = cellfun(@(x) force_column_numeric(x, ...
                            'IgnoreNonVectors', ignoreNonVectors), ...
                    vectors, 'UniformOutput', false);
else
    % Do nothing
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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
