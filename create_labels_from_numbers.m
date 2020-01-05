function labels = create_labels_from_numbers (numbers, varargin)
%% Creates a cell array of labels from an array of numbers with an optional prefix or suffix
% Usage: labels = create_labels_from_numbers (numbers, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       labels = create_labels_from_numbers([3, 7, 5])
%       labels = create_labels_from_numbers(1:3, 'Suffix', ' Mississippi')
%       labels = create_labels_from_numbers(1:3, 'Prefix', 'Husband #')
%       labels = create_labels_from_numbers(1:3, 'Prefix', "Katie ")
%       labels = create_labels_from_numbers(1:3, 'Prefix', 'Katie', 'Delimiter', '_')
%       labels = create_labels_from_numbers(1:3, 'Prefix', 'Make', 'Suffix', 'Wish', 'Delimiter', ' ')
%       labels = create_labels_from_numbers(1:3, 'Prefix', ["Katie ", "Mark ", "Matt "])
%       labels = create_labels_from_numbers(1:3, 'Suffix', {' one', ' two', ' three'})
%
% Outputs:
%       labels     - labels created
%                   specified as a cell array of character vectors
%
% Arguments:
%       numbers     - numbers
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%       varargin    - 'ForceColumnOutput': whether to force output as a column
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Prefix': string to place before each number
%                   must be a string scalar or a character vector
%                       or a cell array of character vectors
%                   default == ''
%                   - 'Suffix': string to place after each number
%                   must be a string scalar or a character vector
%                       or a cell array of character vectors
%                   default == ''
%                   - 'Delimiter': delimiter used
%                   must be a string scalar or a character vector
%                   default == ''
%                   - Any other parameter-value pair for the convert_to_char() function
%
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/convert_to_char.m
%       cd/force_string_end.m
%       cd/force_string_start.m
%       cd/match_dimensions.m
%
% Used by:
%       cd/atfwrite.m
%       cd/combine_variables_across_tables.m
%       cd/compute_all_pulse_responses.m
%       cd/create_default_grouping.m
%       cd/create_row_labels.m
%       cd/create_simulation_output_filenames.m
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/decide_on_filebases.m
%       cd/force_data_as_matrix.m
%       cd/m3ha_compute_gabab_ipsc.m
%       cd/m3ha_extract_component_errors.m
%       cd/m3ha_network_launch.m
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_plot_figure04.m
%       cd/m3ha_plot_figure05.m
%       cd/m3ha_simulate_population.m
%       cd/m3ha_xolotl_plot.m
%       cd/parse_current_family.m
%       cd/plot_chevron.m
%       cd/plot_grouped_histogram.m
%       cd/plot_histogram.m
%       cd/plot_table_parallel.m
%       cd/plot_raster.m
%       cd/plot_spike_density_multiunit.m
%       cd/plot_struct.m
%       cd/plot_swd_raster.m
%       cd/plot_traces.m
%       cd/plot_tuning_curve.m
%       cd/renamevars.m
%       cd/save_all_zooms.m
%       cd/transpose_table.m

% File History:
% 2018-12-17 Created by Adam Lu
% 2019-09-06 Added 'Delimiter' as an optional argument
% 2019-12-22 Now allows prefix and suffix to be cell arrays as well

%% Hard-coded parameters

%% Default values for optional arguments
forceColumnOutputDefault = true;    % force output as a column by default
prefixDefault = '';     % no string to place before each number by default
suffixDefault = '';     % no string to place after each number by default
delimiterDefault = '';  % delimiter is empty string by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'numbers', ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'3d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ForceColumnOutput', forceColumnOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['Prefix must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'Suffix', suffixDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['Suffix must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, numbers, varargin{:});
forceColumnOutput = iP.Results.ForceColumnOutput;
prefix = iP.Results.Prefix;
suffix = iP.Results.Suffix;
delimiter = iP.Results.Delimiter;

% Keep unmatched arguments for the convert_to_char() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Force the numbers as a column vector if requested
if forceColumnOutput
    numbers = numbers(:);
end

% Modify prefix or suffix if necessary
prefix = force_string_end(prefix, delimiter, 'OnlyIfNonempty', true);
suffix = force_string_start(suffix, delimiter, 'OnlyIfNonempty', true);

% Force as a cell array
if ischar(prefix)
    prefix = {prefix};
end
if ischar(suffix)
    suffix = {suffix};
end

% Match the dimensions with numbers
[prefix, suffix] = ...
    argfun(@(x) match_dimensions(x, size(numbers)), prefix, suffix);

%% Do the job
% Create the labels
%   Note: don't use strcat() as it omits spaces in some cases
if iscell(prefix) && iscell(suffix)
    labels = ...
        cellfun(@(x, y, z) [x, convert_to_char(y, otherArguments{:}), z], ...
                    prefix, num2cell(numbers), suffix, 'UniformOutput', false);
else
    labels = ...
        arrayfun(@(x, y, z) [x, convert_to_char(y, otherArguments{:}), z], ...
                    prefix, numbers, suffix, 'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
