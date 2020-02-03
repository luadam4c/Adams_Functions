function varargout = argfun (myFunction, varargin)
%% Applies a function to each input argument
% Usage: varargout = argfun (myFunction, varargin)
% Explanation:
%       This function applies the first argument (a function)
%           to each of the rest of the arguments
%       The first argument becomes the first output,
%           the second argument becomes the second output, and so on ...
%
% Example(s):
%       [a, b] = argfun(@sum, 1:10, magic(3))
%
% Outputs:
%       varargout   - outputs of the function applied to each input argument
%
% Arguments:    
%       myFunction  - a custom function
%                   must be a function handle
%       varargin    - input arguments
%
% Requires:
%       cd/array_fun.m
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/adjust_edges.m
%       cd/align_subplots.m
%       cd/alternate_elements.m
%       cd/annotation_in_plot.m
%       cd/combine_abf_data.m
%       cd/compute_average_pulse_response.m
%       cd/compute_bins.m
%       cd/compute_default_sweep_info.m
%       cd/compute_grouped_histcounts.m
%       cd/compute_peak_decay.m
%       cd/compute_peak_halfwidth.m
%       cd/compute_relative_event_times.m
%       cd/compute_rms_error.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sampsizepwr.m
%       cd/compute_sweep_errors.m
%       cd/compute_lts_errors.m
%       cd/construct_fullpath.m
%       cd/convert_units.m
%       cd/copy_into.m
%       cd/create_average_time_vector.m
%       cd/create_labels_from_numbers.m
%       cd/create_looped_params.m
%       cd/create_pleth_EEG_movies.m
%       cd/create_indices.m
%       cd/extract_subvectors.m
%       cd/filter_and_extract_pulse_response.m
%       cd/find_closest.m
%       cd/find_passive_params.m
%       cd/find_window_endpoints.m
%       cd/m3ha_compute_statistics.m
%       cd/m3ha_extract_component_errors.m
%       cd/m3ha_network_analyze_spikes.m
%       cd/m3ha_network_plot_essential.m
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_neuron_create_simulation_params.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_plot_bar3.m
%       cd/m3ha_plot_figure03.m
%       cd/m3ha_plot_figure04.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/m3ha_network_change_params.m
%       cd/m3ha_network_launch.m
%       cd/m3ha_network_raster_plot.m
%       cd/m3ha_network_single_neuron.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_rank_neurons.m
%       cd/m3ha_simulate_population.m
%       cd/m3ha_xolotl_plot.m
%       cd/match_format_vector_sets.m
%       cd/match_reciprocals.m
%       cd/parse_all_abfs.m
%       cd/parse_assyst_swd.m
%       cd/parse_atf_swd.m
%       cd/parse_current_family.m
%       cd/parse_iox.m
%       cd/parse_ipsc.m
%       cd/parse_lts.m
%       cd/parse_multiunit.m
%       cd/parse_phase_info.m
%       cd/parse_pulse_response.m
%       cd/parse_spike2_mat.m
%       cd/plot_bar.m
%       cd/plot_chevron.m
%       cd/plot_cfit_pulse_response.m
%       cd/plot_error_bar.m
%       cd/plot_fitted_traces.m
%       cd/plot_histogram.m
%       cd/plot_measures.m
%       cd/plot_raw_multiunit.m
%       cd/plot_repetitive_protocols.m
%       cd/plot_relative_events.m
%       cd/plot_traces.m
%       cd/plot_traces_abf.m
%       cd/test_var_difference.m
%       cd/transform_vectors.m
%       cd/xolotl_add_current_injection.m

% File History:
% 2018-10-25 Created by Adam Lu
% 2018-12-17 Now uses create_labels_from_numbers.m
% 2018-12-17 Now uses cellfun instead of structfun
% 

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
addRequired(iP, 'myFunction', ...                  % a custom function
    @(x) validateattributes(x, {'function_handle'}, {'scalar'}));

% Read from the Input Parser
parse(iP, myFunction);

% Count the number of inputs (number of arguments excluding the function handle)
nInputs = nargin - 1;

% If there are more outputs requested than inputs, return error
if nInputs < nargout
    disp('Cannot request more ''outputs'' than provided ''inputs''!!');
    error(create_error_for_nargin(mfilename));
end

%% Do the job
% Extract all fields from the structure
varargout = array_fun(myFunction, varargin, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

myFieldNames = arrayfun(@(x) ['Arg', num2str(x)], 1:nInputs, ...
                        'UniformOutput', false);
myFieldNames = create_labels_from_numbers(1:nInputs, 'Prefix', 'Arg');
myFieldNames = create_labels_from_numbers(1:nInputs);

%       cd/create_labels_from_numbers.m

% Generate field names for the input arguments
%   Note: field names cannot start with a number, so a prefix is necessary
myFieldNames = create_labels_from_numbers(1:nInputs, 'Prefix', 'a');

% Place all arguments in an input structure
%   Note: varargin is a row cell array, so operate along dimension 2
%           to make each element a field
myStructInputs = cell2struct(varargin, myFieldNames, 2);

% Pass the function to all fields in the input structure
myStructOutputs = structfun(myFunction, myStructInputs, ...
                            'UniformOutput', false);

% Extract all fields from the structure
varargout = struct2cell(myStructOutputs);

eval(sprintf('help %s', mfilename));

print_help(mfilename);
return

function check_nargin(nArgInTarget, nArgInActual, functionName)
    
if nArgInActual < nArgInTarget
    error(create_error_for_nargin(functionName));
end

if nargout > nInputs
    error(['Cannot request more outputs than provided inputs, ', ...
            'type ''help %s'' for usage'], mfilename);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
