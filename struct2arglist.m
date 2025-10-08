function argList = struct2arglist (structure, varargin)
%% Converts a scalar structure to an argument list
% Usage: argList = struct2arglist (structure, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       argList     - argument list (parameter value pairs)
%                   specified as a cell array of character vectors
%
% Arguments:
%       structure   - structure with field names as parameters
%                   must be a scalar structure
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_cell.m
%
% Used by:
%       cd/array_fun.m
%       cd/annotation_in_plot.m
%       cd/apply_to_all_cells.m
%       cd/apply_over_cells.m
%       cd/combine_phase_numbers.m
%       cd/compute_bins.m
%       cd/convert_to_rank.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_plot_bar3.m
%       cd/medianfilter.m
%       cd/parse_peaks.m
%       cd/plot_autocorrelogram.m
%       cd/plot_ball_stick.m
%       cd/plot_bar.m
%       cd/plot_calcium_imaging_traces.m
%       cd/plot_correlation_coefficient.m
%       cd/plot_grouped_histogram.m
%       cd/plot_grouped_jitter.m
%       cd/plot_grouped_scatter.m
%       cd/plot_histogram.m
%       cd/plot_traces.m
%       cd/plot_tuning_curve.m
%       cd/plot_violin.m
%       cd/plot_spike_density_multiunit.m
%       cd/solve_function_at_value.m
%       cd/test_normality.m

% File History:
% 2018-12-28 Moved from annotation_in_plot.m
% 2020-03-09 Now returns empty if empty
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'structure', ...
    @(x) validateattributes(x, {'struct'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, structure, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Return empty if empty
if isempty(structure)
    argList = {};
    return
end

% Get all the parameter names
names = fieldnames(structure);

% Get all the parameter values
values = struct2cell(structure);

% Force as a column cell array
argList = force_column_cell(transpose([names, values]), 'ToLinearize', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%