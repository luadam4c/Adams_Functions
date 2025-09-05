function isNumericVector = isnumericvector (x)
%% Returns whether an input is a numeric vector (may be empty)
% Usage: isNumericVector = isnumericvector (x)
% Explanation:
%       Tests whether the input is a numeric vector that may be empty
%
% Example(s):
%       isnumericvector([])
%       isnumericvector(2:20)
%       isnumericvector(magic(3))
%       isnumericvector('sets')
%
% Outputs:
%       isNumericVector - whether the input is a numeric vector (may be empty)
%                       specified as a logical scalar
% Arguments:    
%       x               - an input to check
%
% Requires:
%       cd/isnum.m
%
% Used by:
%       cd/compute_lts_errors.m
%       cd/compute_peak_decay.m
%       cd/compute_peak_halfwidth.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/create_indices.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/find_passive_params.m
%       cd/find_window_endpoints.m
%       cd/iscellnumericvector.m
%       cd/m3ha_neuron_create_sim_params.m
%       cd/match_format_vector_sets.m
%       cd/parse_oscillation.m
%       cd/parse_peaks.m
%       cd/plot_fitted_traces.m
%       cd/plot_histogram.m
%       cd/plot_struct.m
%       cd/plot_traces.m
%       cd/set_axes_properties.m

% File History:
% 2018-10-25 Created by Adam Lu
% 2018-10-28 Vectors can now be empty
% 2019-01-04 Now uses isnum.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
isNumericVector = isnum(x) && (isempty(x) || isvector(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

isNumericVector = isnumeric(x) && (isempty(x) || isvector(x));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
