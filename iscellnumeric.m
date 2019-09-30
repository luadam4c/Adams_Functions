function isCellNumeric = iscellnumeric (x)
%% Returns whether an input is a cell array of numeric arrays
% Usage: isCellNumeric = iscellnumeric (x)
% Explanation:
%       Tests whether the input is a cell array of numeric arrays
%
% Example(s):
%       iscellnumeric({1:10, 2:20})
%       iscellnumeric({'sets', 'lasts'})
%
% Outputs:
%       isCellNumeric   - whether the input is a cell array of numeric arrays
%                       specified as a logical scalar
% Arguments:    
%       x               - an input to check
%
% Requires:
%       cd/isnum.m
%
% Used by:
%       cd/alternate_elements.m
%       cd/compute_baseline_noise.m
%       cd/compute_relative_event_times.m
%       cd/compute_residuals.m
%       cd/compute_rms_error.m
%       cd/compute_sampsizepwr.m
%       cd/create_default_endpoints.m
%       cd/find_pulse_endpoints.m
%       cd/find_pulse_response_endpoints.m
%       cd/find_window_endpoints.m
%       cd/force_column_vector.m
%       cd/force_row_vector.m
%       cd/extract_columns.m
%       cd/m3ha_plot_individual_traces.m
%       cd/match_array_counts.m
%       cd/match_format_vector_sets.m
%       cd/parse_multiunit.m
%       cd/parse_pulse.m
%       cd/parse_pulse_response.m
%       cd/plot_raster.m
%       cd/unique_groups.m

% File History:
% 2018-10-24 Created by Adam Lu
% 2018-12-18 Now accepts 2D arrays as input
% 2019-01-04 Now uses isnum.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

%% Do the job
isCellNumeric = iscell(x) && all(all(all(cellfun(@isnum, x))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

isCellNumeric = iscell(x) && all(all(cellfun(@isnumeric, x)));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
