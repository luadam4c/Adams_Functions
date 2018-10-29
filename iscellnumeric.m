function isCellNumeric = iscellnumeric (x)
%% Returns whether an input is a cell array of numeric arrays
% Usage: isCellNumeric = iscellnumeric (x)
% Explanation:
%       Tests whether the input is a cell array of numeric arrays
% Example(s):
%       iscellnumeric({1:10, 2:20})
%       iscellnumeric({'sets', 'lasts'})
% Outputs:
%       isCellNumeric   - whether the input is a cell array of numeric arrays
%                       specified as a logical scalar
% Arguments:    
%       x               - an input to check
%
% Used by:
%       cd/compute_average_trace.m
%       cd/compute_baseline_noise.m
%       cd/compute_residuals.m
%       cd/compute_rms_error.m
%       cd/find_pulse_endpoints.m
%       cd/find_pulse_response_endpoints.m
%       cd/find_window_endpoints.m
%       cd/force_column_numeric.m
%       cd/force_row_numeric.m
%       cd/extract_columns.m
%       cd/match_array_counts.m
%       cd/match_format_vectors.m
%       cd/parse_pulse.m
%       cd/parse_pulse_response.m
%       cd/plot_raster.m

% File History:
% 2018-10-24 Created by Adam Lu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

%% Do the job
isCellNumeric = iscell(x) && all(cellfun(@isnumeric, x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%