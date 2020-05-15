function isCellNumericVector = iscellnumericvector (x)
%% Returns whether an input is a cell array of numeric vectors (may be empty)
% Usage: isCellNumericVector = iscellnumericvector (x)
% Explanation:
%       Tests whether the input is a cell array of numeric vectors
%
% Example(s):
%       iscellnumericvector({1:10, 2:20})
%       iscellnumericvector({magic(3), 2:20})
%       iscellnumericvector({'sets', 'lasts'})
%
% Outputs:
%       isCellNumericVector   
%                       - whether the input is a cell array of numeric vectors
%                       specified as a logical scalar
% Arguments:    
%       x               - an input to check
%
% Requires:
%       cd/isnumericvector.m
%
% Used by:
%       cd/compute_combined_trace.m
%       cd/compute_lts_errors.m
%       cd/compute_rms_error.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/count_vectors.m
%       cd/create_average_time_vector.m
%       cd/create_default_grouping.m
%       cd/extract_columns.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/medianfilter.m
%       cd/movingaveragefilter.m

% File History:
% 2018-10-25 Created by Adam Lu
% 2018-10-28 Vectors can now be empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
isCellNumericVector = iscell(x) && all(all(all(cellfun(@isnumericvector, x))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
