function isNumericVector = isnumericvector (x)
%% Returns whether an input is a numeric vector (may be empty)
% Usage: isNumericVector = isnumericvector (x)
% Explanation:
%       Tests whether the input is a cell array of numeric vectors
% Example(s):
%       isnumericvector([])
%       isnumericvector(2:20)
%       isnumericvector(magic(3))
%       isnumericvector('sets')
% Outputs:
%       isNumericVector - whether the input is a numeric vector (may be empty)
%                       specified as a logical scalar
% Arguments:    
%       x               - an input to check
%
% Used by:
%       cd/compute_rms_error.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/extract_subvectors.m
%       cd/find_passive_params.m
%       cd/iscellnumericvector.m
%       cd/m3ha_create_simulation_params.m
%       cd/m3ha_plot_individual_traces.m
%       cd/match_format_vectors.m

% File History:
% 2018-10-25 Created by Adam Lu
% 2018-10-28 Vectors can now be empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

%% Do the job
isNumericVector = isnumeric(x) && (isempty(x) || isvector(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%