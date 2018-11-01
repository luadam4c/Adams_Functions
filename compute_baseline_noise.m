function baseNoises = compute_baseline_noise (dataVecs, timeVecs, baseWindows)
%% Computes the baseline noise from a set of data vectors, time vectors and baseline windows
% Usage: baseNoises = compute_baseline_noise (dataVecs, timeVecs, baseWindows)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       baseNoises  - baseline noise value(s)
%                   specified as a numeric vector
% Arguments:
%       dataVecs    - data vector(s)
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric arrays
%       timeVecs    - time vector(s)
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric arrays
%       baseWindows - baseline window(s)
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   must be a numeric array or a cell array of numeric arrays
%
% Requires:
%       cd/compute_rms_error.m
%       cd/find_window_endpoints.m
%       cd/iscellnumeric.m
%
% Used by:    
%       cd/m3ha_plot_individual_traces.m
%       cd/m3ha_run_neuron_once.m

% File History:
% 2018-10-28 Created by Adam Lu
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'dataVecs', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['dataVecs must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'timeVecs', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['timeVecs must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'baseWindows', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['baseWindows must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Read from the Input Parser
parse(iP, dataVecs, timeVecs, baseWindows);

%% Do the job
% Find the indices for the baseline endpoints
endPoints = find_window_endpoints(baseWindows, timeVecs);

% Compute the baseline noise
baseNoises = compute_rms_error(dataVecs, 'Endpoints', endPoints);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%