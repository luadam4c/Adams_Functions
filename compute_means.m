function means = compute_means (vecs, varargin)
%% Computes the mean(s) of vector(s) possibly restricted by endpoint(s)
% Usage: means = compute_means (vecs, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       means       - mean(s) of each vector
%                   specified as a numeric vector 
% Arguments:
%       vecs        - vector(s)
%                   must be a numeric array or a cell array of numeric arrays
%       varargin    - 'Indices': indices for the subvectors to extract 
%                   must be a numeric vector with 2 elements
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == create_indices(endPoints)
%                   - 'Endpoints': endpoints for the subvectors to extract 
%                   must be a numeric vector with 2 elements
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == find_window_endpoints([], vecs)
%                   - 'Windows': value windows to extract 
%                       Note: this assumes that the values are nondecreasing
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric arrays
%                   default == []
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_subvectors.m
%
% Used by:
%       cd/parse_pulse.m
%       cd/parse_pulse_response.m
%
% Related functions:
%       cd/compute_weighted_average.m

% File History:
% 2018-12-17 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
indicesDefault = [];            % set later
endPointsDefault = [];          % set later
windowsDefault = [];            % extract entire trace(s) by default

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
addRequired(iP, 'vecs', ...                  % vectors to extract
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Indices', indicesDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['Indices must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'EndPoints', endPointsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['Windows must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'Windows', windowsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['Windows must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Read from the Input Parser
parse(iP, vecs, varargin{:});

%% Do the job
% Extract subvectors
subVecs = extract_subvectors(vecs, varargin{:});

% Compute the mean(s)
means = cellfun(@mean, subVecs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%