function outVecs = nan_except (inVecs, indices, varargin)
%% Make vectors all NaNs except for selected indices
% Usage: outVecs = nan_except (inVecs, indices, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       nan_except(magic(3), [1, 2])
%       nan_except(magic(3), [1, 3, 2; 2, 1, 3])
%
% Outputs:
%       outVecs     - TODO: Description of outVecs
%                   specified as a TODO
%
% Arguments:
%       inVecs      - TODO: Description of inVecs
%                   must be a TODO
%       indices     - TODO: Description of indices
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/create_empty_match.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
%       cd/match_format_vector_sets.m
%
% Used by:
%       cd/extract_subvectors.m

% File History:
% 2020-05-12 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'inVecs');
addRequired(iP, 'indices');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, inVecs, indices, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Create an empty match of inVecs
outVecs = create_empty_match(inVecs);

% Force as column cell arrays of column vectors
[inVecs, outVecs] = ...
    argfun(@(x) force_column_vector(x, 'ForceCellOutput', true), ...
            inVecs, outVecs);

% Match indices to inVecs
indices = match_format_vector_sets(indices, inVecs);

% Replace seleted indices
outVecs = cellfun(@(a, b, c) replace_selected_indices(a, b, c), ...
                    outVecs, inVecs, indices, 'UniformOutput', false);

% Force as a matrix
outVecs = force_matrix(outVecs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outVecs = replace_selected_indices (outVecs, inVecs, indices)

% Replace with actual values for indices
outVecs(indices, :) = inVecs(indices, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%