function vectorsNew = match_vector_count (vectorsOld, nVecsNew, varargin)
%% Expands or truncates a vector set to match a given number of vectors
% Usage: vectorsNew = match_vector_count (vectorsOld, nVecsNew, varargin)
% Explanation:
%       TODO
% Example(s):
%       match_vector_count([1, 2, 3], 6)
%       match_vector_count([1; 2; 3], 6)
%       match_vector_count([1, 2; 3, 4], 7)
%       match_vector_count({[1, 2]}, 3)
%       match_vector_count({[1; 2]}, 3)
% Outputs:
%       vectorsNew    - set of vectors matched
%                   specified as a numeric, cell or struct array
%
% Arguments:
%       vectorsOld  - set of vectors to match
%                   must be a numeric, cell or struct array
%       nVecsNew    - new number of vectors
%                   must be a positive integer scalar
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_vector.m
%
% Used by:
%       cd/match_time_info.m

% File History:
% 2019-02-20 Modified from match_row_counts.m
% 

%% Hard-coded parameters

%% Default values for optional arguments

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
addRequired(iP, 'vectorsOld');
addRequired(iP, 'nVecsNew', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));

% Add parameter-value pairs to the Input Parser

% Read from the Input Parser
parse(iP, vectorsOld, nVecsNew, varargin{:});

%% Preparation
% Force the vectors to be either columns or a cell array of columns
vectorsOld = force_column_vector(vectorsOld, 'IgnoreNonVectors', true);

% Query the old number of vectors
nVecsOld = count_vectors(vectorsOld);

% If vectorsOld is empty 
%   or if the new number of vectors are the same as the old ones, 
%   just return the old vectors
if isempty(vectorsOld) || nVecsNew == nVecsOld
    vectorsNew = vectorsOld;
    return
end

%% Expand or truncate vector set
if nVecsNew > nVecsOld
    % Expand vector set
    vectorsNew = expand_vector_set(vectorsOld, nVecsOld, nVecsNew);
elseif nVecsNew < nVecsOld
    % Truncate vector set
    vectorsNew = truncate_vector_set(vectorsOld, nVecsNew);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vectorsNew = expand_vector_set(vectorsOld, nVecsOld, nVecsNew)
% Expand vector set

% Compute the factor to expand by
factorToExpand = floor(nVecsNew / nVecsOld);

% Compute the remaining number of vectors to fill after expansion
remainingNVecs = mod(nVecsNew, nVecsOld);

% Expand vector set
if iscell(vectorsOld)
    if iscolumn(vectorsOld)
        % First expand by factorToExpand
        vectorsNew = repmat(vectorsOld, [factorToExpand, 1]);

        % Fill in remaining vectors by the first vectors
        vectorsNew = vertcat(vectorsNew, vectorsOld(1:remainingNVecs));
    else
        % First expand by factorToExpand
        vectorsNew = repmat(vectorsOld, [1, factorToExpand]);

        % Fill in remaining vectors by the first vectors
        vectorsNew = horzcat(vectorsNew, vectorsOld(1:remainingNVecs));
    end
else
    % First expand by factorToExpand
    vectorsNew = repmat(vectorsOld, [1, factorToExpand]);

    % Fill in remaining vectors by the first vectors
    vectorsNew = horzcat(vectorsNew, vectorsOld(:, 1:remainingNVecs));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vectorsNew = truncate_vector_set(vectorsOld, nVecsNew)
% Truncate vector set

if iscell(vectorsOld)
    vectorsNew = vectorsOld(1:nVecsNew);
else
    vectorsNew = vectorsOld(:, 1:nVecsNew);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
