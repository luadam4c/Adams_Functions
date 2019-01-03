function vectors = collapse_identical_vectors (vectors, varargin)
%% Collapses identical vectors into a single one
% Usage: vectors = collapse_identical_vectors (vectors, varargin)
% Explanation:
%       TODO
% Example(s):
%       collapse_identical_vectors(1:5)
%       collapse_identical_vectors({(1:5)', (1:5)'})
%       collapse_identical_vectors([(1:5)', (1:5)'])
%       collapse_identical_vectors([(2:6)', (1:5)'])
% Outputs:
%       vectors     - collapsed set of vectors
% Arguments:
%       vectors     - a set of vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/count_vectors.m
%       cd/create_error_for_nargin.m
%       cd/force_column_numeric.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-01-03 Created by Adam Lu
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
addRequired(iP, 'vectors');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, vectors, varargin{:});
% param1 = iP.Results.param1;

%% Preparation

%% Do the job
% If there is just one vector, return
if count_vectors(vectors) == 1
    return
end

% Force as column cell arrays of column vectors
vectorsAll = force_column_numeric(vectors, 'IgnoreNonVectors', false);

% Get the first vector
vectorTemplate = vectorsAll{1};

% If the same as all vectors, return it
if all(cellfun(@(x) isequal(x, vectorTemplate), vectorsAll))
    vectors = vectorTemplate;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% See if the vector set is a matrix
isMatrix = ismatrix(vectors);

% Recombine as a matrix if it was a matrix before
if isMatrix
    vectors = force_matrix(vectors);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%