function vectors = force_unique_vectors (vectors, varargin)
%% Forces a vector set to have only distinct vectors
% Usage: vectors = force_unique_vectors (vectors, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       vectors     - collapsed set of vectors
% Arguments:
%       vectors     - a set of vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_numeric.m
%       cd/force_matrix.m
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
% See if the vector set is a matrix
isMatrix = ismatrix(vectors);

%% Do the job
% Force as column cell arrays of column vectors
vectors = force_column_numeric(vectors, 'IgnoreNonVectors', false);

% Get unique vectors
vectors = unique(vectors);

% Recombine as a matrix if it was a matrix before
if isMatrix
    vectors = force_matrix(vectors);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%