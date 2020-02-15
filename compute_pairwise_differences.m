function diffVectors = compute_pairwise_differences (vectors, varargin)
%% Computes pairwise differences of vectors
% Usage: diffVectors = compute_pairwise_differences (vectors, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       compute_pairwise_differences(magic(2))
%       compute_pairwise_differences(magic(3))
%       compute_pairwise_differences(magic(4))
%       compute_pairwise_differences(magic(5))
%
% Outputs:
%       diffVectors - difference vectors
%                   specified as an array
%
% Arguments:
%       vectors     - vectors
%                   must be a an array
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_matrix.m
%       cd/vecfun.m
%
% Used by:
%       cd/test_difference.m

% File History:
% 2020-02-15 Created by Adam Lu
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
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, vectors, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Force as a matrix, padding NaNs if needed
vectors = force_matrix(vectors, 'AlignMethod', 'leftAdjustPad');

% Count the number of vectors
nVectors = size(vectors, 2);

% List all pairs of indices
%   Note: nchoosek returns a nPairs x 2 matrix, 
%           so transpose to make a 2 x nPairs matrix
allIndices = transpose(nchoosek(1:nVectors, 2));

% Compute all differences
diffVectors = vecfun(@(x) vectors(:, x(2)) - vectors(:, x(1)), allIndices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%