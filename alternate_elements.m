function outVec = alternate_elements (varargin)
%% Alternate elements between two vectors to create a single vector
% Usage: outVec = alternate_elements (vec1 (opt), vec2 (opt), varargin)
% Explanation:
%       TODO
% Example(s):
%       alternate_elements(1:5, 21:25)
%       alternate_elements('Vectors', {1:5; 21:25; -5:-1})
% Outputs:
%       outVec      - output vector
%                   specified as a numeric vector
% Arguments:
%       vec1        - (opt) input vector 1
%                   must be a numeric vector
%       vec2        - (opt) input vector 2
%                   must be a numeric vector
%       varargin    - 'Vectors': multiple vectors as a cell array
%                   must be a cell array of numeric vectors
%                   default == []
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/argfun.m
%       cd/force_matrix.m
%       cd/iscellnumeric.m
%
% Used by:
%       cd/parse_multiunit.m

% File History:
% 2019-08-13 Adapted from parse_multiunit.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
vec1Default = [];
vec2Default = [];
VectorsDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add optional inputs to the Input Parser
addOptional(iP, 'vec1', vec1Default, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addOptional(iP, 'vec2', vec2Default, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Vectors', VectorsDefault, ...
    @iscellnumeric);

% Read from the Input Parser
parse(iP, varargin{:});
vec1 = iP.Results.vec1;
vec2 = iP.Results.vec2;
vectors = iP.Results.Vectors;

%% Preparation
% Put vectors together
if isempty(vectors)
    vectors = {vec1, vec2};
end

%% Do the job
% Force as a matrix
vectors = force_matrix(vectors, 'AlignMethod', 'leftAdjustPad');

% Transpose it so that each vector is a row
form1 = transpose(vectors);

% Reshape as a single column vector
outVec = reshape(form1, [], 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%