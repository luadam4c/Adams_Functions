function outVec = alternate_elements (varargin)
%% Alternate elements between two vectors to create a single vector
% Usage: outVec = alternate_elements (vec1 (opt), vec2 (opt), varargin)
% Explanation:
%       TODO
% Example(s):
%       alternate_elements(1:5, 21:25)
%       alternate_elements('Vectors', {1:5; 21:25; -5:-1})
%       alternate_elements([], [])
%       alternate_elements([], [], 'ReturnNaNIfEmpty', true)
%
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
%                   - 'ReturnNaNIfEmpty': whether to return NaN 
%                                           if all vectors are empty
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/argfun.m
%       cd/count_vectors.m
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
vectorsDefault = [];
returnNaNIfEmptyDefault = false;

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
addParameter(iP, 'Vectors', vectorsDefault, ...
    @iscellnumeric);
addParameter(iP, 'ReturnNaNIfEmpty', returnNaNIfEmptyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
vec1 = iP.Results.vec1;
vec2 = iP.Results.vec2;
vectors = iP.Results.Vectors;
returnNaNIfEmpty = iP.Results.ReturnNaNIfEmpty;

%% Preparation
% Count the vectors
if returnNaNIfEmpty 
    if ~isempty(vectors)
        nVectors = count_vectors(vectors);
    else
        nVectors = 2;
    end
end

% Put vectors together
if isempty(vectors)
    vectors = {vec1, vec2};
end

%% Do the job
% Force as a matrix
vectors = force_matrix(vectors, 'AlignMethod', 'leftAdjustPad');

% Force NaN values for empty vectors if requested
if returnNaNIfEmpty && isempty(vectors)
    vectors = nan(1, nVectors);
end

% Transpose it so that each vector is a row
form1 = transpose(vectors);

% Reshape as a single column vector
outVec = reshape(form1, [], 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%