function [vec1, vec2] = match_reciprocals (vec1, vec2)
%% Check reciprocals, make them column vectors and generate the reciprocal if empty
% Usage: [vec1, vec2] = match_reciprocals (vec1, vec2)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       vec1        - matched vector 1
%                   specified as a column vector
%       vec2        - matched vector 2
%                   specified as a column vector
% Arguments:    
%       vec1        - vector 1
%                   must be a numeric 3D array
%       vec2        - vector 2
%                   must be a numeric 3D array
%
% Requires:
%       cd/argfun.m
%       cd/force_column_vector.m
%
% Used by:    
%       cd/create_time_vectors.m

% File History:
% 2018-10-25 Created by Adam Lu
% 

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
addRequired(iP, 'vec1', ...
    @(x) validateattributes(x, {'numeric'}, {'3d'}));
addRequired(iP, 'vec2', ...
    @(x) validateattributes(x, {'numeric'}, {'3d'}));

% Read from the Input Parser
parse(iP, vec1, vec2);

%% Preparation
% Force inputs to be column numeric vectors
[vec1, vec2] = argfun(@force_column_vector, vec1, vec2);

%% Do the job
% Match or check based on which one is empty
if isempty(vec1) && isempty(vec2)
    vec1 = 1;
    vec2 = 1;
elseif ~isempty(vec1) && isempty(vec2)
    vec2 = 1 ./ vec1;
elseif isempty(vec1) && ~isempty(vec2)
    vec1 = 1 ./ vec2;
else
    if ~isequal(vec1, 1 ./ vec2)
        vec1 = [];
        vec2 = [];
        fprintf('Reciprocals do not match! Both made empty!\n\n');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

vec1 = vec1(:);
vec2 = vec2(:);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
