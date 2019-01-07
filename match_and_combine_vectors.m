function array = match_and_combine_vectors (vecs1, vecs2, varargin)
%% Match vectors and combine into an array
% Usage: array = match_and_combine_vectors (vecs1, vecs2, varargin)
% Explanation:
%       TODO
% Example(s):
%       match_and_combine_vectors(magic(3), 1)
%       match_and_combine_vectors(magic(3), 1:3)
%       match_and_combine_vectors({1:3, 4:6, 7:9}, 1)
% Outputs:
%       array       - combined array(s)
%                   specified as a numeric array 
%                       or a cell array of numeric arrays
% Arguments:
%       vecs1       - first set of vectors
%                   must be a numeric array, a character array or a cell array
%       vecs2       - second set of vectors
%                   must be a numeric array, a character array or a cell array
%       varargin    - TODO
%
% Requires:
%       cd/match_format_vector_sets.m
%
% Used by:
%       cd/create_indices.m

% File History:
% 2019-01-03 Created by Adam Lu
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vecs1');
addRequired(iP, 'vecs2');

% Read from the Input Parser
parse(iP, vecs1, vecs2, varargin{:});

%% Do the job
% Match vector sets and vectors
[vecs1, vecs2] = match_format_vector_sets(vecs1, vecs2, 'MatchVectors', true);

% Combine vectors
if iscell(vecs1) && iscell(vecs2)
    array = cellfun(@(x, y) horzcat(x, y), vecs1, vecs2, ...
                    'UniformOutput', false);
elseif isnumeric(vecs1) && isnumeric(vecs2)
    array = horzcat(vecs1, vecs2);
else
    error('Vectors don''t match!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%