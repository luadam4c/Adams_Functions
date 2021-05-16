function [indPairs, nPairs] = all_index_pairs (nElements, varargin)
%% List all pairs of indices with or without repetition
% Usage: [indPairs, nPairs] = all_index_pairs (nElements, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       all_index_pairs(3, 'CreateMode', 'all')
%       all_index_pairs(3, 'CreateMode', 'permutation')
%       all_index_pairs(3, 'CreateMode', 'combWithRepeat')
%       all_index_pairs(3, 'CreateMode', 'combination')
%       all_index_pairs(5)
%
% Outputs:
%       indPairs    - indices of pairs, with each row being a pair
%                   specified as a 2-column integer matrix
%       nPairs      - total number of pairs
%                   specified as a positive integer scalar
%
% Arguments:
%       nElements   - total number of elements
%                   must be a positive integer scalar
%       varargin    - 'CreateMode': creation mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'all'           - all possible index pairs
%                       'permutation'   - all permutations
%                       'combWithRepeat'- all combinations with repetition
%                       'combination'   - all combinations without repetition
%                   default == 'combination'
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/crosscorr_profile.m

% File History:
% 2021-05-15 Moved from crosscorr_profile.m
% 

%% Hard-coded parameters
validCreateModes = {'all', 'permutation', 'combination'};

%% Default values for optional arguments
createModeDefault  = 'combination';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'nElements', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'CreateMode', createModeDefault, ...
    @(x) any(validatestring(x, validCreateModes)));

% Read from the Input Parser
parse(iP, nElements, varargin{:});
createMode = validatestring(iP.Results.CreateMode, validCreateModes);

%% Preparation
% Compute the number of total pairs
switch createMode
    case 'all'
        nPairs = nElements ^ 2;
    case 'permutation'
        nPairs = nchoosek(nElements, 2) * factorial(2);
    case 'combination'
        nPairs = nchoosek(nElements, 2);
    otherwise
        error('createMode unrecognized!')
end

%% Do the job
indPairs = zeros(nPairs, 2);
ct = 0;
for i = 1:nElements
    switch createMode
        case 'all'
            withRepetition = true;
            jStart = 1;
        case 'permutation'
            withRepetition = false;
            jStart = 1;
        case 'combWithRepeat'
            withRepetition = true;
            jStart = i;
        case 'combination'
            withRepetition = false;
            jStart = i;
        otherwise
            error('createMode unrecognized!')
    end

    % Only store index pairs from the upper right triangle
    for j = jStart:nElements
        if withRepetition || i ~= j
            % Increment the count
            ct = ct + 1;
            
            % Store the index
            indPairs(ct, 1) = i;
            indPairs(ct, 2) = j;            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%