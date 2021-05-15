function matrix = reorganize_as_matrix (list, indPairs, varargin)
%% Reorganizes a list into a matrix using index pairs
% Usage: matrix = reorganize_as_matrix (list, indPairs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       reorganize_as_matrix([7, 8], [2, 3; 1, 2])
%       reorganize_as_matrix({'cat', 'funny'}, [2, 3; 1, 2])
%
% Outputs:
%       matrix      - matrix returned
%                   specified as a matrix
%
% Arguments:
%       list        - list of things
%                   must be a 1D array
%       indPairs    - index pairs for each item on the list
%                   must be a 2D array with at least 2 columns
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_empty_match.m
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/crosscorr_profile.m

% File History:
% 2021-05-15 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'list');
addRequired(iP, 'indPairs');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, list, indPairs, varargin{:});
param1 = iP.Results.param1;

% Check relationships between arguments
if size(indPairs, 1) < numel(list)
    error('Number of index pairs must be greater than length of list!');
end
if size(indPairs, 2) < 2
    error('Index pairs must have at least 2 columns!');    
end

%% Preparation
% Determine the maximum first index
nRows = max(indPairs(:, 1));

% Determine the maximum second index
nColumns = max(indPairs(:, 2));

%% Do the job
% Initialize matrix
matrix = create_empty_match(list, 'NRows', nRows, 'NColumns', nColumns);

% Fill in the matrix
for iPair = 1:numel(list)
    i = indPairs(iPair, 1);
    j = indPairs(iPair, 2);

    if iscell(matrix)
        matrix{i, j} = list{iPair};
    else
        matrix(i, j) = list(iPair);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%