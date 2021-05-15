function [indPairs, nPairs] = create_index_pairs (nElements, varargin)
%% List all pairs of indices without repetition
% Usage: [indPairs, nPairs] = create_index_pairs (nElements, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       create_index_pairs(3)
%       create_index_pairs(5)
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
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
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

%% Default values for optional arguments
param1Default = [];             % default TODO: Description of param1

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
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, nElements, varargin{:});
param1 = iP.Results.param1;

%% Preparation
% Compute the number of total pairs
nPairs = nElements * (nElements - 1) / 2;

%% Do the job
indPairs = zeros(nPairs, 2);
ct = 0;
for i = 1:nElements
    % Only store index pairs from the upper right triangle
    for j = i:nElements
        if i ~= j
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