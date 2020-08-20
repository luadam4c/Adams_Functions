function vectorsInOrder = extract_in_order (vectorSets, varargin)
%% Reorganizes a set of vectors by the first vectors, the second vectors, ... and so on
% Usage: vectorsInOrder = extract_in_order (vectorSets, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       vectorsInOrder     - TODO: Description of vectorsInOrder
%                   specified as a TODO
%
% Arguments:
%       vectorSets     - TODO: Description of vectorSets
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:

% File History:
% 2019-10-03 Moved from plot_relative_events.m
%% TODO: Merge with extract_columns.m?
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1    % TODO: 1 might need to be changed
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'vectorSets');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, vectorSets, varargin{:});
param1 = iP.Results.param1;

% Keep unmatched arguments for the TODO() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% TODO

%% Do the job
% Force as a column cell array
vectorSets = force_column_cell(vectorSets);

% Count the number of vectors in each set
nVectorsEachSet = count_vectors(vectorSets);

% Compute the maximum
maxNVectorsEachSet = max(nVectorsEachSet);

% Place all first vectors in the first cell, second vectors in the second cell,
%   and so on
vectorsInOrder = ...
    arrayfun(@(x) cellfun(@(y) extract_element_by_index(y, x, true), ...
                            vectorSets, 'UniformOutput', false), ...
            transpose(1:maxNVectorsEachSet), 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function element = extract_element_by_index (vec, index, returnNaNInstead)
%% Extract an element or return empty if the element doesn't exist
% TODO: Merge with extract_elements

if numel(vec) >= index
    if iscell(vec)
        element = vec{index};
    else
        element = vec(index);
    end
else
    if returnNaNInstead
        element = NaN;
    else
        element = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
