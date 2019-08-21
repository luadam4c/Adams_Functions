function boundaryValues = compute_index_boundaries (varargin)
%% Computes boundary values for indices of different groups
% Usage: boundaryValues = compute_index_boundaries (varargin)
% Explanation:
%       TODO
% Example(s):
%       compute_index_boundaries('Grouping', [1 1 1 2 2 3 3 3])
%       compute_index_boundaries('NEachGroup', [3, 2, 3])
% Outputs:
%       boundaryValues  - boundary values
%                       specified as a numeric vector
% Arguments:
%       varargin    - 'Grouping': vector with the same value for each group
%                   must be an array
%                   default == []
%                   - 'NEachGroup': number in each group
%                   must be a positive integer vector
%                   default == []
%                   - Any other parameter-value pair for the unique_groups() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%       cd/unique_groups.m
%
% Used by:
%       cd/combine_data_from_same_slice.m

% File History:
% 2019-08-21 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
groupingDefault = [];
nEachGroupDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Grouping', groupingDefault);
addParameter(iP, 'NEachGroup', nEachGroupDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive', 'integer'}));

% Read from the Input Parser
parse(iP, varargin{:});
grouping = iP.Results.Grouping;
nEachGroup = iP.Results.NEachGroup;

% Keep unmatched arguments for the unique_groups() function
otherArguments = struct2arglist(iP.Unmatched);

% Check relationships between arguments
if ~isempty(grouping) && ~isempty(nEachGroup)
    disp('Grouping will be ignored since NEachGroup is provided!');
    return
elseif isempty(grouping) && isempty(nEachGroup)
    disp('One of Grouping or NEachGroup must be provided!');
    return
end

%% Preparation
% Count the number in each group
if isempty(nEachGroup)
    [~, nEachGroup] = unique_groups(grouping, otherArguments{:});
end

% Count the number of groups
nGroups = numel(nEachGroup);

% Count the number of boundaries
nBoundaries = nGroups - 1;

%% Do the job
% Get the index of the last sweep for each phase
iLastEachGroup = cumsum(nEachGroup);

% Compute the phase boundaries
boundaryValues = iLastEachGroup(1:nBoundaries) + 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%