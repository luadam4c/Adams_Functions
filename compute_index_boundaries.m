function boundaryValues = compute_index_boundaries (varargin)
%% Computes boundary values for indices of different groups, assuming the groups are all consecutive in the array
% Usage: boundaryValues = compute_index_boundaries (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       compute_index_boundaries('NEachGroup', 1)
%       compute_index_boundaries('Grouping', [1 1 1 2 2 3 3 3])
%       compute_index_boundaries('Grouping', [NaN NaN 1 1 NaN NaN 2 2 NaN NaN 3 3 NaN NaN])
%       compute_index_boundaries('Grouping', [NaN NaN 1 1 NaN NaN 2 2 NaN NaN 3 3 NaN NaN], 'TreatNaNsAsGroup', false)
%       compute_index_boundaries('Grouping', {[1 1 1 2 2 3 3 3], [2; 2; 4; 4]})
%       compute_index_boundaries('NEachGroup', [3, 2, 3])
%
% Outputs:
%       boundaryValues  - boundary values
%                       specified as a numeric vector
% Arguments:
%       varargin    - 'Grouping': vector with the same value for each group
%                   must be an array
%                   default == []
%                   - 'NEachGroup': number in each group
%                   must be a positive integer vector or a cell array of them
%                   default == []
%                   - 'TreatNaNsAsGroup': whether to treat NaNs 
%                                           as a distinct group
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for the unique_groups() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%       cd/unique_groups.m
%
% Used by:
%       cd/combine_abf_data.m
%       cd/plot_tuning_curve.m

% File History:
% 2019-08-21 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
groupingDefault = [];
nEachGroupDefault = [];
treatNaNsAsGroupDefault = true;       % ignore NaN by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Grouping', groupingDefault);
addParameter(iP, 'NEachGroup', nEachGroupDefault);
addParameter(iP, 'TreatNaNsAsGroup', treatNaNsAsGroupDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'scalar'}));

% Read from the Input Parser
parse(iP, varargin{:});
grouping = iP.Results.Grouping;
nEachGroup = iP.Results.NEachGroup;
treatNaNsAsGroup = iP.Results.TreatNaNsAsGroup;

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

%% Do the job
% Count the number in each group, preserving the order of groups
if isempty(nEachGroup) || (~treatNaNsAsGroup && ~isempty(grouping))
    % Count the number in each group, preserving the order of groups
    [uniqueGroups, nEachGroup] = ...
        unique_groups(grouping, 'PreserveOrder', true, otherArguments{:});
else
    uniqueGroups = [];
end

% Compute index boundary values
if iscellnumeric(nEachGroup)
    boundaryValues = ...
        cellfun(@(x, y) compute_boundaries_helper(x, y, treatNaNsAsGroup), ...
                uniqueGroups, nEachGroup, 'UniformOutput', false);
else
    boundaryValues = compute_boundaries_helper (uniqueGroups, ...
                                                nEachGroup, treatNaNsAsGroup);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function boundaryValues = compute_boundaries_helper (uniqueGroups, ...
                                                nEachGroup, treatNaNsAsGroup)

% Count the number of groups
nGroups = numel(nEachGroup);

% Count the number of boundaries
nBoundaries = nGroups - 1;

% Get the index of the last sweep for each phase
iLastEachGroup = cumsum(nEachGroup);

% Compute the phase boundaries
boundaryValues = iLastEachGroup(1:nBoundaries) + 0.5;

% Remove boundaries between NaNs and actual groups if requested
if ~treatNaNsAsGroup && ~iscell(uniqueGroups)
    % Test whether each group is NaN
    isNaNGroup = isnan(uniqueGroups);

    % Remove the last boundary if it is followed by a group of NaNs
    if isNaNGroup(end)
        boundaryValues(end) = [];
    end

    % If there is anything left, remove any boundary preceded by a group of NaNs
    if length(boundaryValues) > 1
        boundaryValues(isNaNGroup(1:end-1)) = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%