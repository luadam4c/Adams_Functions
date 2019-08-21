function varargout = unique_groups (grouping, varargin)
%% Retrieves the unique groups and counts the number in each group
% Usage: [uniqueGroups, nEachGroup] = unique_groups (grouping, varargin)
% Explanation:
%       TODO
% Example(s):
%       [uG, nEG] = unique_groups([1 1 1 2 2 3 3 3])
%       [uG, nEG] = unique_groups([1 NaN 2 1 1 NaN 3 NaN 2 3])
%       [uG, nEG] = unique_groups([1 NaN 2 1 1 NaN 3 NaN 2 3], 'IgnoreNaN', true)
% Outputs:
%       uniqueGroups    - the unique groups
%                       specified as an array with the same type as grouping
%       nEachGroup      - number in each group
%                       specified as a numeric vector
% Arguments:
%       grouping    - a vector with the same value for all elements of each group
%                   must be an array
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the unique_custom() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%       cd/ismatch.m
%       cd/unique_custom.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-08-21 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'grouping');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, grouping, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the unique_custom() function
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
% Get the unique groups
uniqueGroups = unique_custom(grouping, 'stable', otherArguments{:});

% Count the number of elements in each group
if nargout >= 2
    if iscell(uniqueGroups)
        nEachGroup = cellfun(@(x) sum(ismatch(grouping, x)), uniqueGroups);
    else
        nEachGroup = arrayfun(@(x) sum(ismatch(grouping, x)), uniqueGroups);
    end
end

%% Output results
varargout{1} = uniqueGroups;
if nargout >= 2
    varargout{2} = nEachGroup;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%