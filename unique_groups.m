function varargout = unique_groups (grouping, varargin)
%% Retrieves the unique groups and counts the number in each group
% Usage: [uniqueGroups, nEachGroup] = unique_groups (grouping, setOrder (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [uG, nEG] = unique_groups([3; 3; 3; 2; 2; 1])
%       [uG, nEG] = unique_groups([3; 3; 3; 2; 2; 1], 'sorted')
%       [uG, nEG] = unique_groups([1 1 1 2 2 3 3 3])
%       [uG, nEG] = unique_groups([1 1 2 2 3 3 1 1])
%       [uG, nEG] = unique_groups([1 1 2 2 3 3 1 1], 'PreserveOrder', true)
%       [uG, nEG] = unique_groups([1 NaN 2 1 1 NaN 3 NaN 2 3], 'IgnoreNaN', true)
%    TODO:   [uG, nEG] = unique_groups([1 NaN 2 1 1 NaN 3 NaN 2 3], 'TreatNaNAsEqual', false)
%       [uG, nEG] = unique_groups([NaN NaN 2 2 1 1 NaN NaN 2 2 NaN NaN])
%       [uG, nEG] = unique_groups([NaN NaN 2 2 1 1 NaN NaN 2 2 NaN NaN], 'PreserveOrder', true)
%       [uG, nEG] = unique_groups({[1 1 1 2 2 3 3 3], [2; 4; 4; 2]})
%
% Outputs:
%       uniqueGroups    - the unique groups
%                       specified as an array with the same type as grouping
%       nEachGroup      - number in each group
%                       specified as a numeric vector
%
% Arguments:
%       grouping    - a vector with the same value for all elements of each group
%                   must be an array
%       setOrder    - (opt) setOrder for the unique() function
%                   must be an unambiguous, case-insensitive match to one of: 
%                       - 'sorted'
%                       - 'stable'
%                   default == 'stable'
%       varargin    - 'PreserveOrder': whether to preserve order of groups
%                       Note: this will cause nonconsecutive groups
%                               to become separate groups
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'IgnoreNaN': whether to ignore NaNs
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatNanAsEqual': whether to treat all NaN values
%                                           as the same
%                       Note: If 'IgnoreNaN' == true, 
%                           'TreatNanAsEqual' has no effect
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for the unique() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_vector.m
%       cd/force_row_vector.m
%       cd/iscellnumeric.m
%       cd/ismatch.m
%       cd/struct2arglist.m
%       cd/unique_custom.m
%
% Used by:
%       cd/test_difference.m
%       cd/test_var_difference.m

% File History:
% 2019-08-21 Created by Adam Lu
% 

%% Hard-coded parameters
validSetOrders = {'sorted', 'stable'};

%% Default values for optional arguments
setOrderDefault = 'stable';
preserveOrderDefault = false;   % do not preserve order by default
ignoreNanDefault = false;       % do not ignore NaN by default
treatNanAsEqualDefault = true;  % treat all NaN values equal by default

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

% Add optional inputs to the Input Parser
addOptional(iP, 'setOrder', setOrderDefault, ...
    @(x) any(validatestring(x, validSetOrders)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PreserveOrder', preserveOrderDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'scalar'}));
addParameter(iP, 'IgnoreNaN', ignoreNanDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'scalar'}));
addParameter(iP, 'TreatNanAsEqual', treatNanAsEqualDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'scalar'}));

% Read from the Input Parser
parse(iP, grouping, varargin{:});
setOrder = validatestring(iP.Results.setOrder, validSetOrders);
preserveOrder = iP.Results.PreserveOrder;
ignoreNan = iP.Results.IgnoreNaN;
treatNanAsEqual = iP.Results.TreatNanAsEqual;

% Keep unmatched arguments for the unique() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
if isnumeric(grouping) && ~isrow(grouping) && ~iscolumn(grouping)
    grouping = force_column_vector(grouping);
end

%% Do the job
% Do this recursively if grouping is a cell array of numeric vectors
if iscellnumeric(grouping)
    [uniqueGroups, nEachGroup] = ...
        cellfun(@(x) unique_groups_helper(x, setOrder, preserveOrder, ...
                                ignoreNan, treatNanAsEqual, otherArguments), ...
                grouping, 'UniformOutput', false);
else
    [uniqueGroups, nEachGroup] = ...
        unique_groups_helper(grouping, setOrder, preserveOrder, ...
                            ignoreNan, treatNanAsEqual, otherArguments);
end

%% Output results
varargout{1} = uniqueGroups;
if nargout >= 2
    varargout{2} = nEachGroup;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = unique_groups_helper(grouping, setOrder, preserveOrder, ...
                                    ignoreNan, treatNanAsEqual, otherArguments)

% Get the unique groups
if preserveOrder
    % Remember whether column or row
    if isrow(grouping)
        isRow = true;
    else
        isRow = false;
    end

    % Force as a column vector
    grouping = grouping(:);

    if ignoreNan
        % Remove NaNs from the data
        grouping = grouping(~isnan(grouping));
    end

    % Count the number of elements
    nElements = numel(grouping);

    if nElements == 0
        if iscell(grouping)
            uniqueGroups = {};
        else
            uniqueGroups = [];
        end
        nEachGroup = 0;
    elseif nElements == 1
        uniqueGroups = grouping;
        nEachGroup = 1;
    else
        % Create the before vector
        befores = grouping(1:end-1);

        % Create the after vector
        afters = grouping(2:end);

        % Decide whether each index starting from the 2nd is to be removed
        if iscell(grouping)
            toRemove2ToEnd = ...
                cellfun(@(x, y) isequal_to_use(x, y, treatNanAsEqual), ...
                        befores, afters);
        else
            toRemove2ToEnd = ...
                arrayfun(@(x, y) isequal_to_use(x, y, treatNanAsEqual), ...
                        befores, afters);
        end

        % Always keep the first index
        toKeep = [true; ~toRemove2ToEnd];

        % Find the indices to keep
        indToKeep = find(toKeep);

        % Compute the number in each group
        nEachGroup = diff([indToKeep; nElements]);

        % Get the unique groups
        uniqueGroups = grouping(toKeep);

        if isRow
            uniqueGroups = force_row_vector(uniqueGroups);
            nEachGroup = force_row_vector(nEachGroup);
        end
    end
else
    uniqueGroups = unique_custom(grouping, setOrder, 'IgnoreNan', ignoreNan, ...
                        'TreatNanAsEqual', treatNanAsEqual, otherArguments{:});

    % Count the number of elements in each group
    % TODO: This doesn't work when treatNanAsEqual is false
    if nargout >= 2
        if iscell(uniqueGroups)
            nEachGroup = cellfun(@(x) sum(ismatch(grouping, x)), uniqueGroups);
        else
            nEachGroup = arrayfun(@(x) sum(ismatch(grouping, x)), uniqueGroups);
        end
    end
end

%% Output results
varargout{1} = uniqueGroups;
if nargout >= 2
    varargout{2} = nEachGroup;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isEqual = isequal_to_use(x, y, treatNanAsEqual)

if treatNanAsEqual
    isEqual = isequaln(x, y);
else
    isEqual = isequal(x, y);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%