function ranks = convert_to_rank (array, varargin)
%% Creates a positive integer array from an array showing the ranks of each element
% Usage: ranks = convert_to_rank (array, varargin)
% Explanation:
%       TODO
% Example(s):
%       convert_to_rank([20, 3, 3, 10, 10, 20, 3])
%       convert_to_rank({'dog'; 'cat'; 'cat'; 'dog'})
%       convert_to_rank(["dog", "cat", "cat", "dog"])
%       convert_to_rank({'dog'; 'cat'}, 'Ranked', ["dog"; "cat"])
%       convert_to_rank(["dog", "fly", "cat"], 'Ranked', {'dog'; 'cat'})
% Outputs:
%       ranks        - ranks of each element
%                   specified as a positive integer array
% Arguments:
%       array       - an array that can be passed into unique()
%       varargin    - 'RankedElements': ranked elements
%                   must be an array of the same type as array
%                   default == unique(array)
%                   - Any other parameter-value pair for the unique() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/find_ind_str_in_cell.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/plot_grouped_histogram.m

% File History:
% 2019-01-09 Created by Adam Lu
% 

%% Hard-coded parameters
validSearchModes = {'exact', 'substrings', 'regexp'};

%% Default values for optional arguments
rankedElementsDefault = [];         % set later
searchModeDefault = 'exact';        % search exact matches by default
ignoreCaseDefault = false;          % don't ignore case by default

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
addRequired(iP, 'array');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'RankedElements', rankedElementsDefault);
addParameter(iP, 'SearchMode', searchModeDefault, ...   % the search mode
    @(x) any(validatestring(x, validSearchModes)));
addParameter(iP, 'IgnoreCase', ignoreCaseDefault, ...   % whether to ignore case
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, array, varargin{:});
rankedElements = iP.Results.RankedElements;
searchMode = validatestring(iP.Results.SearchMode, validSearchModes);
ignoreCase = iP.Results.IgnoreCase;

% Keep unmatched arguments for the unique() function
otherArguments = iP.Unmatched;
% otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Set default ranks of elements
if isempty(rankedElements)
    rankedElements = unique(array, otherArguments);
end

% Check type compatibility of array and rankedElements
% TODO: cellstr and string should be made compatible
% if ~strcmp(class(array), class(rankedElements))
%     error('array and rankedElements do not have the same type!!')
% end

%% Do the job
% Find the 
if iscellstr(array) || isstring(array)
    if iscell(array)
        % Find the index in rankedElements for each element in the cell array
        ranks = cellfun(@(x) find_ind_str_in_cell(x, rankedElements, ...
                                'MaxNum', 1, 'ReturnNan', true, ...
                                'SearchMode', searchMode, ...
                                'IgnoreCase', ignoreCase), array);
    else
        % Find the index in rankedElements for each element in the string array
        ranks = arrayfun(@(x) find_ind_str_in_cell(x, rankedElements, ...
                                'MaxNum', 1, 'ReturnNan', true, ...
                                'SearchMode', searchMode, ...
                                'IgnoreCase', ignoreCase), array);
    end
else
    % Find the indices in rankedElements for each element in array
    %   Note: If not found, zero will be returned
    [~, ranks] = ismember(array, rankedElements);
end

% Make all zeros NaNs instead
ranks(ranks == 0) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%