function index = find_index_in_array (element, array, varargin)
%% Returns the index for element(s) in an array
% Usage: index = find_index_in_array (element, array, varargin)
% Explanation:
%       TODO
% Example(s):
%       find_index_in_array({'dog'; 'cat'}, ["dog"; "fly"; "cat"])
%       find_index_in_array(["dog", "fly", "cat"], {'dog'; 'cat'})
% Outputs:
%       index       - the index in an array with matching element
%                   specified as a positive integer array (may contain NaN)
% Arguments:
%       element     - element(s) of an array
%       array       - an array
%       varargin    - 'SearchMode': the search mode
%                   must be an unambiguous, case-insensitive match to one of:
%                       'exact'         - str must be identical to 
%                                           an element in cellArray
%                       'substrings'    - str can be a substring or 
%                                           a cell array of substrings
%                       'regexp'        - str is considered a regular expression
%                   if searchMode is 'exact' or 'regexp', 
%                       str cannot be a cell array
%                   default == 'substrings'
%                   - 'IgnoreCase': whether to ignore differences in letter case
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/find_ind_str_in_cell.m
%
% Used by:
%       cd/convert_to_rank.m

% File History:
% 2019-01-09 Created by Adam Lu
% 

%% Hard-coded parameters
validSearchModes = {'exact', 'substrings', 'regexp'};

%% Default values for optional arguments
searchModeDefault = 'exact';        % search exact matches by default
ignoreCaseDefault = false;          % don't ignore case by default

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
addRequired(iP, 'element');
addRequired(iP, 'array');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SearchMode', searchModeDefault, ...   % the search mode
    @(x) any(validatestring(x, validSearchModes)));
addParameter(iP, 'IgnoreCase', ignoreCaseDefault, ...   % whether to ignore case
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, element, array, varargin{:});
searchMode = validatestring(iP.Results.SearchMode, validSearchModes);
ignoreCase = iP.Results.IgnoreCase;

%% Preparation
% Check type compatibility of element and array
% TODO: cellstr and string should be made compatible
% if ~strcmp(class(element), class(array))
%     error('element and array do not have the same type!!')
% end

%% Do the job
if ischar(element) || iscellstr(element) || isstring(element)
    % Make sure element is either a string or a cell array
    if ischar(element)
        element = {element};
    end

    % Use either cellfun or arrayfun
    if iscell(element)
        % Find the index in array for each element in element
        index = cellfun(@(x) find_ind_str_in_cell(x, array, ...
                                'MaxNum', 1, 'ReturnNan', true, ...
                                'SearchMode', searchMode, ...
                                'IgnoreCase', ignoreCase), element);
    else
        % Find the index in array for each element in element
        index = arrayfun(@(x) find_ind_str_in_cell(x, array, ...
                                'MaxNum', 1, 'ReturnNan', true, ...
                                'SearchMode', searchMode, ...
                                'IgnoreCase', ignoreCase), element);
    end
else
    % Find the index in array for each element in element
    %   Note: If not found, zero will be returned
    [~, index] = ismember(element, array);

    % Make all zeros NaNs instead
    index(index == 0) = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%