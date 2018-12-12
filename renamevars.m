function table = renamevars (table, prevNames, newNames, varargin)
%% Rename variable(s) in a table
% Usage: table = renamevars (table, prevNames, newNames, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       table       - new table
%                   specified as a 2-D table
% Arguments:
%       table       - old table
%                   must be a 2-D table
%       prevNames   - previous variable names
%                   must be a string/character array or 
%                       a cell array of strings/character arrays
%       newNames    - new variable names
%                   must be a string/character array or 
%                       a cell array of strings/character arrays
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
%                   - 'MaxNum': maximum number of indices to find
%                   must be a positive integer
%                   default == numel(cellArray)
%
% Requires:
%       cd/find_ind_str_in_cell.m
%       cd/ispositiveintegerscalar.m
%
% Used by:
%       cd/load_params.m

% File History:
% 2018-12-12 Created by Adam Lu

%% Hard-coded constants
validSearchModes = {'exact', 'substrings', 'regexp'};

%% Default values for optional arguments
searchModeDefault = 'substrings';       % default search mode
ignoreCaseDefault = false;              % whether to ignore case by default
maxNumDefault = [];                     % will be changed to numel(cellArray)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'table', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addRequired(iP, 'prevNames', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
                ['prevNames must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));
addRequired(iP, 'newNames', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
                ['prevNames must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SearchMode', searchModeDefault, ...   % the search mode
    @(x) any(validatestring(x, validSearchModes)));
addParameter(iP, 'IgnoreCase', ignoreCaseDefault, ...   % whether to ignore case
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MaxNum', maxNumDefault, ...       % maximum number of indices
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'MaxNum must be either empty or a positive integer scalar!'));

% Read from the Input Parser
parse(iP, table, prevNames, newNames, varargin{:});
searchMode = validatestring(iP.Results.SearchMode, validSearchModes);
ignoreCase = iP.Results.IgnoreCase;
maxNum = iP.Results.MaxNum;

% Check relationships between arguments
% TODO: prevNames and newNames must be compatible

%% Preparation
% Match the row counts of prevNames and newNames
%   and force them into column cell arrays
[prevNames, newNames] = ...
    match_format_vectors(prevNames, newNames, 'ForceCellOutputs', true);

%% Do the job
% Get all the variable names
variableNames = table.Properties.VariableNames;

% Loop through all variables to change
for iVar = 1:numel(prevNames)
    % Get this previous name
    prevName = prevNames{iVar};

    % Find the matching variable number(s)
    varNumber = find_ind_str_in_cell(prevName, variableNames, ...
                                  'SearchMode', searchMode, ...
                                  'IgnoreCase', ignoreCase, 'MaxNum', maxNum);

    % Rename the column
    if ~isempty(varNumber)
        % Get this new name as an element
        newName = newNames{iVar};

        % Count the number of matching variables
        nMatches = length(varNumber)

        % Construct variable names
        if nMatches > 1
            % If there are more than one matching variables
            %   append an iteration to avoid columns with the same names
            newNameCell = arrayfun(@(x) [newName, '_', num2str(x)], ...
                                    transpose(1:nMatches), ...
                                    'UniformOutput', false);
        else
            % Just use the new variable name
            newNameCell = {newName};
        end

        % Change the variable names
        table.Properties.VariableNames(varNumber) = newNameCell;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Force prevNames and newNames into column cell arrays
[prevNames, newNames] = ...
    argfun(@force_column_cell, prevNames, newNames);

newNameCell = newNames(iVar);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%