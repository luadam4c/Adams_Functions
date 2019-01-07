function [hasSameAttributes, rowIndices] = ...
                has_same_attributes (table, rowNames, attributes, varargin)
%% Returns all row indices that have the same attributes (column values) combination as a given set of row names
% Usage: [hasSameAttributes, rowIndices] = ...
%               has_same_attributes (table, rowNames, attributes, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       hasSameAttributes   - whether each row has matching attributes
%                           specified as a logical vector
%       rowIndices  - row indices with matching attributes
%                   specified as a positive integer vector
% Arguments:
%       table       - a table
%                   must be a 2D table
%       rowNames    - row names to match attributes
%                   must be a cell array of character vectors or a string array
%       attributes  - variable names for the attributes to match
%                   must be a cell array of character vectors or a string array
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
%
% Used by:
%       cd/m3ha_select_sweeps_to_fit.m

% File History:
% 2018-11-19 Adapted from code in m3ha_select_sweeps_to_fit.m
% 2018-12-06 Changed find_rows_with_same_attributes() -> has_same_attributes()
% 2018-12-06 Now outputs a logical vector
% 

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
addRequired(iP, 'rowNames', ...
    @(x) assert(iscellstr(x) || isstring(x), ...
                ['rowNames must be a cell array of character arrays', ...
                'or a string array!']));
addRequired(iP, 'attributes', ...
    @(x) assert(iscellstr(x) || isstring(x), ...
                ['attributes must be a cell array of character arrays', ...
                'or a string array!']));

% Read from the Input Parser
parse(iP, table, rowNames, attributes, varargin{:});

%% Preparation
% Count the total number of rows
nRows = height(table);

% Extract all the row names
allRowNames = table.Properties.RowNames;

% Extract just the attributes of interest
attributesOfInterest = table{:, attributes};

% Initialize all rows to not have the same attributes
hasSameAttributes = false(nRows, 1);

%% Do the job
% Loop through all rows but only doing any detection if needed
for iRow = 1:nRows
    % If this row is already included, skip it
    if hasSameAttributes(iRow)
        continue
    end

    % Get the current row name
    rowNameThis = allRowNames(iRow);

    % If this row is not one to match attributes for, skip it
    if ~ismember(rowNameThis, rowNames)
        continue
    end

    % Get the attribute values for this row
    attributesThis = attributesOfInterest(iRow, :);

    % Test whether each row has the same attribute values as this row
    hasSameAttributesToThisRow = ...
        ismember(attributesOfInterest, attributesThis, 'rows');

    % Update all rows with the same attributes values as this row
    hasSameAttributes = hasSameAttributes | hasSameAttributesToThisRow;
end

% Find the corresponding row indices
rowIndices = find(hasSameAttributes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Extract just the attributes of interest
tableOfInterest = table(:, attributes);

% Get the attribute values
attributesThis = tableOfInterest{rowNameThis, :};

rowfun(@(x) isequal(x, attributesThis), tableOfInterest);

% Initialize a vector for storing row indices
rowIndices = [];

if ismember(iRow, rowIndices)
    continue
end

% Find the rows with the same attribute values
indSameAttributes = find(hasSameAttributesToThisRow);

% Add to rowIndices
rowIndices = [rowIndices; indSameAttributes];

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%