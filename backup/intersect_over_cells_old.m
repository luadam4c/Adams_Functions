function intersection = intersect_over_cells (cellArray, varargin)
%% Apply the intersect function over all contents of a cell array
% Usage: intersection = intersect_over_cells (cellArray, varargin)
% Explanation:
%       TODO
%
% Examples:
%       intersect_over_cells([2; 3])
%       intersect_over_cells({'ab', 'bc', 'bg'})
%       intersect_over_cells({[1 2], [2; 3]})
% Outputs:    
%       intersection    - intersection of contents
%
% Arguments:
%       cellArray       - a cell array of arrays that will be intersected;
%                           if just an array, return the array
%                       must be a cell array of input arrays that 
%                           can be recognized by the built-in intersect function
%
% Used by:
%       TODO

% File History:
%   2018-01-10 Created
%   2018-05-08 Fixed the case if cellArray is not a cell array
%   2018-05-08 Changed tabs to spaces and limited width to 80
%   2018-08-17 Renamed intersectm -> intersect_over_cells
%   2019-08-21 Now uses apply_over_cells.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'cellArray');

% Read from the Input Parser
parse(iP, cellArray, varargin{:});

%% Return the array if not a cell array
if ~iscell(cellArray)
    intersection = cellArray;
    return;
end

%% Return the intersection of the contents of the cell array
% Count the number of elements
nElements = numel(cellArray);

% If there are no elements, return an empty matrix
if nElements == 0
    intersection = [];
    return
end

% Initialize the intersection
intersection = cellArray{1};

% Iterate over all elements
if nElements > 1
    for iElement = 2:nElements
        intersection = intersect(intersection, cellArray{iElement}, setOrder);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%       cd/find_in_strings.m
if isempty(cellArray)
    error('Argument cannot be empty!');
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
