function unionOfContents = union_over_cells (cellArray, varargin)
%% Apply the union function over all contents of a cell array
% Usage: unionOfContents = union_over_cells (cellArray, varargin)
% Outputs:    
%       unionOfContents - union of contents
%
% Arguments:
%       cellArray       - a cell array of arrays that will be unioned;
%                           if just an array, return the array
%                       must be a cell array of input arrays that 
%                           can be recognized by the built-in union function
%
%
% Used by:
%   /media/adamX/m3ha/optimizer4gabab/compare_and_plot_across_conditions.m
%
% File History:
%   2018-08-17 Created by Adam Lu

%% Hard-coded parameters
validSetOrders = {'sorted', 'stable'};

%% Default values for optional arguments
setOrderDefault = 'sorted';

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

% TODO: Add required inputs to an Input Parser
% Add required inputs to the Input Parser
addRequired(iP, 'cellArray', @iscell)

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SetOrder', setOrderDefault, ...
    @(x) any(validatestring(x, validSetOrders)));

% Read from the Input Parser
parse(iP, cellArray, varargin{:});
setOrder = validatestring(iP.Results.SetOrder, validSetOrders);

if isempty(cellArray)
    error('Argument cannot be empty!');
end

%% Return the array if not a cell array
if ~iscell(cellArray)
    unionOfContents = cellArray;
    return;
end

%% Return the union of the contents of the cell array
% Count the number of elements
nElements = numel(cellArray);

% If there are no elements, return an empty matrix
if nElements == 0
    unionOfContents = [];
    return
end

% Initialize the union
unionOfContents = cellArray{1};

% Iterate over all elements
if nElements > 1
    for iElement = 2:nElements
        unionOfContents = union(unionOfContents, cellArray{iElement}, setOrder);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}