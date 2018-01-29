function cellArray = vec2cell (vecORarray, classVector, varargin)
%% Reorganize a vector or array into a cell array of partial vectors/arrays according to class
% Usage: cellArray = vec2cell (vecORarray, classVector, varargin)
% Outputs:  
%       cellArray   - the cell array of vectors divided by class
%                   specified as a cell array of numeric vectors
% Arguments:
%       vecORarray  - vector(s) to reorganize, each vector is a column
%                   must be a nonempty array
%       classVector - vector of classes of each element
%                   must be a numeric vector with 
%                       the same number of rows as vecORarray
%       varargin    - 'SetOrder' - a flag indicating the 
%                       order of vectors/arrays in cellArray
%                   TODO: annotate and implement this
%
% Requires:
%
% Used by:    
%       /media/adamX/RTCl/tuning_curves.m
%       /media/adamX/m3ha/network_model/tuning_curves.m
%
% File History:
% 2017-05-12 Created by Adam Lu
% 2017-10-31 Now accepts an array (each column is a vector to be separated) 
%               as the first argument
% TODO: Add a parameter 'SetOrder' that can take the values 
%           'sorted' (ascending order), 'stable' (original order)
%           with default being 'sorted'
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error('Not enough input arguments, type ''help vec2cell'' for usage');
end

% Add required inputs to an Input Parser
iP = inputParser;
addRequired(iP, 'vecORarray', ...       % vector(s) to reorganize
    @(x) validateattributes(x, {'numeric', 'cell', 'struct'}, {'nonempty'}));
addRequired(iP, 'classVector', ...      % vector of classes of each element
    @(x) validateattributes(x, {'numeric'}, {'vector', 'nonempty'}));

% Add parameter-value pairs to the Input Parser

% Read from the Input Parser
parse(iP, vecORarray, classVector, varargin{:});

% Check relationships between arguments
if size(vecORarray, 1) ~= length(classVector)
    error(['classVector length must be the same as ', ...
            'the number of rows in the array!']);
end

%% Find the unique classes
uniqueClasses = unique(classVector);    
                        % this also puts the classes in ascending order
                        % TODO: allow for original order
                        %       Hint: see documentation for unique() 
numClasses = numel(uniqueClasses);

%% Create cell array of vectors for each class
%% TODO: Test if this works for cell arrays and struct arrays
cellArray = cell(1, numClasses);
for c = 1:numClasses
    cellArray{c} = vecORarray(classVector == uniqueClasses(c), :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
