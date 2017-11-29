function cellarray = vec2cell (vector, classvector, varargin)
%% Reorganize a vector or array into a cell array of partial vectors/arrays according to class
% Usage: cellarray = vec2cell (vector, classvector, varargin)
% Outputs:  
%       cellarray   - the cell array of vectors divided by class
%                   specified as a cell array of numeric vectors
% Arguments:
%       vector      - vector(s) to reorganize, each vector is a column
%                   must be a nonempty array
%       classvector - vector of classes of each element
%                   must be a numeric vector with the same length as vector
%
% Requires:
%
% Used by:    
%        /media/adamX/RTCl/tuning_curves.m
%
% File History:
% 2017-05-12 Created by Adam Lu
% 2017-10-31 Now accepts an array (each column is a vector to be separated) 
%               as the first argument
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error('Not enough input arguments, type ''help vec2cell'' for usage');
end

% Add required inputs to an Input Parser
iP = inputParser;
addRequired(iP, 'vector', ...           % vector(s) to reorganize
    @(x) validateattributes(x, {'numeric', 'cell', 'struct'}, {'nonempty'}));
addRequired(iP, 'classvector', ...      % vector of classes of each element
    @(x) validateattributes(x, {'numeric'}, {'vector', 'nonempty'}));

% Add parameter-value pairs to the Input Parser

% Read from the Input Parser
parse(iP, vector, classvector, varargin{:});

% Check relationships between arguments
if size(vector, 1) ~= length(classvector)
    error(['classvector length must be the same as ', ...
            'the number of rows in the array!']);
end

%% Find the unique classes 
uniqueclasses = unique(classvector);    
                        % this also puts the classes in ascending order 
                        %% TODO: 'Order' == 'original'
numclasses = numel(uniqueclasses);

%% Create cell array of vectors for each class
%% TODO: Not sure if this works for cell arrays and struct arrays
cellarray = cell(1, numclasses);
for c = 1:numclasses
    cellarray{c} = vector(classvector == uniqueclasses(c), :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
