function vector = piecelinspace (nodes, nPoints, varargin)
%% Generates a piece-wise linear vector from nodes and number of points
% Usage: vector = piecelinspace (nodes, nPoints, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       piecelinspace([4, 1, 6], 10)
%       piecelinspace([4; 1; 6], 10)
%
% Outputs:    
%       vector      - vector spanning nodes
%                   must be a numeric vector
%
% Arguments:    
%       nodes       - nodes that must be included
%                   must be a numeric vector
%       nPoints     - total number of points in the vector
%                   must be a positive integer scalar
%
% Requires:
%       cd/distribute_balls_into_boxes.m
%       cd/force_row_vector.m
%
% Used by:
%       cd/m3ha_compute_gabab_ipsc.m
%       cd/m3ha_parse_sweep_settings.m
%       /media/adamX/m3ha/data_dclamp/CountSweeps.m

% File History:
%   2018-08-06 - Created by Adam Lu
%   2020-01-05 - Now returns a column vector if nodes is a column vector
% TODO: Add 'ForceRowOutput' as an optional argument

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'nodes', ...        % nodes that must be included
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'nPoints', ...      % total number of points in the vector
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Read from the Input Parser
parse(iP, nodes, nPoints, varargin{:});

%% Preparation
% Save row or column status and force as a row vector
if iscolumn(nodes)
    wasColumn = true;
    nodes = force_row_vector(nodes);
else
    wasColumn = false;
end

% Count the number of nodes
nNodes = numel(nodes);

% Make sure the number of points are at least the number of nodes
if nPoints < nNodes
    fprintf('The number of points can''t be smaller than the number of nodes!');
    return
end

%% Generate the vector
% If there is only one node, generate a vector with the same value repeated
%   nPoints times
if nNodes == 1
    vector = nodes(1) * ones(1, nPoints);
    return;
end

% Compute the number of segments
nSegments = nNodes - 1;

% Compute the number of points not on a node (these will be interpolated)
nPointsToInterpolate = nPoints - nNodes;

% Compute the number of points for interpolating each segment
nPointsEachSegment = ...
    distribute_balls_into_boxes(nPointsToInterpolate, nSegments, ...
                                'Evenly', true) + 2;

% Generate linearly interpolated vectors for each segment in a cell array
allSegments = arrayfun(@(x, y, z) linspace(x, y, z), ...
                      nodes(1:end-1), nodes(2:end), nPointsEachSegment, ...
                      'UniformOutput', false);

% Take away the last point of all segments
allSegmentsTrimmed = cellfun(@(x) x(1:end-1), allSegments, ...
                             'UniformOutput', false);

% Concatenate the vectors in the cell array
vectorWithoutEnd = cat(2, allSegmentsTrimmed{:});

% Add the very last node back
vector = [vectorWithoutEnd, nodes(end)];

% Force as a column vector
if wasColumn
    vector = vector(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%