function vector = piecelinspace (nodes, nPoints)
%% Generates a piece-wise linear vector from nodes and number of points
% Usage: vector = piecelinspace (nodes, nPoints)
%
% Used by:
%       /media/adamX/m3ha/data_dclamp/CountSweeps.m
%
% File History:
%   2018-08-02 - Created by Adam Lu

% Count the number of nodes
nNodes = length(nodes);

% Make sure the number of points are at least the number of nodes
if nPoints < nNodes
    fprintf('The number of points can''t be smaller than the number of nodes!');
    return
end

% Compute the number of segments
nSegments = nNodes - 1;

% Compute the number of points not on a node (these will be interpolated)
nPointsToInterpolate = nPoints - nNodes;

% Partition the points to interpolate into nSegments bins
nPointsToInterpolateEachSegment = ...
    distribute(nPointsToInterpolate, nSegments);

%