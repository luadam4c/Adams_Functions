function vector = piecelinspace (nodes, nPoints)
%% Generates a piece-wise linear row vector from nodes and number of points
% Usage: vector = piecelinspace (nodes, nPoints)
%
% Requires:
%       /home/Matlab/Adams_Functions/distribute_balls_into_boxes.m
%
% Used by:
%       /media/adamX/m3ha/data_dclamp/CountSweeps.m

% File History:
%   2018-08-06 - Created by Adam Lu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:   

%}