function plot_arrow (pointFrom, pointTo)
%% Draws an arrow from pointFrom to pointTo
% Usage: plot_arrow (pointFrom, pointTo)
% Arguments:
%       pointFrom   - coordinates for origin point
%                   must be a numeric vector with 2 elements
%       pointTo     - coordinates for destination point
%                   must be a numeric vector with 2 elements
%
% Used by:

% File History:
% 2017-02-06 Moved from /media/adamX/Computational\ Neuroscience/
%                       Week\ 7/Week\ 7\ Quiz/Quiz_7_Pr_5.m (not really)
% 2020-05-31 Renamed as plot_arrow.m
% TODO: Input Parser

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot a velocity vector without automatic scaling
quiver(pointFrom(1), pointFrom(2), ...
        pointTo(1) - pointFrom(1), pointTo(2) - pointFrom(2), 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%