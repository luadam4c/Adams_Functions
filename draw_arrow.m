function draw_arrow(p1, p2)
%% Draw an arrow from p1 to p2
% USAGE: draw_arrow(p1, p2)
% Arguments:    p1, p2 must be 2x1 numeric arrays
%
% Used by:
%
%
% 2017-02-06 Moved from /media/adamX/Computational\ Neuroscience/Week\ 7/Week\ 7\ Quiz/Quiz_7_Pr_5.m (not really)
% 2017-05-21 Renamed drawArrow() -> draw_arrow()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot a velocity vector without automatic scaling
quiver(p1(1), p1(2), p2(1)-p1(1), p2(2)-p1(2), 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%