function stdErr = stderr(X)
% Calculate the standard error of the mean
% Usage: stdErr = stderr(X)
%
% Used by:
%
% 2017-12-14 Created

stdErr = std(X)./sqrt(length(X));
