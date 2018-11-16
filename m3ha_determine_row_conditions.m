function rowConditions = m3ha_determine_row_conditions (rowmode, colmode, attemptNumber, pp, ggNew)
%% Determine what each row would be
% Usage: rowConditions = m3ha_determine_row_conditions (rowmode, colmode, attemptNumber, pp, ggNew)
% 
% Used by:
%       cd/singleneuronfitting42.m and later versions
%
% File History:
% 2017-05-20 Moved from singleneuronfitting2.m
% 2017-05-22 Changed line width and indentation
% 2018-11-15 Moved to Adams_Functions

global outparams

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rowmode == 1 || ...
    (colmode == 1 && attemptNumber <= 2) || ...
    (colmode == 2 && attemptNumber <= 3)
    rowmode = 1;
    % Each row is one pharm condition
    rowConditions = pp;            % pharm condition of each row
    nRows = length(pp);
elseif rowmode == 2
    % Each row is a pharm condition paired with a g incr
    rowConditions = zeros(length(pp)*length(ggNew), 2);
    nRows = length(pp)*length(ggNew);
    for k = 1:nRows
        % Find pharm condition of each row (index in pp)
        rowConditions(k, 1) = floor((k-1)/length(ggNew)) + 1;    

        % Find g-incr condition of each row (index in ggNew)
        rowConditions(k, 2) = mod(k-1, length(ggNew)) + 1;        
    end
else
    error('row mode undefined!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%