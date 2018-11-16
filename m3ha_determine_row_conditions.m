function rowConditions = ...
                m3ha_determine_row_conditions (rowMode, colMode, attemptNumber)
%% Determine the conditions for each row
% Usage: rowConditions = ...
%               m3ha_determine_row_conditions (rowMode, colMode, attemptNumber)
% 
% Explanation: 
%       TODO
% Example(s):
%       TODO
% Outputs:
%       TODO
%
% Arguments:
%       TODO
%
% Requires:
%
% Used by:
%       cd/singleneuronfitting42.m and later versions

% File History:
% 2017-05-20 Moved from singleneuronfitting2.m
% 2017-05-22 Changed line width and indentation
% 2018-11-15 Moved to Adams_Functions
% 2018-11-15 Improved documentation

%% Hard-coded parameters
% The following must be consistent with dclampDataExtractor.m
gIncrAll = [100; 200; 400]; % possible conductance amplitude scaling percentages
pCondAll = [1; 2; 3; 4];    % possible pharm conditions 
                            %   1 - Control
                            %   2 - GAT1 Block
                            %   3 - GAT3 Block
                            %   4 - Dual Block

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Count the number of possible conductance amplitude scaling percentages
nGIncr = length(gIncrAll);

% Count the number of possible pharm conditions 
nPCond = length(pCondAll);

%% Do the job
if rowMode == 1 || ...
    (colMode == 1 && attemptNumber <= 2) || ...
    (colMode == 2 && attemptNumber <= 3)
    % Each row is one pharm condition
    rowConditions = pCondAll;            % pharm condition of each row
elseif rowMode == 2
    % Each row is a pharm condition paired with a g incr
    nRows = nPCond * nGIncr;

    % Find each row condition
    rowConditions = zeros(nRows, 2);
    for iRow = 1:nRows
        % Find pharm condition of each row (index in pCondAll)
        rowConditions(iRow, 1) = floor((iRow-1) / nGIncr) + 1;    

        % Find g-incr condition of each row (index in gIncrAll)
        rowConditions(iRow, 2) = mod(iRow - 1, nGIncr) + 1;        
    end
else
    error('row mode undefined!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

global outparams
rowMode = 1;

    nRows = nPCond;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%