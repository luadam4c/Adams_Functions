function rowConditions = m3ha_determine_row_conditions (rowMode, colMode, ...
                                        attemptNumber, gCondToFit, pCondToFit)
%% Determine the conditions for each row
% Usage: rowConditions = m3ha_determine_row_conditions (rowMode, colMode, ...
%                                       attemptNumber, gCondToFit, pCondToFit)
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
%       cd/m3ha_select_raw_traces.m

% File History:
% 2017-05-20 Moved from singleneuronfitting2.m
% 2017-05-22 Changed line width and indentation
% 2018-11-15 Moved to Adams_Functions
% 2018-11-15 Improved documentation

%% Hard-coded parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Count the number of possible conductance amplitude scaling percentages
nGIncr = length(gCondToFit);

% Count the number of possible pharm conditions 
nPCond = length(pCondToFit);

%% Do the job
if rowMode == 1 || ...
    (colMode == 1 && attemptNumber <= 2) || ...
    (colMode == 2 && attemptNumber <= 3)
    % Each row is one pharm condition
    rowConditions = pCondToFit;
elseif rowMode == 2
    % Each row is a pharm condition paired with a g incr
    nRows = nPCond * nGIncr;

    % Find each row condition
    rowConditions = zeros(nRows, 2);
    for iRow = 1:nRows
        % Find pharm condition of each row (index in pCondToFit)
        rowConditions(iRow, 1) = floor((iRow-1) / nGIncr) + 1;    

        % Find g-incr condition of each row (index in gCondToFit)
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

% The following must be consistent with dclampDataExtractor.m
gCondToFit = [100; 200; 400]; % possible conductance amplitude scaling percentages
pCondToFit = [1; 2; 3; 4];    % possible pharm conditions 
                            %   1 - Control
                            %   2 - GAT1 Block
                            %   3 - GAT3 Block
                            %   4 - Dual Block

%       cd/singleneuronfitting42.m and later versions

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%