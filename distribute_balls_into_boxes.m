function [allWays, nWays] = distribute_balls_into_boxes (nBalls, varargin)
%% Returns the ways and number of ways to distribute identical/discrete balls into identical/discrete boxes
% Usage: [allWays, nWays] = distribute_balls_into_boxes (nBalls, varargin)
%
% Outputs:
%       allWays     - all ways of distribution
%                   specified as a numeric array with nBoxes columns 
%                       (each row is a way of distribution)
%       nWays       - number of ways of distribution
%                   specified as a positive integer scalar
% Arguments:    
%       nBalls      - number of balls to place in boxes
%                   must be a positive integer scalar
%       nBoxes      - (opt) number of boxes to accept balls
%                   must be a positive integer scalar
%       varargin    - 'DistinctBalls': whether balls are distinct
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'DistinctBoxes': whether boxes are distinct
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'NonemptyBoxes': whether boxes must be nonempty
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Evenly': whether to distribute as evenly as possible
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/isemptycell.m
%
% Used by:
%       cd/load_examples.m
%       cd/piecelinspace.m

% File History:
%   2018-08-02 - Created by Adam Lu
%   TODO: Implement the cases where balls and/or boxes are distinct

%% Hard-coded parameters

%% Default values for optional arguments
nBoxesDefault = nBalls;             % default number of boxes is the number of balls
distinctBallsDefault = false;       % balls are identical by default
distinctBoxesDefault = false;       % boxes are identical by default
nonemptyBoxesDefault = false;       % boxes can be empty by default
evenlyDefault = false;              % not always evenly distributed by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'nBalls', ...                  % number of balls
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'nBoxes', nBoxesDefault, ...   % number of boxes
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'DistinctBalls', distinctBallsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'DistinctBoxes', distinctBoxesDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'NonemptyBoxes', nonemptyBoxesDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Evenly', evenlyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, nBalls, varargin{:});
nBoxes = iP.Results.nBoxes;
distinctBalls = iP.Results.DistinctBalls;
distinctBoxes = iP.Results.DistinctBoxes;
nonemptyBoxes = iP.Results.NonemptyBoxes;
evenly = iP.Results.Evenly;

%% Check if the distribution is possible
if nonemptyBoxes && nBalls < nBoxes
    fprintf('It is impossible to put %d balls into %d nonempy boxes!\n', ...
            nBalls, nBoxes);
    allways = {};
    nWays = 0;
    return;
end

%% Distribute balls!
if distinctBalls && distinctBoxes && nonemptyBoxes && evenly
    error('Not implemented yet!');
elseif distinctBalls && distinctBoxes && nonemptyBoxes && ~evenly
    error('Not implemented yet!');
elseif distinctBalls && distinctBoxes && ~nonemptyBoxes && evenly
    error('Not implemented yet!');
elseif distinctBalls && distinctBoxes && ~nonemptyBoxes && ~evenly
    error('Not implemented yet!');
elseif distinctBalls && ~distinctBoxes && nonemptyBoxes && evenly
    error('Not implemented yet!');
elseif distinctBalls && ~distinctBoxes && nonemptyBoxes && ~evenly
    error('Not implemented yet!');
elseif distinctBalls && ~distinctBoxes && ~nonemptyBoxes && evenly
    error('Not implemented yet!');
elseif distinctBalls && ~distinctBoxes && ~nonemptyBoxes && ~evenly
    error('Not implemented yet!');
elseif ~distinctBalls && distinctBoxes && nonemptyBoxes && evenly
    error('Not implemented yet!');
elseif ~distinctBalls && distinctBoxes && nonemptyBoxes && ~evenly
    error('Not implemented yet!');
elseif ~distinctBalls && distinctBoxes && ~nonemptyBoxes && evenly
    error('Not implemented yet!');
elseif ~distinctBalls && distinctBoxes && ~nonemptyBoxes && ~evenly
    error('Not implemented yet!');
elseif ~distinctBalls && ~distinctBoxes
    % Balls and boxes are all identical, 
    %   so we need to partition nBalls into nBoxes numbers
% TODO: Try a version starting from the evenly position

    % Start a counter
    ct = 0;

    % Initialize a partition
    if evenly
        % Divide nBalls by nBoxes and get the quotient and remainder
        quotient = floor(nBalls / nBoxes);
        remainder = mod(nBalls, nBoxes);

        % Distribute as evenly as possible
        partition = quotient * ones(1, nBoxes);
        partition(1:remainder) = quotient + 1;

        % There is no need for further computations
        partitionPaths = [];
    elseif nonemptyBoxes
        % First place one ball in each box
        partition = ones(1, nBoxes);

        % Then place the remaining balls in the first box
        partition(1) = nBalls - nBoxes + 1;
        partitionPaths = partition;
    else
        % Initialize a partition with all balls in the first box
        partition = zeros(1, nBoxes);
        partition(1) = nBalls;
        partitionPaths = partition;
    end

    % Place it in allWays output
    ct = ct + 1;
    allWays(ct, :) = partition;

    % Continue to move a ball from a box to one of its right neighbors 
    %   (starting from the left-most right neighbor), as long as the box has 
    %   at least two more balls than that right neighbor.
    %   If there is more than one move possible, split into two different paths.
    %   If a path has no more moves, remove it
    %   Stop when there are no more paths
    while ~isempty(partitionPaths)
        % Use the last partition path if possible
        partition = partitionPaths(end, :);

        % For each box, look for the left-most right neighbor 
        %   with at least two fewer balls
        indices = num2cell(1:(nBoxes - 1));
        idxReceiveAll = cellfun(@(x) x + find(partition((x + 1):end) <= ...
                                              partition(x) - 2, 1, 'first'), ...
                                indices, 'UniformOutput', false);

        % Determine whether each box has a ball to be moved
        hasBallToMove = ~isemptycell(idxReceiveAll);

        % Determine whether each box's immediate right neighbor has
        %   the same number of balls
        rightNeighborTheSame = ...
            cellfun(@(x) partition(x + 1) == partition(x), indices);

        % Find the indices of all boxes that has a ball to be moved
        %   but does not have an immediate right neighbor with the 
        %   same number of balls
        indToMove = find(hasBallToMove & ~rightNeighborTheSame);

        % Count the number of boxes that has a ball to be moved
        nToMove = length(indToMove);

        % For each possible move, perform it and add it as a partition path
        for iToMove = 1:nToMove
            % Copy the previous partition
            partitionThis = partition;

            % Get the index of the box to move
            idxToMove = indToMove(iToMove);

            % Get the index of the box to receive
            idxReceive = idxReceiveAll{idxToMove};

            % Move a ball from the box to move to the box to receive
            partitionThis(idxToMove) = partitionThis(idxToMove) - 1;
            partitionThis(idxReceive) = partitionThis(idxReceive) + 1;

            % Check if this partition already exists
            % TODO

            % Increment the counter and store the new partition in allWays
            ct = ct + 1;
            allWays(ct, :) = partitionThis;

            % Update the partition in partitionPaths
            partitionPaths(end + iToMove - 1, :) = partitionThis;
        end

        % If nToMove is zero, reset everything
        if nToMove == 0
            % Remove the partition from the partition paths
            partitionPaths(end, :) = [];
        end
    end

    % Collapse identical partitions
    allWays = unique(allWays, 'rows');

    % Record the number of ways
    nWays = size(allWays, 1);
else
    error('Problem with code logic!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:   

% For each right neighbor, check if it has at least two fewer balls
%   than the box in question
canReceive = partition(idxToMove + 1:end) <= partition(idxToMove) - 2;

% If a box can be received, move a ball
if any(canReceive)
    % Choose the left-most such right neighbor
    idxReceive = idxToMove + find(canReceive, 1, 'first');

    % Start with the second box from the right for a move-ball-attempt
    idxToMove = nBoxes - 1;

    % Move a ball from the box to move to the box to receive
    partition(idxToMove) = partition(idxToMove) - 1;
    partition(idxReceive) = partition(idxReceive) + 1;

    % Increment the counter and store the new partition in allWays
    ct = ct + 1;
    allWays{ct, 1} = partition;

    % Update the partition in partitionPaths
    partitionPaths(1, :) = partition;

    % Reset to the second box from the right
    idxToMove = nBoxes - 1;
else            
    % Decrement the index of the box for a move-ball-attempt
    idxToMove = idxToMove - 1;
end

%       allWays     - all ways of distribution
%                   specified as a column cell array of row numeric vectors

ct = ct + 1;
allWays{ct, 1} = partition;

% Collapse identical partitions
allWays = ...
    cellfun(@str2num, ...
            unique(cellfun(@num2str, allWays, ...
                           'UniformOutput', false), ...
                   'stable'), ...
            'UniformOutput', false);

% Record the number of ways
nWays = numel(allWays);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%