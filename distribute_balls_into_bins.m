function [allWays, nWays] = distribute_balls_into_boxes (nBalls, nBoxes, varargin)
%% Returns the ways and number of ways to distribute identical/discrete balls into identical/discrete boxes
% Usage: [allWays, nWays] = distribute_balls_into_boxes (nBalls, nBoxes, varargin)
%
% Arguments:    
%       nBalls      - number of balls to place in boxes
%                   must be a positive integer scalar
%       nBoxes       - number of boxes to accept balls
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
%
% Used by:
%       /home/Matlab/Adams_Functions/piecelinspace.m
%
% File History:
%   2018-08-02 - Created by Adam Lu
%   TODO: Implement the cases where balls and/or boxes are distinct

%% Hard-coded parameters

%% Default values for optional arguments
distinctBallsDefault = false;           % balls are identical by default
distinctBoxesDefault = false;           % boxes are identical by default
nonemptyBoxesDefault = false;           % boxes can be empty by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'nBalls', ...                  % number of balls
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addRequired(iP, 'nBoxes', ...                  % number of boxes
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'DistinctBalls', distinctBallsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'DistinctBoxes', distinctBoxesDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'NonemptyBoxes', nonemptyBoxesDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, nBalls, nBoxes, varargin{:});
distinctBalls = iP.Results.DistinctBalls;
distinctBoxes = iP.Results.DistinctBoxes;
nonemptyBoxes = iP.Results.NonemptyBoxes;

%% Distribute balls!
if distinctBalls && distinctBoxes && nonemptyBoxes
    error('Not implemented yet!');
elseif distinctBalls && distinctBoxes && ~nonemptyBoxes
    error('Not implemented yet!');
elseif distinctBalls && ~distinctBoxes && nonemptyBoxes
    error('Not implemented yet!');
elseif distinctBalls && ~distinctBoxes && ~nonemptyBoxes
    error('Not implemented yet!');
elseif ~distinctBalls && distinctBoxes && nonemptyBoxes
    error('Not implemented yet!');
elseif ~distinctBalls && distinctBoxes && ~nonemptyBoxes
    error('Not implemented yet!');
elseif ~distinctBalls && ~distinctBoxes && nonemptyBoxes
    error('Not implemented yet!');
elseif ~distinctBalls && ~distinctBoxes && ~nonemptyBoxes
    % Balls and boxes are all identical, 
    %   so we need to partition nBalls into nBoxes numbers

    % Compute the number 

else
    error('Problem with code logic!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:   

%}