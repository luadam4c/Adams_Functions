function nearestOdd = find_nearest_odd (realNumber, varargin)
%% Returns the nearest odd integer to real number(s)
% Usage: nearestOdd = find_nearest_odd (realNumber, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       nearestOdd  - the nearest odd integer returned
%                   specified as an odd integer vector
% Arguments:    
%       realNumber  - real number(s)
%                   must be a numeric vector
%       varargin    - 'Direction': rounding direction
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'nearest'   - round to 'nearest'
%                       'down'      - always round down
%                       'up'        - always round up
%                   default == 'round'
%
% Used by:
%       cd/medianfilter.m
%       cd/movingaveragefilter.m

% File History:
% 2018-10-12 Created by Adam Lu & ShinShin Nien
% 2018-12-18 Now accepts a vector as an argument
% 

%% Hard-coded parameters
validDirections = {'nearest', 'down', 'up'};

%% Default values for optional arguments
directionDefault = 'nearest';

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
addRequired(iP, 'realNumber', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Direction', directionDefault, ...
    @(x) any(validatestring(x, validDirections)));

% Read from the Input Parser
parse(iP, realNumber, varargin{:});
direction = validatestring(iP.Results.Direction, validDirections);

%% Do the job
% Get the sign(s) of the real number(s)
sgn = sign(realNumber);

% Take the absolute value
absValue = abs(realNumber);

% Find the nearest positive odd number
absOdd = arrayfun(@(x) find_nearest_positive_odd(x, direction), absValue);

% Restore the sign
nearestOdd = sgn * absOdd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function absOdd = find_nearest_positive_odd (absValue, direction)
%% Returns the nearest positive odd number

if absValue < 1
    % If less than 1, use 1
    absOdd = 1;
else
    % Otherwise, transform -> round -> transform back
    switch direction
    case 'nearest'
        absOdd = round(( absValue - 1 ) / 2 ) * 2 + 1;
    case 'down'
        absOdd = floor(( absValue - 1 ) / 2 ) * 2 + 1;
    case 'up'
        absOdd = ceil(( absValue - 1 ) / 2 ) * 2 + 1;
    otherwise
        error('direction unrecognized!!');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
