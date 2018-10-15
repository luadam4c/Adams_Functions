function nearestOdd = nearest_odd (realNumber, varargin)
%% Returns the nearest odd integer of a real number
% Usage: nearestOdd = nearest_odd (realNumber, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       nearestOdd  - the nearest odd integer returned
%                   specified as an odd integer
% Arguments:    
%       realNumber  - a real number
%                   must be a number
%       varargin    - 'Direction': rounding direction
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'nearest' - round to 'nearest'
%                       'down' - always round down
%                       'up'  - always round up
%                   default == 'round'
%
% Used by:    
%       cd/find_passive_params.m

% File History:
% 2018-10-12 Created by ShinShin Nien
% 

%% Hard-coded parameters
validDirections = {'nearest', 'down', 'up'};

%% Default values for optional arguments
directionDefault = 'nearest';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'realNumber', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Direction', directionDefault, ...
    @(x) any(validatestring(x, validDirections)));

% Read from the Input Parser
parse(iP, realNumber, varargin{:});
direction = validatestring(iP.Results.Direction, validDirections);

%% Do the job
% Get the sign of the real number
sgn = sign(realNumber);

% Take the absolute value
absValue = abs(realNumber);

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

% Restore the sign
nearestOdd = sgn * absOdd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%