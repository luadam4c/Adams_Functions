function slope = compute_slope(x, y, idxStart, idxEnd)
%% Computes the slope given two vectors and two indices
% Usage: slope = compute_slope(x, y, idxStart, idxEnd)
% Outputs:
%       slope       - dy/dx
%                   specified as a numeric scalar
% Arguments:    
%       x           - the x vector
%                   must be a numeric vector
%       y           - the y vector
%                   must be a numeric vector
%       idxStart    - the starting index
%                   must be a numeric positive scalar
%       idxEnd      - the ending index
%                   must be a numeric positive scalar
%
% Used by:

% File History:
% 2018-06-11 Created by Adam Lu
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 4
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'x', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'y', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'idxStart', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addRequired(iP, 'idxEnd', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Read from the Input Parser
parse(iP, x, y, idxStart, idxEnd);

%% Perform job
% Compute the slope
slope = (y(idxEnd) - y(idxStart)) / (x(idxEnd) - x(idxStart));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%