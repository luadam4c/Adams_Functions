function [result, distToTarget] = find_nearest_multiple (base, target, varargin)
%% Finds the nearest integer multiple of base to target and the distance to between it and target
% Usage: [result, distToTarget] = find_nearest_multiple (base, target, varargin)
% Explanation:
%       TODO
% Example(s):
%       [res, dist] = find_nearest_multiple(3, 5)
%       [res, dist] = find_nearest_multiple(4, 10)
%       [res, dist] = find_nearest_multiple([1, 2], [4, 4])
%       [res, dist] = find_nearest_multiple([1, 2], [5, 4])
% Outputs:
%       result      - nearest multiple of base to target
%                   specified as a numeric vector
%       distToTarget- distance between the nearest multiple and target
%                   specified as a numeric scalar
% Arguments:
%       base        - base for the multiple
%                   must be a numeric vector
%       target      - target to find nearest to
%                   must be a numeric vector
%       varargin    - 'RelativeToHalfBase': whether to specify distance as
%                                       relative to half the length of base
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/parse_multiunit.m

% File History:
% 2019-03-15 Created by Adam Lu
% 2019-03-17 Now uses dot()

%% Hard-coded parameters

%% Default values for optional arguments
relativeToHalfBaseDefault = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'base', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'target', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'RelativeToHalfBase', relativeToHalfBaseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, base, target, varargin{:});
relativeToHalfBase = iP.Results.RelativeToHalfBase;

% Check relationships between arguments
if numel(base) ~= numel(target)
    error('base and target must have equal number of elements!');
end

%% Preparation
% Find the magnitude of the base vector
baseNorm = norm(base);

%% Do the job
if isscalar(base)
    % Round to the nearest integer after division
    nearestMultiple = round(target / base);

    % Retrieve the nearest multiple
    result = nearestMultiple * base;

    % Compute the distance between result and target 
    distToTarget = norm(result - target);
else
    % Find the projection of target onto the base vector
    projection = dot(target, base) / baseNorm;

    % Find the two possible candidates
    cand1 = ceil(projection / baseNorm) * base;
    cand2 = floor(projection / baseNorm) * base;

    % Compute the distances to target
    distToTarget1 = norm(cand1 - target);
    distToTarget2 = norm(cand2 - target);

    % Compute the 
    if distToTarget1 < distToTarget2
        result = cand1;
        distToTarget = distToTarget1;
    else
        result = cand2;
        distToTarget = distToTarget2;
    end
end

if relativeToHalfBase
    % Compute half magnitude of base
    halfBaseNorm = baseNorm / 2;

    % Compute the distance relative to halfBaseNorm
    distToTarget = distToTarget / halfBaseNorm;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

projection = sum(target .* base) / baseNorm;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%