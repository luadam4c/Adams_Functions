function isAllWithinBounds = ...
                check_within_bounds (values, lowerBounds, upperBounds, varargin)
%% Checks whether all values are within bounds and print the ones that aren't
% Usage: isAllWithinBounds = ...
%               check_within_bounds (values, lowerBounds, upperBounds, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       isAllWithinBounds     - whether all values are within bounds
% Arguments:
%       values      - values to check
%                   must be a numeric vector
%       lowerBounds - lower bounds
%                   must be a numeric vector
%       upperBounds - upper bounds
%                   must be a numeric vector
%       varargin    - 'SuppressOutput': whether to suppress standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Used by:
%       cd/update_param_values.m

% File History:
% 2018-11-14 Created by Adam Lu

%% Default values for optional arguments
suppressOutputDefault = false;          % whether to suppress standard output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'values', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'lowerBounds', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'upperBounds', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SuppressOutput', suppressOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, values, lowerBounds, upperBounds, varargin{:});
suppressOutput = iP.Results.SuppressOutput;

%% Do the job
% Return whether each parameter name exist in the table
isWithinBound = values <= upperBounds & values >= lowerBounds;

% If any of the parameters does not exist, print the parameter(s)
if ~all(isWithinBound)
    % Print message(s)
    if ~suppressOutput
        for iElement = find(~isWithinBound)
            % Get the current member
            thisElement = values(iElement);
            thisLowerBound = lowerBounds(iElement);
            thisUpperBound = upperBounds(iElement);

            % Print message
            fprintf('%g is not within the range [%g, %g]!\n', ...
                    thisElement, thisLowerBound, thisUpperBound);
        end
    end

    % Return false
    isAllWithinBounds = false;
else
    % Return true
    isAllWithinBounds = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%