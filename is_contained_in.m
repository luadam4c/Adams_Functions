function isSubset = is_contained_in (setSmall, setLarge, varargin)
%% Checks whether all elements of the first set are elements of the second set and print the ones that aren't
% Usage: isSubset = is_contained_in (setSmall, setLarge, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       isSubset    - whether all elements of setSmall are part of setLarge
%
% Arguments:
%       setSmall    - the first set (to be checked whether contained)
%                   must be recognized by the ismember() function
%       setLarge    - the second set (to be checked whether contains)
%                   must be recognized by the ismember() function
%       varargin    - 'SuppressOutput': whether to suppress standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Used by:
%       cd/update_param_values.m

% File History:
% 2018-11-14 Created by Adam Lu
% TODO: Rename as is_subset.m

%% Default values for optional arguments
suppressOutputDefault = false;          % whether to suppress standard output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SuppressOutput', suppressOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
suppressOutput = iP.Results.SuppressOutput;

%% Do the job
% Return whether each parameter name exist in the table
existsInSet = ismember(setSmall, setLarge);

% If any of the parameters does not exist, print the parameter(s)
if ~all(existsInSet)
    % Print message(s)
    if ~suppressOutput
        for iElement = find(~existsInSet)
            % Get the current member
            if iscell(setSmall)
                thisElement = setSmall{iElement};
            else
                thisElement = num2str(setSmall(iElement));
            end

            % Print message
            fprintf('%s does not exist in %s!\n', thisElement, inputname(2));
        end
    end

    % Return false
    isSubset = false;
else
    % Return true
    isSubset = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
