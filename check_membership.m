function isAllMembers = check_membership (members, setToCheck, varargin)
%% Checks whether all elements of the first set are elements of the second set and print the ones that aren't
% Usage: isAllMembers = check_membership (members, setToCheck, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       isAllMembers    - whether all elements of members are part of set
% Arguments:
%       members     - values that should be a member of the set
%                   must be recognized by the ismember() function
%       setToCheck  - set of data values
%                   must be recognized by the ismember() function
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
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
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
existsInSet = ismember(members, setToCheck);

% If any of the parameters does not exist, print the parameter(s)
if ~all(existsInSet)
    % Print message(s)
    if ~suppressOutput
        for iElement = find(~existsInSet)
            % Get the current member
            if iscell(members)
                thisElement = members{iElement};
            else
                thisElement = num2str(members(iElement));
            end

            % Print message
            fprintf('%s does not exist in %s!\n', thisElement, inputname(2));
        end
    end

    % Return false
    isAllMembers = false;
else
    % Return true
    isAllMembers = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%