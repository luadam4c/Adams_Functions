function newVec = trim_nans (oldVec, varargin)
%% Removes leading and trailing NaNs from vector(s)
% Usage: newVec = trim_nans (oldVec, trimMethod, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       trim_nans([NaN NaN 3 2 NaN 3 1 NaN NaN])
%       trim_nans([NaN; NaN; 3; 2; NaN; 3; 1; NaN; NaN])
%       TODO: trim_nans([NaN NaN 3 2; NaN 3 1 NaN])
%
% Outputs:
%       newVec      - trimmed vector(s)
%                   specified as a TODO
%
% Arguments:
%       oldVec      - original vector(s)
%                   must be a TODO
%       trimMethod  - (opt) method for trimming
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'leading'  - trim leading NaNs
%                       'trailing' - trim trailing NaNs
%                       'flanks'   - trim leading and trailing NaNs
%                       'all'      - don't trim
%                       'none'     - remove all NaNs
%                   default == 'flanks'
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/parse_phase_info.m

% File History:
% 2019-08-22 Created by Adam Lu
% TODO: Expand to many vectors
% TODO: Option to trim rows
% TODO: Option to trim columns

%% Hard-coded parameters
validTrimMethods = {'leading', 'trailing', 'flanks', 'all', 'none'};

%% Default values for optional arguments
trimMethodDefault = 'flanks';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'oldVec');

% Add optional inputs to the Input Parser
addOptional(iP, 'trimMethod', trimMethodDefault, ...
    @(x) any(validatestring(x, validTrimMethods)));

% Read from the Input Parser
parse(iP, oldVec, varargin{:});
trimMethod = validatestring(iP.Results.trimMethod, validTrimMethods);

% Keep unmatched arguments for the TODO() function
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
% Test whether each element is NaN
isNaN = isnan(oldVec);

% Find the first index that's not NaN
idxFirstToLeave = find(~isNaN, 1, 'first');

% Find the last index that's not NaN
idxLastToLeave = find(~isNaN, 1, 'last');

% Extract new vector
switch trimMethod
    case 'leading'
        newVec = oldVec(idxFirstToLeave:end);
    case 'trailing'
        newVec = oldVec(1:idxLastToLeave);
    case 'flanks'
        newVec = oldVec(idxFirstToLeave:idxLastToLeave);
    case 'all'
        newVec = oldVec(~isNaN);
    case 'none'
        newVec = oldVec;
    otherwise
        error('trimMethod unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%