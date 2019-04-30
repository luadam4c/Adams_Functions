function [y, ia, ic] = unique_custom (x, varargin)
%% Returns the unique values in x, optionally without NaN
% Usage: [y] = unique_custom (x, varargin)
% Explanation:
%       This function removes NaN values not identified by the default unique() function.
% Example(s):
%       unique_custom([3, NaN, 3, 5, NaN], 'IgnoreNaN', true)
%       [y, ia, ic] = unique_custom([3 NaN 3 5 NaN])
%               = unique([3 NaN 3 5 NaN])
%               = unique_custom([3 NaN 3 5 NaN], 'IgnoreNaN', false, 'TreatNanAsEqual', true)
%               = unique_custom([3 NaN 3 5 NaN], 'IgnoreNaN', false, 'TreatNanAsEqual', false)
%               = unique_custom([3 NaN 3 5 NaN], 'TreatNanAsEqual', true)
%       [y, ia, ic] = unique_custom([3 NaN 3 5 NaN], 'IgnoreNaN', true)
%                   = unique_custom([3 NaN 3 5 NaN], 'IgnoreNaN', true, 'TreatNanAsEqual', true)
%       [y, ia, ic] = unique_custom([3 NaN 3 5 NaN], 'IgnoreNaN', true, 'TreatNanAsEqual', false)
%       
% Outputs:
%       y           - All unique values in x
%                   specified as a array
% Arguments:
%       x           - Matrix to check unique values
%                   must be a array
%       varargin    - 'IgnoreNaN': whether to include NaN as distinct elements
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatNanAsEqual': whether to treat all NaN values
%                                           as the same
%                       Note: If 'IgnoreNaN' == false, 
%                           'TreatNanAsEqual' has no effect
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for the unique() function
%
% Requires:
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-04-01 BT - Adapted from https://www.mathworks.com/matlabcentral/
%                         answers/42561-treating-nan-as-a-
%                         unique-value-instead-of-as-a-distinct#answer_52371
% 

%% Hard-coded parameters

%% Default values for optional arguments
ignoreNanDefault = false;  	    % do not ignore NaN by default
treatNanAsEqualDefault = false; 	% do not treat all NaN values equal by default

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
addRequired(iP, 'x', ...                  % array to be operated on
    @(z) validateattributes(z, {'char', 'string', 'cell', 'numeric'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'IgnoreNaN', ignoreNanDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'scalar'}));
addParameter(iP, 'TreatNanAsEqual', treatNanAsEqualDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'scalar'}));

% Read from the Input Parser
parse(iP, x, varargin{:});
ignoreNan = iP.Results.IgnoreNaN;
treatNanAsEqual = iP.Results.TreatNanAsEqual;

% Keep unmatched arguments for the unique_custom() function
otherArguments = struct2arglist(iP.Unmatched);

% Check relationships between arguments

%% Preparation
% Initial unique matrix y of x
[y, ia, ic] = unique(x, otherArguments{:});

%% Do the job

% Ignoring NaN
if ignoreNan
    % All NaN indices
    indNaN = isnan(y(1:end));

    % Delete all NaN elements in y, ia, and ic
    y(indNaN) = [];
    ia(indNaN) = [];
    ic(indNaN) = [];

    return;
end

% Preserving one NaN
if treatNanAsEqual
    % NaN indices, does not include last NaN if present
    indNaN = isnan(y(1:end-1));

    % Delete all NaN elements in y, ia, and ic
    y(indNaN) = [];
    ia(indNaN) = [];
    ic(indNaN) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%