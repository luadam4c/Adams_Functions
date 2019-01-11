function varargout = find_custom (X, varargin)
%% Same as find() but takes custom parameter-value pairs
% Usage: varargout = find_custom (X, varargin)
% Explanation:
%       TODO
% Outputs:
%       varargout   - TODO: Description of k
%                   specified as a TODO
% TODO: Implement varargout
% Arguments:    
%       X           - Input array
%                   must be a scalar, vector, matrix, or multidimensional array
%       n           - (opt) Number of nonzero elements to find
%                   must be a positive integer scalar
%                   default == []
%       direction   - (opt) Search direction
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'first' - Search forward from the beginning
%                       'last'  - Search backward from the end
%                   default == 'first'
%       varargin    - 'ReturnNan': Return NaN instead of empty if nothing found
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/compute_peak_decay.m
%       cd/compute_peak_halfwidth.m
%       cd/find_in_strings.m
%       /home/Matlab/Kojis_Functions/find_directional_events.m
%       /home/Matlab/minEASE/gui_examine_events.m
%
% File History:
% 2017-06-01 Created by Adam Lu

%% Default values for optional arguments
nDefault  = [];             % default number of nonzero elements to find
directionDefault = 'first'; % default search direction
returnNanDefault = false;   % whether to return NaN instead of empty 
                            %   if nothing found by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add required inputs to an Input Parser
addRequired(iP, 'X', ...                        % Input array
    @(x) validateattributes(x, {'numeric', 'logical', 'char'}, {}));

% Add optional inputs to the Input Parser
addOptional(iP, 'n', nDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addOptional(iP, 'direction', directionDefault, ...
    @(x) any(validatestring(x, {'first', 'last'})));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ReturnNan', returnNanDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, X, varargin{:});
n = iP.Results.n;
direction = validatestring(iP.Results.direction, {'first', 'last'});
returnNan = iP.Results.ReturnNan;

%% Apply find() based on nargout
if nargout == 1

    % Apply find() based on input
    if ~isempty(n)
        k = find(X, n, direction);
    else
        k = find(X);
    end

    % If returnNan is true, change empty matrices into NaN
    if returnNan && isempty(k)
        k = NaN;
    end

    % Return output
    varargout{1} = k;

elseif nargout >= 2

    % Apply find() based on input
    if ~isempty(n)
        [row, col, v] = find(X, n, direction);
    else
        [row, col, v] = find(X);
    end

    % Return output
    varargout{1} = row;
    varargout{2} = col;
    varargout{3} = v;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
