function varargout = apply_or_return (toApply, myFunction, varargin)
%% Apply a function if a condition is true, or return the original argument(s)
% Usage: varargout = apply_or_return (toApply, myFunction, varargin)
% Explanation:
%       This function applies a function to arguments only
%           if the first argument is true or 1.
%       Otherwise, the function simply returns the input arguments
% Example(s):
%       a = apply_or_return(true, @sum, 1:10)
%       b = apply_or_return(false, @sum, 1:10)
%       c = cellfun(@(x) apply_or_return(~isempty(x), @mean, x), {[], 1:5}, ...
%               'UniformOutput', false)
% Outputs:
%       varargout   - whatever the function outputs
% Arguments:
%       toApply     - whether to apply the function
%                   must be a function handle
%       myFunction  - a custom function
%                   must be a function handle
%       varargin    - input arguments
%
% Used by:    
%       cd/match_format_vectors.m
%       cd/m3ha_import_raw_traces.m

% File History:
% 2018-10-31 Created by Adam Lu
% 

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

% Add required inputs to the Input Parser
addRequired(iP, 'toApply', ...              % whether to apply the function
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addRequired(iP, 'myFunction', ...           % a custom function
    @(x) validateattributes(x, {'function_handle'}, {'scalar'}));

% Read from the Input Parser
parse(iP, toApply, myFunction);

%% Do the job
if toApply
    % Get the number of output arguments that will be returned by the function
    nOutputsPossible = nargout(myFunction);

    % Check the number of outputs requested
    nargoutchk(0, nOutputsPossible);

    % Get the number of outputs requested
    nOutputsRequested = nargout;

    % Initialize a cell array with the requested number of outputs
    varargout = cell(nOutputsRequested, 1);

    % Place the outputs in the cell array
    [varargout{:}] = myFunction(varargin{:});
else
    varargout = varargin;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%