function varargout = apply_to_all_cells (myFunction, inputs, varargin)
%% Applies a function to inputs, separately to each cell if a inputs is a cell array
% Usage: varargout = apply_to_all_cells (myFunction, inputs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       apply_to_all_cells(@(x) x .* 3, {[1, 2], {3, 4}})
%
% Outputs:
%       outputs     - outputs
%                   specified as the same type as inputs
%
% Arguments:
%       myFunction  - a custom function that takes one required argument
%                   must be a function handle
%       inputs      - inputs
%                   must be an array
%       varargin    - 'OptArg': optional argument that's 
%                                   not a parameter-value pair 
%                   default == []
%                   - optional parameter-value pairs for myFunction
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/compare_events_pre_post_stim.m
%       cd/convert_units.m
%       cd/m3ha_find_decision_point.m

% File History:
% 2020-01-22 Created by Adam Lu
% 2020-06-29 Now uses varargout

%% Hard-coded parameters

%% Default values for optional arguments
optArgDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'myFunction', ...           % a custom function
    @(x) validateattributes(x, {'function_handle'}, {'scalar'}));
addRequired(iP, 'inputs');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OptArg', optArgDefault);

% Read from the Input Parser
parse(iP, myFunction, inputs, varargin{:});
optArg = iP.Results.OptArg;

% Keep unmatched arguments for myFunction
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
if iscell(inputs)
    [varargout{1:nargout}] = ...
        cellfun(@(x) apply_to_all_cells(myFunction, x, ...
                                    'OptArg', optArg, otherArguments{:}), ...
                        inputs, 'UniformOutput', false);
else
    if ~isempty(optArg)
        [varargout{1:nargout}] = myFunction(inputs, optArg, otherArguments{:});
    else
        [varargout{1:nargout}] = myFunction(inputs, otherArguments{:});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%