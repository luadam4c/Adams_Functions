function flag = set_default_flag(flag, varargin)
%% Sets the default flag if empty according to an optional auxFlag
% Usage: flag = set_default_flag(flag, auxFlag (opt))
% Explanation:
%       If flag is already set, nothing is done
%       If flag is empty, it is set to be true if auxFlag is true,
%           but false if auxFlag is false.
%       If auxFlag is not provided, it is false.
% Example(s):
%       set_default_flag([], 0)
%       set_default_flag([], true)
%       set_default_flag([], [])
% Outputs:
%       flag        - flag that has been set
%                   specified as a logical scalar
% Arguments:
%       flag        - flag that will be set if empty
%                   must be either empty or a logical scalar
%       auxFlag     - (opt) auxiliary flag that sets the default flag
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/parse_multiunit.m
%       cd/plot_measures.m

% File History:
% 2019-07-22 Pulled from parse_multiunit.m
% 2019-07-22 Made auxFlag an optional argument
% 

%% Hard-coded parameters

%% Default values for optional arguments
auxFlagDefault = false;

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
addRequired(iP, 'flag');

% Add optional inputs to the Input Parser
addOptional(iP, 'auxFlag', auxFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, flag, varargin{:});
auxFlag = iP.Results.auxFlag;

%% Do the job
if isempty(flag)
    if auxFlag
        flag = true;
    else
        flag = false;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%