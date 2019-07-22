function flag = set_default_flag(flag, varargin)
%% Sets the default flag if empty according to an optional auxFlag
% Usage: flag = set_default_flag(flag, auxFlag (opt))
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       flag        - flag that has been set
%                   specified as a logical scalar
% Arguments:
%       flag        - flag that will be set if empty
%                   must be either empty or a logical scalar
%       auxFlag     - auxiliary flag that sets the default flag
%                   must be numeric/logical 1 (true) or 0 (false)
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/parse_multiunit.m

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

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'auxFlag', auxFlagDefault, ...
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