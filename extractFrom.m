function strOut = extractFrom (strIn, startStr, varargin)
%% Extract the part of string(s) from a starting substring
% Usage: strOut = extractFrom (strIn, startStr, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       extractFrom('apple_orange_cat', 'or')
%
% Outputs:
%       strOut     - TODO: Description of strOut
%                   specified as a TODO
%
% Arguments:
%       strIn       - TODO: Description of strIn
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/m3ha_plot_figure07.m

% File History:
% 2020-04-19 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'strIn');
addRequired(iP, 'startStr');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, strIn, startStr, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
strOut = strcat(startStr, extractAfter(strIn, startStr));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%