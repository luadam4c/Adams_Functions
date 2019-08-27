function examine_figure_objects (varargin)
%% Examines figure objects
% Usage: examine_figure_objects (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
%
% Arguments:
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%       /TODO:dir/TODO:file
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-08-27 Adapted from https://www.mathworks.com/matlabcentral/answers/5814-plot-move-data-on-top
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the TODO() function
% otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
set(gcf,'windowbuttonmotionfcn', 'uistack(hittest,''top'')');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%