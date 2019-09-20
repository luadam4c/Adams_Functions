function figHandle = update_figure_for_corel (figHandle, varargin)
%% Update figure to be journal-friendly
% Usage: figHandle = update_figure_for_corel (figHandle, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       figHandle     - TODO: Description of figHandle
%                   specified as a TODO
%
% Arguments:
%       figHandle     - TODO: Description of figHandle
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       ~/Adams_Functions/create_error_for_nargin.m
%       ~/Adams_Functions/struct2arglist.m
%       /TODO:dir/TODO:file
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-19 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'figHandle');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, figHandle, varargin{:});
param1 = iP.Results.param1;

% Keep unmatched arguments for the TODO() function
otherArguments = struct2arglist(iP.Unmatched);

% Check relationships between arguments
% TODO

%% Preparation
% TODO

%% Do the job
% Find all axes in the figure
ax = findall(figHandle, 'type', 'axes');

% Count the number of axes
nAx = numel(ax);

% Make all ticks go outward
set(ax, 'TickDir', 'out');
set(ax, 'TickDirMode', 'manual');

% Make all rulers linewidth 1
set(ax, 'LineWidth', 1);

% Remove boxes
for iAx = 1:nAx
    box(ax(iAx), 'off');
end

% Remove x axis ticks
set(ax, 'XTick', []);

% Remove y axis ticks
set(ax, 'YTick', []);


%% Output results
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

for iAx = 1:nAx
    ax(iAx).XAxis.LineWidth = 1;
    ax(iAx).YAxis.LineWidth = 1;
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%